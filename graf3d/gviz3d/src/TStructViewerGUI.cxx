// @(#)root/gviz3d:$Id$
// Author: Tomasz Sosnicki   18/09/09

/************************************************************************
* Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

#include "TStructViewerGUI.h"
#include "TStructViewer.h"
#include "TStructNodeEditor.h"
#include "TStructNodeProperty.h"
#include "TStructNode.h"
#include <TCanvas.h>
#include <RQ_OBJECT.h>
#include <TGLLogicalShape.h>
#include <TGLPhysicalShape.h>
#include <TGLWidget.h>
#include <TGButtonGroup.h>
#include <TGSplitter.h>
#include <TList.h>
#include <TClass.h>
#include <TDataMember.h>
#include <TExMap.h>
#include <TPolyLine3D.h>
#include <TGTab.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TROOT.h>
#include <TApplication.h>

ClassImp(TStructViewerGUI);

//________________________________________________________________________
//////////////////////////////////////////////////////////////////////////
//
// TStructViewerGUI is main window of TStructViewer. It provides graphical
// interface. In the window we can find panel with tabs and frame with
// GLViewer. Tab "Info" serves information about node and is used to naviagate
// backward and forward. Second tab "Options" is used to set few options
// such as links visibility, scaling method or setting a pointer.
// Last tab "Editor" is tab when the TStructNodeEditor is placed.
//
//////////////////////////////////////////////////////////////////////////

TGeoMedium* TStructViewerGUI::fgMedium = NULL;
UInt_t      TStructViewerGUI::fgCounter = 0;

////////////////////////////////////////////////////////////////////////////////
/// Constructs window with "w" as width, "h" as height and given parent "p". Argument "parent" is a pointer to TStructViewer which contains this GUI.
/// This constructor build window with all controls, build map with colors, init OpenGL Viewer and create TGeoVolumes.

TStructViewerGUI::TStructViewerGUI(TStructViewer* parent, TStructNode* nodePtr, TList* colors, const TGWindow *p,UInt_t w,UInt_t h)
   : TGMainFrame(p, w, h, kHorizontalFrame)
{
   fParent = parent;
   fNodePtr = nodePtr;

   fMaxSlices = 10;
   fMouseX = 0;
   fMouseY = 0;
   fSelectedObject = NULL;
   fMaxRatio = 0;
   fColors = colors;

   if (!gGeoManager) new TGeoManager("tmp","tmp");
   if (!fgMedium) {
      fgMedium = new TGeoMedium("MED",1,new TGeoMaterial("Mat", 26.98,13,2.7));
   }

   SetCleanup(kDeepCleanup);
   //////////////////////////////////////////////////////////////////////////
   // layout
   //////////////////////////////////////////////////////////////////////////
   TGVerticalFrame* leftFrame = new TGVerticalFrame(this, 200, 200, kFixedWidth);
   this->AddFrame(leftFrame, new TGLayoutHints(kFixedWidth, 1, 1, 1, 1));
   TGTab* tabs = new TGTab(leftFrame);
   TGLayoutHints* expandX = new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 5,5,5,5);
   //////////////////////////////////////////////////////////////////////////
   // INFO
   //////////////////////////////////////////////////////////////////////////
   TGCompositeFrame* infoFrame = tabs->AddTab("Info");
   TGGroupFrame* fInfoMenu = new TGGroupFrame(infoFrame, "Info");
   fNodeNameLabel = new TGLabel(fInfoMenu, "Name:");
   fInfoMenu->AddFrame(fNodeNameLabel, expandX);
   fNodeTypelabel = new TGLabel(fInfoMenu, "Type:");
   fInfoMenu->AddFrame(fNodeTypelabel, expandX);
   fMembersCountLabel = new TGLabel(fInfoMenu, "Members:");
   fInfoMenu->AddFrame(fMembersCountLabel, expandX);
   fAllMembersCountLabel = new TGLabel(fInfoMenu, "All members:");
   fInfoMenu->AddFrame(fAllMembersCountLabel, expandX);
   fLevelLabel = new TGLabel(fInfoMenu, "Level:");
   fInfoMenu->AddFrame(fLevelLabel, expandX);
   fSizeLabel = new TGLabel(fInfoMenu, "Size:");
   fInfoMenu->AddFrame(fSizeLabel, expandX);
   fTotalSizeLabel = new TGLabel(fInfoMenu, "Total size:");
   fInfoMenu->AddFrame(fTotalSizeLabel, expandX);
   infoFrame->AddFrame(fInfoMenu, expandX);

   //////////////////////////////////////////////////////////////////////////
   // OPTIONS
   //////////////////////////////////////////////////////////////////////////
   TGCompositeFrame* options = tabs->AddTab("Options");

   fShowLinksCheckButton = new TGCheckButton(options, "Show links");
   fShowLinksCheckButton->Connect("Toggled(Bool_t)", "TStructViewerGUI", this, "ShowLinksToggled(Bool_t)");
   options->AddFrame(fShowLinksCheckButton);
   fShowLinksCheckButton->SetOn();

   TGVButtonGroup* scaleByGroup = new TGVButtonGroup(options, "Scale by");
   fScaleBySizeButton = new TGRadioButton(scaleByGroup, "Size");
   fScaleBySizeButton->Connect("Clicked()", "TStructViewerGUI", this, "ScaleByChangedSlot()");
   fScaleBySizeButton->SetOn();
   fScaleByMembersButton = new TGRadioButton(scaleByGroup, "Members count");
   fScaleByMembersButton->Connect("Clicked()", "TStructViewerGUI", this, "ScaleByChangedSlot()");
   options->AddFrame(scaleByGroup, expandX);

   TGHorizontalFrame* defaultColorFrame = new TGHorizontalFrame(options);
   options->AddFrame(defaultColorFrame, expandX);
   TGLabel* defColorlabel = new TGLabel(defaultColorFrame, "Default color");
   defaultColorFrame->AddFrame(defColorlabel, expandX);
   TGColorSelect* defColorSelect = new TGColorSelect(defaultColorFrame, GetDefaultColor()->GetPixel());
   defColorSelect->Connect("ColorSelected(Pixel_t)", "TStructViewerGUI", this, "ColorSelectedSlot(Pixel_t)");
   defaultColorFrame->AddFrame(defColorSelect);

   TGHorizontalFrame* boxHeightFrame = new TGHorizontalFrame(options);
   options->AddFrame(boxHeightFrame, expandX);
   TGLabel* boxHeightLabel = new TGLabel(boxHeightFrame, "Box height:");
   boxHeightFrame->AddFrame(boxHeightLabel, expandX);
   fBoxHeightEntry = new TGNumberEntry(boxHeightFrame, 0.1);
   fBoxHeightEntry->SetLimits(TGNumberEntry::kNELLimitMin, 0.01);
   fBoxHeightEntry->Connect("ValueSet(Long_t)", "TStructViewerGUI", this, "BoxHeightValueSetSlot(Long_t)");
   boxHeightFrame->AddFrame(fBoxHeightEntry);

   TGHorizontalFrame* levelDistanceFrame = new TGHorizontalFrame(options);
   options->AddFrame(levelDistanceFrame, expandX);
   TGLabel* lvlDistLabel = new TGLabel(levelDistanceFrame, "Distance between levels");
   levelDistanceFrame->AddFrame(lvlDistLabel, expandX);
   fLevelDistanceEntry = new TGNumberEntry(levelDistanceFrame, 1.1);
   fLevelDistanceEntry->SetLimits(TGNumberEntry::kNELLimitMin, 0.01);
   fLevelDistanceEntry->Connect("ValueSet(Long_t)", "TStructViewerGUI", this, "LevelDistValueSetSlot(Long_t)");
   levelDistanceFrame->AddFrame(fLevelDistanceEntry);

   fAutoRefesh = new TGCheckButton(options, "Auto refresh");
   fAutoRefesh->SetOn();
   fAutoRefesh->Connect("Toggled(Bool_t)", "TStructViewerGUI", this, "AutoRefreshButtonSlot(Bool_t)");
   options->AddFrame(fAutoRefesh, expandX);

   TGLabel* pointerLabel = new TGLabel(options, "Pointer:");
   options->AddFrame(pointerLabel, expandX);
   fPointerTextEntry = new TGTextEntry(options, "0x0000000");
   options->AddFrame(fPointerTextEntry, expandX);
   TGLabel* fPointerTypeLabel = new TGLabel(options, "Pointer Type:");
   options->AddFrame(fPointerTypeLabel, expandX);
   fPointerTypeTextEntry = new TGTextEntry(options, "TObject");
   options->AddFrame(fPointerTypeTextEntry, expandX);
   TGTextButton* setPointerButton = new TGTextButton(options, "Set pointer");
   setPointerButton->Connect("Clicked()", "TStructViewerGUI", this, "SetPointerButtonSlot()");
   options->AddFrame(setPointerButton, expandX);

   //////////////////////////////////////////////////////////////////////////
   // EDITOR
   //////////////////////////////////////////////////////////////////////////
   TGCompositeFrame* editTab = tabs->AddTab("Editor");
   fEditor = new TStructNodeEditor(fColors, editTab);
   fEditor->Connect("Update(Bool_t)", "TStructViewerGUI", this, "Update(Bool_t)");
   editTab->AddFrame(fEditor, expandX);

   leftFrame->AddFrame(tabs, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 1,1,1,1));

   TGVSplitter* splitter = new TGVSplitter(this);
   splitter->SetFrame(leftFrame, true);
   this->AddFrame(splitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   //////////////////////////////////////////////////////////////////////////
   // NAVIGATE
   //////////////////////////////////////////////////////////////////////////
   fUndoButton = new TGTextButton(leftFrame, "Undo");
   fUndoButton->Connect("Clicked()", "TStructViewerGUI", this, "UndoButtonSlot()");
   fUndoButton->SetEnabled(false);
   leftFrame->AddFrame(fUndoButton, expandX);

   fRedoButton = new TGTextButton(leftFrame, "Redo");
   fRedoButton->Connect("Clicked()", "TStructViewerGUI", this, "RedoButtonSlot()");
   fRedoButton->SetEnabled(false);
   leftFrame->AddFrame(fRedoButton, expandX);

   TGTextButton* resetCameraButton = new TGTextButton(leftFrame, "Reset camera");
   leftFrame->AddFrame(resetCameraButton, expandX);
   resetCameraButton->Connect("Clicked()", "TStructViewerGUI", this, "ResetButtonSlot()");

   TGTextButton* updateButton = new TGTextButton(leftFrame, "Update");
   updateButton->Connect("Clicked()", "TStructViewerGUI", this, "UpdateButtonSlot()");
   leftFrame->AddFrame(updateButton, expandX);

   TGTextButton* quitButton = new TGTextButton(leftFrame, "Quit");
   leftFrame->AddFrame(quitButton, expandX);
   quitButton->Connect("Clicked()", "TApplication", gApplication, "Terminate()");

   fTopVolume = gGeoManager->MakeBox("TOPVolume", fgMedium,100, 100, 100);
   gGeoManager->SetTopVolume(fTopVolume);
   gGeoManager->SetNsegments(40);

   fCanvas = new TCanvas("c", "c", 0, 0);
   // drawing after creating canvas to avoid drawing in default canvas
   fGLViewer = new TGLEmbeddedViewer(this, fCanvas);
   AddFrame(fGLViewer->GetFrame(), new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,10));
   fGLViewer->PadPaint(fCanvas);
   fGLViewer->Connect("MouseOver(TGLPhysicalShape*)", "TStructViewerGUI", this, "MouseOverSlot(TGLPhysicalShape*)");
   fGLViewer->GetGLWidget()->Connect("ProcessedEvent(Event_t*)", "TStructViewerGUI", this, "GLWidgetProcessedEventSlot(Event_t*))");
   fGLViewer->Connect("DoubleClicked()", "TStructViewerGUI", this, "DoubleClickedSlot()");
   fGLViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
   Update();
   fGLViewer->SetResetCamerasOnUpdate(false);

   SetWindowName("Struct Viewer");
   MapSubwindows();
   this->SetWMSizeHints(w, h, 2000, 2000, 0, 0);
   Resize(GetDefaultSize());
   MapWindow();

   fToolTip = new TGToolTip(0, 0, "ToolTip", 500);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TStructViewerGUI::~TStructViewerGUI()
{
   delete fCanvas;
}

////////////////////////////////////////////////////////////////////////////////
/// Activated when user chage condition

void TStructViewerGUI::AutoRefreshButtonSlot(Bool_t on)
{
   if (on) {
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emmited when user changes height of boxes

void TStructViewerGUI::BoxHeightValueSetSlot(Long_t /* h */)
{
   if(fAutoRefesh->IsOn()) {
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recursive method to calculating nodes posistion in 3D space

void TStructViewerGUI::CalculatePosistion(TStructNode* parent)
{
   // choose scaling method
   if (fScaleBySizeButton->GetState() == kButtonDown) {
      TStructNode::SetScaleBy(kSize);
   } else if (fScaleByMembersButton->GetState() == kButtonDown) {
      TStructNode::SetScaleBy(kMembers);
   }
   Float_t ratio = (Float_t)((parent->GetLevel()+1.0) / parent->GetLevel());

   // changing the angle between parent object and daughters
   // if center of parent is 0 that is real piramid
   parent->SetWidth(1);
   parent->SetHeight(1);
   parent->SetX(-parent->GetWidth()/2);
   parent->SetY(-parent->GetHeight()/2);

   fMaxRatio = parent->GetVolumeRatio();

   // sorting list of members by size or number of members
   parent->GetMembers()->Sort(kSortDescending);
   Divide(parent->GetMembers(), (parent->GetX()) *ratio, (parent->GetX() + parent->GetWidth())* ratio, (parent->GetY())* ratio, (parent->GetY() + parent->GetHeight())*ratio);

   // sclale all the objects
   Scale(parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if all of nodes can be displayed on scene. Hides redendant nodes.

void TStructViewerGUI::CheckMaxObjects(TStructNode* parent)
{
   UInt_t object = 0;

   TList queue;
   queue.Add(parent);
   TStructNode* node;

   while ((node = (TStructNode*) queue.First() )) {
      object++;

      if (object > fNodePtr->GetMaxObjects() || node->GetLevel() - fNodePtr->GetLevel() >= fNodePtr->GetMaxLevel()) {
         break;
      }

      node->SetVisible(true);

      queue.AddAll(node->GetMembers());
      queue.RemoveFirst();

      fVisibleObjects.Add(node);
   }

   TIter it(&fVisibleObjects);
   TStructNode* member;
   while ((node = (TStructNode*) it() )) {
      if(node->GetLevel() - fNodePtr->GetLevel() == fNodePtr->GetMaxLevel()-1 && node->GetMembersCount() > 0) {
         node->SetCollapsed(true);
         continue;
      }

      TIter memIt(node->GetMembers());
      while ((member = (TStructNode*) memIt() )) {
         if(member->IsVisible() == false) {
            node->SetCollapsed(true);
            break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete window

void TStructViewerGUI::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for default color selsect.
/// Sets default colot to "pixel"

void TStructViewerGUI::ColorSelectedSlot(Pixel_t pixel)
{
   TStructNodeProperty* prop = GetDefaultColor();
   if(prop) {
      prop->SetColor(pixel);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Divides rectangle where the outlining box is placed.

void TStructViewerGUI::Divide(TList* list, Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
   if (list->GetSize() > 1) { // spliting node into two lists
      ULong_t sum1 = 0, sum = 0;

      TStructNode* node;
      TList list1, list2;
      TIter it(list);

      while((node = (TStructNode*) it() )) {
         sum += node->GetVolume();
      }
      it.Reset();
      while((node = (TStructNode*) it() )) {
         if(sum1 >= sum/2.0) {
            list2.Add(node);
         } else {
            sum1 += node->GetVolume();
            list1.Add(node);
         }
      }

      if (!sum) return;
      Float_t ratio = (float)sum1/sum;

      Float_t width = x2 - x1;
      Float_t height = y2 - y1;
      if (width < height) { // vertical split
         Float_t split = y1 + ratio * height;
         Divide(&list1, x1, x2, y1, split);
         Divide(&list2, x1, x2, split, y2);
      } else { // horizontal
         Float_t split = x1 + ratio * width;
         Divide(&list1, x1, split, y1, y2);
         Divide(&list2, split, x2, y1, y2);
      }
   } else if (list->GetSize() == 1) { // divide place to node
      TStructNode* node = (TStructNode*)(list->First());

      node->SetWidth(x2 - x1);
      node->SetHeight(y2 - y1);
      node->SetX(x1);
      node->SetY(y1);

      if (node->GetVolumeRatio() > fMaxRatio) {
         fMaxRatio = node->GetVolumeRatio();
      }

      Float_t ratio = (Float_t)((node->GetLevel()+1.0)/node->GetLevel());
      node->GetMembers()->Sort(kSortDescending);
      Divide(node->GetMembers(), x1*ratio, x2*ratio, y1*ratio, y2*ratio);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Activated when user double click on objects on 3D scene. Sets clicked node to top node
/// and updates scene with camers reset.

void TStructViewerGUI::DoubleClickedSlot()
{
   if (fSelectedObject) {
      if(fSelectedObject == fNodePtr) {
         return;
      }

      fUndoList.Add(fNodePtr);
      fNodePtr = fSelectedObject;
      fUndoButton->SetEnabled(true);

      Update(kTRUE);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Check limits and draws nodes and links

void TStructViewerGUI::Draw(Option_t* /*option*/)
{
   fVolumes.Clear();
   CheckMaxObjects(fNodePtr);

   CalculatePosistion(fNodePtr);
   DrawVolumes(fNodePtr);

   if(fShowLinksCheckButton->GetState() == kButtonDown) {
      DrawLink(fNodePtr);
   }

   UnCheckMaxObjects();
}

////////////////////////////////////////////////////////////////////////////////
/// Recursive method to draw links

void TStructViewerGUI::DrawLink(TStructNode* parent)
{
   if(parent->GetLevel() - fNodePtr->GetLevel() >= fNodePtr->GetMaxLevel()) {
      return;
   }

   if(parent->IsCollapsed()) {
      return;
   }

   TIter it(parent->GetMembers());
   TStructNode* node;
   while((node = (TStructNode*) it())) {
      TPolyLine3D *l = new TPolyLine3D(2);
      l->SetPoint(0 ,node->GetCenter(), node->GetMiddle(), -(node->GetLevel() * fLevelDistanceEntry->GetNumber()));
      l->SetPoint(1 ,parent->GetCenter(), parent->GetMiddle(), -(parent->GetLevel() * fLevelDistanceEntry->GetNumber()));

      l->SetLineColor(GetColor(node));
      l->SetLineWidth(1);
      l->Draw();

      if(!node->IsCollapsed()) {
         DrawLink(node);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates and draws TGeoVolume from given "node"

void TStructViewerGUI::DrawNode(TStructNode* node)
{
   TGeoVolume* vol;

   /*if(node->IsCollapsed())
   {
   //float r = (node->GetWidth() < node->GetHeight() ? 0.5 * node->GetWidth() : 0.5 * node->GetHeight());
   //vol = gGeoManager->MakeTorus(node->GetName(),TStructNode::GetMedium(), 0.75*r, 0, r/4);

   vol = gGeoManager->MakeBox(TString(node->GetName()) + "up",TStructNode::GetMedium(), 0.45*node->GetWidth(), 0.45*node->GetHeight(), (node->GetWidth() < node->GetHeight() ? 0.45 * node->GetWidth() : 0.45 * node->GetHeight()));
   Double_t max = TMath::Max(0.22 * node->GetWidth(), 0.22 * node->GetHeight());
   TGeoVolume* subvol = gGeoManager->MakeTrd2(node->GetName(), TStructNode::GetMedium(), 0, 0.45 * node->GetWidth(), 0, 0.45 * node->GetHeight(), max);
   subvol->SetLineColor(GetColor(node));
   subvol->SetNumber((Int_t)node);
   TGeoTranslation* subtrans = new TGeoTranslation("subtranslation", 0, 0, -max);
   vol->AddNodeOverlap(subvol, 1, subtrans);

   subvol = gGeoManager->MakeTrd2(TString(node->GetName()) + "down", TStructNode::GetMedium(), 0.45 * node->GetWidth(), 0, 0.45 * node->GetHeight(), 0, max);
   subvol->SetLineColor(GetColor(node));
   subvol->SetNumber((Int_t)node);
   subtrans = new TGeoTranslation("subtranslation", 0, 0, max);
   vol->AddNodeOverlap(subvol, 1, subtrans);
   }
   else*/ if(node->GetNodeType() == kCollection) {
      vol = gGeoManager->MakeBox(Form("%s_%d", node->GetName(), fgCounter++), fgMedium, 0.45*node->GetWidth(), 0.45*node->GetHeight(), fBoxHeightEntry->GetNumber());
      // subboxes
      Float_t slices = (Float_t)(node->GetMembersCount());
      if (slices > fMaxSlices) {
         slices = (Float_t)fMaxSlices;
      }

      for (Float_t i = -(slices-1)/2; i < slices/2; i++) {
         TGeoVolume* sub = gGeoManager->MakeBox(Form("%s_%d", node->GetName(), fgCounter++), fgMedium,0.45*node->GetWidth() * 0.7 / slices, 0.45*node->GetHeight(), fBoxHeightEntry->GetNumber());
         sub->SetLineColor(GetColor(node));
         fVolumes.Add((Long_t)sub, (Long_t)node);
         TGeoTranslation* subtrans = new TGeoTranslation("subtranslation", i * node->GetWidth() / slices, 0, 0);
         vol->AddNodeOverlap(sub, 1, subtrans);
      }
   } else {
      vol = gGeoManager->MakeBox(Form("%s_%d", node->GetName(), fgCounter++), fgMedium, 0.45*node->GetWidth(), 0.45*node->GetHeight(), fBoxHeightEntry->GetNumber());
   }

   vol->SetLineColor(GetColor(node));
   vol->SetLineWidth(1);

   TGeoTranslation* trans = new TGeoTranslation("translation", node->GetCenter(), node->GetMiddle(), -(node->GetLevel() * fLevelDistanceEntry->GetNumber()));
   fVolumes.Add((Long_t)vol, (Long_t)node);

   fTopVolume->AddNode(vol,1, trans);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursive method to draw GeoVolumes

void TStructViewerGUI::DrawVolumes(TStructNode* parent)
{
   if(parent->GetLevel() - fNodePtr->GetLevel() >= fNodePtr->GetMaxLevel()) {
      return;
   }

   DrawNode(parent);

   if(parent->IsCollapsed()) {
      return;
   }

   TIter nextVis(parent->GetMembers());
   TStructNode* node;
   while((node = (TStructNode*)nextVis())) {
      DrawVolumes(node);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to property associated with node "node". If property is not found
/// then it returns default property

TStructNodeProperty* TStructViewerGUI::FindNodeProperty(TStructNode* node)
{
   TIter it(fColors);
   TStructNodeProperty* prop;
   while ((prop = (TStructNodeProperty*) it() )) {
      TString propName(prop->GetName());
      if (propName.EndsWith("+")) {

         if (TClass* cl = TClass::GetClass(node->GetTypeName())) {
            propName.Remove(propName.Length()-1, 1);
            if (cl->InheritsFrom(propName.Data())) {
               return prop;
            }
         }
      } else {
         if (propName == TString(node->GetTypeName())) {
            return prop;
         }
      }
   }

   return (TStructNodeProperty*)fColors->Last();
}

//________________________________________________________________________`
TCanvas* TStructViewerGUI::GetCanvas()
{
   // Returns canvas used to keep TGeoVolumes

   return fCanvas;
}
////////////////////////////////////////////////////////////////////////////////
/// Returns color form fColors for given "node"

Int_t TStructViewerGUI::GetColor(TStructNode* node)
{
   TStructNodeProperty* prop = FindNodeProperty(node);
   if (prop) {
      return prop->GetColor().GetNumber();
   }

   return 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default color for nodes

TStructNodeProperty* TStructViewerGUI::GetDefaultColor()
{
   return ((TStructNodeProperty*)(fColors->Last()));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if links are visible, otherwise return false.

Bool_t TStructViewerGUI::GetLinksVisibility() const
{
   if (fShowLinksCheckButton->GetState() == kButtonDown) {
      return true;
   } else {
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns top node pointer

TStructNode* TStructViewerGUI::GetNodePtr() const
{
   return fNodePtr;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle events. Sets fMouseX and fMouseY when user move a mouse over viewer and hides ToolTip

void TStructViewerGUI::GLWidgetProcessedEventSlot(Event_t* event)
{
   switch (event->fType) {
      case kMotionNotify:
         fMouseX = event->fXRoot + 15;
         fMouseY = event->fYRoot + 15;
         break;

      case kButtonPress:
         fToolTip->Hide();
         if (fSelectedObject) {
            UpdateLabels(fSelectedObject);
            fEditor->SetModel(fSelectedObject);
         }
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emmited when user changes distance between levels

void TStructViewerGUI::LevelDistValueSetSlot(Long_t /* dist */)
{
   if(fAutoRefesh->IsOn()) {
      Update(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// MouseOver slot. Activated when user out mouse over object on scene.
/// Sets ToolTip and updates labels

void TStructViewerGUI::MouseOverSlot(TGLPhysicalShape* shape)
{
   fToolTip->Hide();
   fSelectedObject = NULL;
   if (shape && shape->GetLogical()) {
      fSelectedObject =  (TStructNode*)(shape->GetLogical()->ID());
      if (fSelectedObject) {
         if (fSelectedObject->IsA()->InheritsFrom(TPolyLine3D::Class())) {
            fSelectedObject = NULL;
            return;
         }
         Long_t shapeID  = (Long_t)(shape->GetLogical()->ID());
         Long_t volValue = (Long_t)fVolumes.GetValue(shapeID);
         fSelectedObject = (TStructNode*)volValue;

         fToolTip->SetText(TString(fSelectedObject->GetName()) + "\n" + fSelectedObject->GetTypeName());
         fToolTip->SetPosition(fMouseX, fMouseY);
         fToolTip->Reset();
         UpdateLabels(fSelectedObject);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Activated when user click Redo button. Repeat last Undo action.

void TStructViewerGUI::RedoButtonSlot()
{
   fUndoList.Add(fNodePtr);
   fUndoButton->SetEnabled(true);
   fNodePtr = (TStructNode*) fRedoList.Last();
   fRedoList.RemoveLast();
   if (!fRedoList.First()) {
      fRedoButton->SetEnabled(false);
   }
   Update(kTRUE);
   UpdateLabels(fNodePtr);
}

////////////////////////////////////////////////////////////////////////////////
/// Resets camera

void TStructViewerGUI::ResetButtonSlot()
{
   fGLViewer->UpdateScene();
   fGLViewer->ResetCurrentCamera();
}

////////////////////////////////////////////////////////////////////////////////
/// Recursive method to scaling all modes on scene. We have to scale nodes to get real ratio between nodes.
/// Uses fMaxRatio.

void TStructViewerGUI::Scale(TStructNode* parent)
{
   // newRatio = sqrt(ratio/maxratio)
   Float_t newRatio = (Float_t)(TMath::Sqrt(parent->GetRelativeVolumeRatio()/fMaxRatio));
   // set left top conner in the center
   parent->SetX(parent->GetX() + parent->GetWidth()/2);
   parent->SetY(parent->GetY() + parent->GetHeight()/2);
   // set new size
   Float_t min = (Float_t)TMath::Min(parent->GetWidth(), parent->GetHeight());
   parent->SetWidth(parent->GetWidth() * newRatio);
   parent->SetHeight(parent->GetHeight() * newRatio);
   // fit the ratio -> height to width
   Float_t sqrt = (Float_t)(TMath::Sqrt(parent->GetWidth() * parent->GetHeight()));
   // it's a square
   if (min > sqrt) {
      parent->SetWidth(sqrt);
      parent->SetHeight(sqrt);
   } else { // it's rectangle
      if (parent->GetHeight() > parent->GetWidth()) {
         parent->SetWidth(min);
         parent->SetHeight(sqrt * sqrt / min);
      } else {
         parent->SetWidth(sqrt * sqrt / min);
         parent->SetHeight(min);
      }
   }
   // move left top corner
   parent->SetX(parent->GetX() - parent->GetWidth()/2);
   parent->SetY(parent->GetY() - parent->GetHeight()/2);

   // scale others nodes
   TStructNode* node;
   TIter it(parent->GetMembers());
   while ((node = (TStructNode*) it() )) {
      Scale(node);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets top node pointer and updates view

void TStructViewerGUI::SetNodePtr(TStructNode* val)
{
   fNodePtr = val;
   Update(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets links visibility to "visible"

void TStructViewerGUI::SetLinksVisibility(Bool_t visible)
{
   if (visible) {
      fShowLinksCheckButton->SetState(kButtonDown);
   } else {
      fShowLinksCheckButton->SetState(kButtonUp);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets pointer given in fPointerTestEntry to the main pointer

void TStructViewerGUI::SetPointerButtonSlot()
{
   void* obj = (void*)gROOT->ProcessLine(fPointerTextEntry->GetText());
   fParent->SetPointer(obj, fPointerTypeTextEntry->GetText());
}

////////////////////////////////////////////////////////////////////////////////
/// Changes links visibility and refresh view.

void TStructViewerGUI::ShowLinksToggled(Bool_t /*on*/)
{
   if (fAutoRefesh->IsOn()) {
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Shows hidden nodes

void TStructViewerGUI::UnCheckMaxObjects()
{
   TStructNode* node;
   TIter it(&fVisibleObjects);

   while ((node = (TStructNode*) it() )) {
      node->SetCollapsed(false);
      node->SetVisible(false);
   }

   fVisibleObjects.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Updates view. Clear all the nodes, call draw function and update scene. Doesn't reset camera.

void TStructViewerGUI::Update(Bool_t resetCamera)
{
   if (!fNodePtr) {
      return;
   }

   fCanvas->GetListOfPrimitives()->Clear();
   fTopVolume->ClearNodes();
   Draw();
   fCanvas->GetListOfPrimitives()->Add(fTopVolume);
   fGLViewer->UpdateScene();

   if(resetCamera) {
      fGLViewer->ResetCurrentCamera();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update button slot. Updates scene

void TStructViewerGUI::UpdateButtonSlot()
{
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh information in labels when user put mouse over object

void TStructViewerGUI::UpdateLabels(TStructNode* node)
{
   fNodeNameLabel->SetText(node->GetName());
   fNodeTypelabel->SetText(node->GetTypeName());

   TString name = "Members: ";
   name += node->GetMembersCount();
   fMembersCountLabel->SetText(name);
   name = "All members: ";
   name += node->GetAllMembersCount();
   fAllMembersCountLabel->SetText(name);
   name = "Level: ";
   name += node->GetLevel();
   fLevelLabel->SetText(name);
   name = "Size: ";
   name += node->GetSize();
   fSizeLabel->SetText(name);
   name = "Total size: ";
   name += node->GetTotalSize();
   fTotalSizeLabel->SetText(name);
}

////////////////////////////////////////////////////////////////////////////////
/// UndoButton Slot. Activated when user press Undo button. Restore last top node pointer.

void TStructViewerGUI::UndoButtonSlot()
{
   fRedoList.Add(fNodePtr);
   fRedoButton->SetEnabled(true);
   fNodePtr = (TStructNode*) fUndoList.Last();
   fUndoList.RemoveLast();
   if (!fUndoList.First()) {
      fUndoButton->SetEnabled(false);
   }
   Update(kTRUE);
   UpdateLabels(fNodePtr);
}

////////////////////////////////////////////////////////////////////////////////
/// Activated when user press radio button

void TStructViewerGUI::ScaleByChangedSlot()
{
    if (fAutoRefesh->IsOn()) {
       Update();
    }
}
