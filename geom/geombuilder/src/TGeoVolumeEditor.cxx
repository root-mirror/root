// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoVolumeEditor
\ingroup Geometry_builder

Editor for geometry volumes and assemblies of volumes. Besides the volume
name and line attributes, a TGeoVolume has the following editable categories
split vertically by a shutter:

  - Properties: one can edit the shape and medium components from here. It is
    also possible to change the existing ones.
  - Daughters: the main category allowing defining, editing, removing or
    positioning daughter volumes inside the current edited volume. To add a
    daughter, one needs to select first a volume and a matrix. Currently no check
    is performed if the daughter volume creates an extrusion (illegal for tracking).
    To remove or change the position of an existing daughter, one should simply
    select the desired daughter from the combo box with the existing ones, then
    simply click the appropriate button.
  - Visualization: One can set the visibility of the volume and of its daughters,
    set the visibility depth and the view type. Selecting "All" will draw the
    volume and all visible daughters down to the selected level starting from the
    edited volume. Selecting "Leaves" will draw just the deepest daughters within
    the selected visibility level, without displaying the containers, while "Only"
    will just draw the edited volume.
  - Division: The category becomes active only if there are no daughters of the
    edited volume added by normal positioning (e.g. from <Daughters> category). The
    minimum allowed starting value for the selected division axis is automatically
    selected, while the slicing step is set to 0 - meaning that only the number
    of slices matter.
*/

#include "TGeoVolumeEditor.h"
#include "TGeoVolume.h"
#include "TGeoPatternFinder.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGTab.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGShutter.h"
#include "TG3DLine.h"
#include "TGeoTabManager.h"
#include "TGedEditor.h"

ClassImp(TGeoVolumeEditor);

enum ETGeoVolumeWid {
   kVOL_NAME, kVOL_TITLE, kVOL_SHAPE_SELECT, kVOL_MEDIA_SELECT, kVOL_NODE_SELECT,
   kVOL_VOL_SELECT, kVOL_MATRIX_SELECT, kVOL_EDIT_SHAPE, kVOL_EDIT_MEDIUM, kVOL_NODEID,
   kVOL_APPLY, kVOL_CANCEL, kVOL_UNDO, kVOL_VISLEVEL, kVOL_DIVSTART, kVOL_DIVEND,
   kVOL_DIVSTEP, kVOL_DIVN, kCAT_GENERAL, kCAT_DAUGHTERS, kCAT_DIVISION, kCAT_VIS,
   kDIV_NAME
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for volume editor.

TGeoVolumeEditor::TGeoVolumeEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fGeometry = 0;
   fVolume   = 0;

   fIsModified = kFALSE;
   fIsAssembly = kFALSE;
   fIsDivided = kFALSE;

   // TGShutter for categories
   fCategories = new TGShutter(this, kSunkenFrame);
   TGCompositeFrame *container, *f1;
   Pixel_t color;
   TGLabel *label;

   // General settings
   TGShutterItem *si = new TGShutterItem(fCategories, new TGHotString("Properties"),kCAT_GENERAL);
   container = (TGCompositeFrame*)si->GetContainer();
   container->SetBackgroundColor(GetDefaultFrameBackground());
   fCategories->AddItem(si);

   // TextEntry for volume name
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Volume name"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   fVolumeName = new TGTextEntry(container, "", kVOL_NAME);
   fVolumeName->SetDefaultSize(135, fVolumeName->GetDefaultHeight());
   fVolumeName->SetToolTipText("Enter the volume name");
   container->AddFrame(fVolumeName, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 3, 1, 2, 5));

   // Current shape
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Shape and medium"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 10, 0));
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame);
   fSelectedShape = 0;
   fLSelShape = new TGLabel(f1, "Select shape");
   gClient->GetColorByName("#0000ff", color);
   fLSelShape->SetTextColor(color);
   fLSelShape->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelShape, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelShape = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kVOL_SHAPE_SELECT);
   fBSelShape->SetToolTipText("Replace with one of the existing shapes");
   fBSelShape->Associate(this);
   f1->AddFrame(fBSelShape, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditShape = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditShape, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditShape->SetToolTipText("Edit selected shape");
   fEditShape->Associate(this);
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 2, 0, 0));

   // Current medium
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame);
   fSelectedMedium = 0;
   fLSelMedium = new TGLabel(f1, "Select medium");
   gClient->GetColorByName("#0000ff", color);
   fLSelMedium->SetTextColor(color);
   fLSelMedium->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelMedium, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelMedium = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kVOL_MEDIA_SELECT);
   fBSelMedium->SetToolTipText("Replace with one of the existing media");
   fBSelMedium->Associate(this);
   f1->AddFrame(fBSelMedium, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMedium = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditMedium, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMedium->SetToolTipText("Edit selected medium");
   fEditMedium->Associate(this);
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 2, 0, 0));

   // List of daughters
   si = new TGShutterItem(fCategories, new TGHotString("Daughters"),kCAT_DAUGHTERS);
   container = (TGCompositeFrame*)si->GetContainer();
   container->SetBackgroundColor(GetDefaultFrameBackground());
   fCategories->AddItem(si);

   // Existing daughters
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Existing daughters"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame | kRaisedFrame);
   fNodeList = new TGComboBox(f1, kVOL_NODE_SELECT);
   fNodeList->Resize(100, fVolumeName->GetDefaultHeight());
   fNodeList->Associate(this);
   f1->AddFrame(fNodeList, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));
   // Buttons for editing matrix and removing node
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame | kSunkenFrame | kFixedWidth);
   fEditMatrix = new TGTextButton(f1, "Position");
   f1->AddFrame(fEditMatrix, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   fEditMatrix->SetToolTipText("Edit the position of selected node");
   fEditMatrix->Associate(this);
   fRemoveNode = new TGTextButton(f1, "Remove");
   f1->AddFrame(fRemoveNode, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   fRemoveNode->SetToolTipText("Remove the selected node. Cannot undo !)");
   fRemoveNode->Associate(this);
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

   // Adding daughters
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Add daughter"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 10, 0));

   // Select from existing volumes
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame | kFixedWidth);
   fSelectedVolume = 0;
   fLSelVolume = new TGLabel(f1, "Select volume");
   gClient->GetColorByName("#0000ff", color);
   fLSelVolume->SetTextColor(color);
   fLSelVolume->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelVolume, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelVolume = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kVOL_VOL_SELECT);
   fBSelVolume->SetToolTipText("Select one of the existing volumes");
   fBSelVolume->Associate(this);
   f1->AddFrame(fBSelVolume, new TGLayoutHints(kLHintsRight, 1, 1, 2, 2));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

   // Matrix selection for nodes
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame | kFixedWidth);
   fSelectedMatrix = 0;
   fLSelMatrix = new TGLabel(f1, "Select matrix");
   gClient->GetColorByName("#0000ff", color);
   fLSelMatrix->SetTextColor(color);
   fLSelMatrix->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelMatrix, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelMatrix = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kVOL_MATRIX_SELECT);
   fBSelMatrix->SetToolTipText("Select one of the existing matrices");
   fBSelMatrix->Associate(this);
   f1->AddFrame(fBSelMatrix, new TGLayoutHints(kLHintsRight, 1, 1, 2, 2));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

   // Copy number
   f1 = new TGCompositeFrame(container, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Node id"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fCopyNumber = new TGNumberEntry(f1, 0., 5, kVOL_NODEID);
   fCopyNumber->SetNumStyle(TGNumberFormat::kNESInteger);
   fCopyNumber->SetNumAttr(TGNumberFormat::kNEANonNegative);
   fCopyNumber->Resize(20,fCopyNumber->GetDefaultHeight());
   TGTextEntry *nef = (TGTextEntry*)fCopyNumber->GetNumberEntry();
   nef->SetToolTipText("Enter node copy number");
   fCopyNumber->Associate(this);
   f1->AddFrame(fCopyNumber, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   fAddNode = new TGTextButton(f1, "Add");
   f1->AddFrame(fAddNode, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   fAddNode->Associate(this);
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

   // Visualization
   si = new TGShutterItem(fCategories, new TGHotString("Visualization"),kCAT_VIS);
   container = (TGCompositeFrame*)si->GetContainer();
   container->SetBackgroundColor(GetDefaultFrameBackground());
   fCategories->AddItem(si);

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Visibility"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
//   gClient->GetColorByName("#ff0000", color);
//   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame | kDoubleBorder);
   fBVis[0] = new TGCheckButton(f1, "Volume");
   fBVis[1] = new TGCheckButton(f1, "Nodes");
   f1->AddFrame(fBVis[0], new TGLayoutHints(kLHintsLeft, 2, 2, 0 ,0));
   f1->AddFrame(fBVis[1], new TGLayoutHints(kLHintsRight, 2, 2, 0 ,0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Depth"), new TGLayoutHints(kLHintsLeft, 2, 2, 4, 0));
//   gClient->GetColorByName("#0000ff", color);
//   label->SetTextColor(color);
   fEVisLevel = new TGNumberEntry(f1, 0, 5, kVOL_VISLEVEL);
   fEVisLevel->SetNumStyle(TGNumberFormat::kNESInteger);
   fEVisLevel->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEVisLevel->Resize(40,fEVisLevel->GetDefaultHeight());
   nef = (TGTextEntry*)fEVisLevel->GetNumberEntry();
   nef->SetToolTipText("Set visibility level here");
   fEVisLevel->SetNumber(3);
   fEVisLevel->Associate(this);
   f1->AddFrame(fEVisLevel, new TGLayoutHints(kLHintsLeft, 2, 2, 0 ,0));
   fBAuto = new TGCheckButton(f1,"Auto");
   f1->AddFrame(fBAuto, new TGLayoutHints(kLHintsRight, 0, 0, 2, 0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   TString stitle = "View";
   TGButtonGroup *bg = new TGVButtonGroup(container, stitle);
   fBView[0] = new TGRadioButton(bg, "All");
   fBView[1] = new TGRadioButton(bg, "Leaves");
   fBView[2] = new TGRadioButton(bg, "Only");
   bg->SetRadioButtonExclusive();
   bg->Show();
   container->AddFrame(bg, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame | kDoubleBorder);
   fBRaytrace = new TGCheckButton(f1,"Raytrace");
   f1->AddFrame(fBRaytrace, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Division
   si = new TGShutterItem(fCategories, new TGHotString("Division"),kCAT_DIVISION);
   container = (TGCompositeFrame*)si->GetContainer();
   container->SetBackgroundColor(GetDefaultFrameBackground());
   fCategories->AddItem(si);
   // TextEntry for division name
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Division name"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   fDivName = new TGTextEntry(container, new TGTextBuffer(50), kDIV_NAME);
   fDivName->Resize(135, fVolumeName->GetDefaultHeight());
   fDivName->SetToolTipText("Enter the volume name");
   container->AddFrame(fDivName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
   // Axis selection
   stitle = "Axis";
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   bg = new TGVButtonGroup(f1, stitle);
   fBDiv[0] = new TGRadioButton(bg, "Axis 1");
   fBDiv[1] = new TGRadioButton(bg, "Axis 2");
   fBDiv[2] = new TGRadioButton(bg, "Axis 3");
   bg->Insert(fBDiv[0]);
   bg->Insert(fBDiv[1]);
   bg->Insert(fBDiv[2]);
   bg->SetRadioButtonExclusive();
   f1->AddFrame(bg, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
   fApplyDiv = new TGTextButton(f1, "Apply");
   fApplyDiv->SetToolTipText("Apply new division settings");
   f1->AddFrame(fApplyDiv, new TGLayoutHints(kLHintsRight, 0, 2, 30, 0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
   // Division range
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Division parameters"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
//   gClient->GetColorByName("#ff0000", color);
//   label->SetTextColor(color);
   container->AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "From"), new TGLayoutHints(kLHintsLeft, 2, 2, 4, 0));
//   gClient->GetColorByName("#0000ff", color);
//   label->SetTextColor(color);
   fEDivFrom = new TGNumberEntry(f1, 0, 5, kVOL_DIVSTART);
//   fEDivFrom->SetNumStyle(TGNumberFormat::kNESInteger);
//   fEDivFrom->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDivFrom->Resize(100,fEDivFrom->GetDefaultHeight());
   nef = (TGTextEntry*)fEDivFrom->GetNumberEntry();
   nef->SetToolTipText("Set start value");
   fEDivFrom->Associate(this);
   f1->AddFrame(fEDivFrom, new TGLayoutHints(kLHintsRight, 2, 2, 0 ,0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Step"), new TGLayoutHints(kLHintsLeft, 2, 2, 4, 0));
//   gClient->GetColorByName("#0000ff", color);
//   label->SetTextColor(color);
   fEDivStep = new TGNumberEntry(f1, 0, 5, kVOL_DIVSTEP);
//   fEDivFrom->SetNumStyle(TGNumberFormat::kNESInteger);
   fEDivStep->SetNumAttr(TGNumberFormat::kNEANonNegative);
   fEDivStep->Resize(100,fEDivStep->GetDefaultHeight());
   nef = (TGTextEntry*)fEDivStep->GetNumberEntry();
   nef->SetToolTipText("Set division step");
   fEDivStep->Associate(this);
   f1->AddFrame(fEDivStep, new TGLayoutHints(kLHintsRight, 2, 2, 0 ,0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   f1 = new TGCompositeFrame(container, 155, 10, kHorizontalFrame |kFixedWidth);
   f1->AddFrame(label = new TGLabel(f1, "Nslices"), new TGLayoutHints(kLHintsLeft, 2, 2, 4, 0));
//   gClient->GetColorByName("#0000ff", color);
//   label->SetTextColor(color);
   fEDivN = new TGNumberEntry(f1, 0, 5, kVOL_DIVN);
   fEDivN->SetNumStyle(TGNumberFormat::kNESInteger);
   fEDivN->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDivN->Resize(100,fEDivN->GetDefaultHeight());
   nef = (TGTextEntry*)fEDivN->GetNumberEntry();
   nef->SetToolTipText("Set number of slices");
   fEDivN->Associate(this);
   f1->AddFrame(fEDivN, new TGLayoutHints(kLHintsRight, 2, 2, 0 ,0));
   container->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));


   fCategories->Resize(163,340);
   AddFrame(fCategories, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 0, 0, 4, 4));

   fCategories->Layout();
   fCategories->SetDefaultSize(GetDefaultWidth(), GetDefaultHeight());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVolumeEditor::~TGeoVolumeEditor()
{
   TGCompositeFrame *cont;
   cont = (TGCompositeFrame*)fCategories->GetItem("Properties")->GetContainer();
   TGeoTabManager::Cleanup(cont);
   fCategories->GetItem("Properties")->SetCleanup(0);
   cont = (TGCompositeFrame*)fCategories->GetItem("Daughters")->GetContainer();
   TGeoTabManager::Cleanup(cont);
   fCategories->GetItem("Daughters")->SetCleanup(0);
   cont = (TGCompositeFrame*)fCategories->GetItem("Visualization")->GetContainer();
   TGeoTabManager::Cleanup(cont);
   fCategories->GetItem("Visualization")->SetCleanup(0);
   cont = (TGCompositeFrame*)fCategories->GetItem("Division")->GetContainer();
   TGeoTabManager::Cleanup(cont);
   fCategories->GetItem("Division")->SetCleanup(0);

   delete fBView[0]; delete fBView[1]; delete fBView[2];
   delete fBDiv [0]; delete fBDiv [1]; delete fBDiv [2];
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TGeoVolumeEditor::ConnectSignals2Slots()
{
   fVolumeName->Connect("TextChanged(const char *)", "TGeoVolumeEditor", this, "DoVolumeName()");
   fDivName->Connect("TextChanged(const char *)", "TGeoVolumeEditor", this, "DoDivName()");
   fEditMedium->Connect("Clicked()", "TGeoVolumeEditor", this, "DoEditMedium()");
   fEditShape->Connect("Clicked()", "TGeoVolumeEditor", this, "DoEditShape()");
   fEditMatrix->Connect("Clicked()", "TGeoVolumeEditor", this, "DoEditMatrix()");
   fAddNode->Connect("Clicked()", "TGeoVolumeEditor", this, "DoAddNode()");
   fRemoveNode->Connect("Clicked()", "TGeoVolumeEditor", this, "DoRemoveNode()");
   fBSelShape->Connect("Clicked()", "TGeoVolumeEditor", this, "DoSelectShape()");
   fBSelMedium->Connect("Clicked()", "TGeoVolumeEditor", this, "DoSelectMedium()");
   fBSelVolume->Connect("Clicked()", "TGeoVolumeEditor", this, "DoSelectVolume()");
   fBSelMatrix->Connect("Clicked()", "TGeoVolumeEditor", this, "DoSelectMatrix()");
   fBVis[0]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoVisVolume()");
   fBVis[1]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoVisDaughters()");
   fBAuto->Connect("Clicked()", "TGeoVolumeEditor", this, "DoVisAuto()");
   fEVisLevel->Connect("ValueSet(Long_t)", "TGeoVolumeEditor", this, "DoVisLevel()");
   fBView[0]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoViewAll()");
   fBView[1]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoViewLeaves()");
   fBView[2]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoViewOnly()");
   fBDiv[0]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoDivSelAxis()");
   fBDiv[1]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoDivSelAxis()");
   fBDiv[2]->Connect("Clicked()", "TGeoVolumeEditor", this, "DoDivSelAxis()");
   fEDivFrom->Connect("ValueSet(Long_t)", "TGeoVolumeEditor", this, "DoDivFromTo()");
   fEDivStep->Connect("ValueSet(Long_t)", "TGeoVolumeEditor", this, "DoDivStep()");
   fEDivN->Connect("ValueSet(Long_t)", "TGeoVolumeEditor", this, "DoDivN()");
   fBRaytrace->Connect("Clicked()", "TGeoVolumeEditor", this, "DoRaytrace()");
   fApplyDiv->Connect("Clicked()", "TGeoVolumeEditor", this, "DoApplyDiv()");
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to the picked volume.

void TGeoVolumeEditor::SetModel(TObject* obj)
{
   if (obj == 0 || !obj->InheritsFrom(TGeoVolume::Class())) {
      SetActive(kFALSE);
      return;
   }
   fVolume = (TGeoVolume*)obj;
   fGeometry = fVolume->GetGeoManager();
   const char *vname = fVolume->GetName();
   fVolumeName->SetText(vname);
   fSelectedShape = fVolume->GetShape();
   if (fSelectedShape) fLSelShape->SetText(fSelectedShape->GetName());
   fSelectedMedium = fVolume->GetMedium();
   if (fSelectedMedium) fLSelMedium->SetText(fSelectedMedium->GetName());

   fNodeList->RemoveEntries(0, fNodeList->GetNumberOfEntries()+1);
   TIter next2(fVolume->GetNodes());
   TGeoNode *node;
   Int_t icrt = 0;
   while ((node=(TGeoNode*)next2()))
      fNodeList->AddEntry(node->GetName(), icrt++);
   fNodeList->Select(0);
   fCopyNumber->SetNumber(fVolume->GetNdaughters()+1);
   if (!fVolume->GetNdaughters() || fVolume->GetFinder()) {
      fEditMatrix->SetEnabled(kFALSE);
      fRemoveNode->SetEnabled(kFALSE);
   } else {
      fEditMatrix->SetEnabled(kTRUE);
      fRemoveNode->SetEnabled(kTRUE);
   }
   if (!fSelectedVolume) fAddNode->SetEnabled(kFALSE);
   if (fVolume->IsAssembly()) {
      fBSelShape->SetEnabled(kFALSE);
      fBSelMedium->SetEnabled(kFALSE);
   }
   fBVis[0]->SetState((fVolume->IsVisible())?kButtonDown:kButtonUp);
   fBVis[1]->SetState((fVolume->IsVisibleDaughters())?kButtonDown:kButtonUp);
   fBView[0]->SetState((fVolume->IsVisContainers())?kButtonDown:kButtonUp, kTRUE);
   fBView[1]->SetState((fVolume->IsVisLeaves())?kButtonDown:kButtonUp, kTRUE);
   fBView[2]->SetState((fVolume->IsVisOnly())?kButtonDown:kButtonUp, kTRUE);
   fBRaytrace->SetState((fVolume->IsRaytracing())?kButtonDown:kButtonUp);
   fBAuto->SetState((fGeometry->GetVisLevel())?kButtonUp:kButtonDown);
   fEVisLevel->SetNumber(fGeometry->GetVisLevel());
   fApplyDiv->SetEnabled(kFALSE);
   if ((!fVolume->GetFinder() && fVolume->GetNdaughters()) || fVolume->IsAssembly()) {
      fCategories->GetItem("Division")->GetButton()->SetEnabled(kFALSE);
   } else {
      fCategories->GetItem("Division")->GetButton()->SetEnabled(kTRUE);
      Double_t start=0., step=0., end = 0.;
      Int_t ndiv = 2, iaxis = 1;
      TString axis_name;
      for (Int_t i=0; i<3; i++) {
         axis_name = fVolume->GetShape()->GetAxisName(i+1);
         fBDiv[i]->SetText(axis_name);
      }

      if (fVolume->GetFinder()) {
         fDivName->SetText(fVolume->GetNode(0)->GetVolume()->GetName());
         iaxis = fVolume->GetFinder()->GetDivAxis();
         start = fVolume->GetFinder()->GetStart();
         step = fVolume->GetFinder()->GetStep();
         ndiv = fVolume->GetFinder()->GetNdiv();
      } else {
         fDivName->SetText("Enter name");
         fSelectedShape->GetAxisRange(iaxis,start,end);
         step = 0;
      }
      fBDiv[iaxis-1]->SetState(kButtonDown, kTRUE);
      fEDivFrom->SetNumber(start);
      fEDivStep->SetNumber(step);
      fEDivN->SetNumber(ndiv);
   }

   if (fInit) ConnectSignals2Slots();
   SetActive();
   if (GetParent()==fTabMgr->GetVolumeTab()) fTab->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Add editors to fGedFrame and exclude TLineEditor.

void TGeoVolumeEditor::ActivateBaseClassEditors(TClass* cl)
{
   fGedEditor->ExcludeClassEditor(TAttFill::Class());
   TGedFrame::ActivateBaseClassEditors(cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Modify volume name.

void TGeoVolumeEditor::DoVolumeName()
{
   fVolume->SetName(fVolumeName->GetText());
}

////////////////////////////////////////////////////////////////////////////////
/// Select a new shape.

void TGeoVolumeEditor::DoSelectShape()
{
   TGeoShape *shape = fSelectedShape;
   new TGeoShapeDialog(fBSelShape, gClient->GetRoot(), 200,300);
   fSelectedShape = (TGeoShape*)TGeoShapeDialog::GetSelected();
   if (fSelectedShape) fLSelShape->SetText(fSelectedShape->GetName());
   else fSelectedShape = shape;
}

////////////////////////////////////////////////////////////////////////////////
/// Select a new medium.

void TGeoVolumeEditor::DoSelectMedium()
{
   TGeoMedium *med = fSelectedMedium;
   new TGeoMediumDialog(fBSelMedium, gClient->GetRoot(), 200,300);
   fSelectedMedium = (TGeoMedium*)TGeoMediumDialog::GetSelected();
   if (fSelectedMedium) fLSelMedium->SetText(fSelectedMedium->GetName());
   else fSelectedMedium = med;
}

////////////////////////////////////////////////////////////////////////////////
/// Select a matrix for positioning.

void TGeoVolumeEditor::DoSelectMatrix()
{
   TGeoMatrix *matrix = fSelectedMatrix;
   new TGeoMatrixDialog(fBSelMatrix, gClient->GetRoot(), 200,300);
   fSelectedMatrix = (TGeoMatrix*)TGeoMatrixDialog::GetSelected();
   if (fSelectedMatrix) fLSelMatrix->SetText(fSelectedMatrix->GetName());
   else fSelectedMatrix = matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Select a daughter volume.

void TGeoVolumeEditor::DoSelectVolume()
{
   TGeoVolume *vol = fSelectedVolume;
   new TGeoVolumeDialog(fBSelVolume, gClient->GetRoot(), 200,300);
   fSelectedVolume = (TGeoVolume*)TGeoVolumeDialog::GetSelected();
   if (fSelectedVolume) fLSelVolume->SetText(fSelectedVolume->GetName());
   else fSelectedVolume = vol;
   if (fSelectedVolume)
      fAddNode->SetEnabled(kTRUE);
}


////////////////////////////////////////////////////////////////////////////////
/// Edit the shape of the volume.

void TGeoVolumeEditor::DoEditShape()
{
   fTabMgr->GetShapeEditor(fVolume->GetShape());
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the medium of the volume.

void TGeoVolumeEditor::DoEditMedium()
{
   fTabMgr->GetMediumEditor(fVolume->GetMedium());
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the position of the selected node.

void TGeoVolumeEditor::DoEditMatrix()
{
   if (!fVolume->GetNdaughters()) return;
   Int_t i = fNodeList->GetSelected();
   if (i<0) return;
   fTabMgr->GetMatrixEditor(fVolume->GetNode(i)->GetMatrix());
}

////////////////////////////////////////////////////////////////////////////////
/// Add a daughter.

void TGeoVolumeEditor::DoAddNode()
{
   if (!fSelectedVolume || fVolume->GetFinder()) return;
   Int_t icopy = fCopyNumber->GetIntNumber();
   fVolume->AddNode(fSelectedVolume, icopy, fSelectedMatrix);
   Int_t nd = fVolume->GetNdaughters();
   fNodeList->AddEntry(fVolume->GetNode(nd-1)->GetName(), nd-1);
   fNodeList->Select(nd-1);
   fCopyNumber->SetNumber(nd+1);
   if (fSelectedMatrix) fEditMatrix->SetEnabled(kTRUE);
   fRemoveNode->SetEnabled(kTRUE);
   fGeometry->SetTopVisible();
   fEditMatrix->SetEnabled(kTRUE);
   fRemoveNode->SetEnabled(kTRUE);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a daughter.

void TGeoVolumeEditor::DoRemoveNode()
{
   if (!fVolume->GetNdaughters() || fVolume->GetFinder()) {
      fRemoveNode->SetEnabled(kFALSE);
      fEditMatrix->SetEnabled(kFALSE);
      return;
   }
   Int_t i = fNodeList->GetSelected();
   if (i<0) return;
   fVolume->RemoveNode(fVolume->GetNode(i));
   fNodeList->RemoveEntries(0, fNodeList->GetNumberOfEntries()+1);
   TIter next(fVolume->GetNodes());
   TGeoNode *node;
   i = 0;
   while ((node=(TGeoNode*)next()))
      fNodeList->AddEntry(node->GetName(), i++);
   fNodeList->Select(0);
   fCopyNumber->SetNumber(fVolume->GetNdaughters()+1);
   if (!fVolume->GetNdaughters()) {
      fRemoveNode->SetEnabled(kFALSE);
      fEditMatrix->SetEnabled(kFALSE);
      fCategories->GetItem("Division")->GetButton()->SetEnabled(kTRUE);
      Double_t start=0., step=0., end=0.;
      Int_t ndiv = 2, iaxis = 1;
      fSelectedShape->GetAxisRange(iaxis,start,end);
      step = end-start;
      fBDiv[iaxis-1]->SetState(kButtonDown, kTRUE);
      fEDivFrom->SetNumber(start);
      fEDivStep->SetNumber(step);
      fEDivN->SetNumber(ndiv);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting volume visible/invisible.

void TGeoVolumeEditor::DoVisVolume()
{
   Bool_t on = (fBVis[0]->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (fVolume->IsVisible() == on) return;
   fVolume->SetVisibility(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting daughters visible/invisible.

void TGeoVolumeEditor::DoVisDaughters()
{
   Bool_t on = (fBVis[1]->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (fVolume->IsVisibleDaughters() == on) return;
   fVolume->VisibleDaughters(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting visibility depth auto.

void TGeoVolumeEditor::DoVisAuto()
{
   Bool_t on = (fBAuto->GetState()==kButtonDown)?kTRUE:kFALSE;
   if ((fGeometry->GetVisLevel()==0) == on) return;
   if (on) fGeometry->SetVisLevel(0);
   else    fGeometry->SetVisLevel(fEVisLevel->GetIntNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for visibility level.

void TGeoVolumeEditor::DoVisLevel()
{
   fBAuto->SetState(kButtonUp);
   fGeometry->SetVisLevel(fEVisLevel->GetIntNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for viewing volume and containers.

void TGeoVolumeEditor::DoViewAll()
{
   Bool_t on = (fBView[0]->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (!on) return;
   if (fVolume->IsVisContainers() == on) return;
   if (fVolume->IsRaytracing()) {
      fVolume->Raytrace(kFALSE);
      fBRaytrace->SetState(kButtonUp);
   }
   fVolume->SetVisContainers(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for viewing last leaves only.

void TGeoVolumeEditor::DoViewLeaves()
{
   Bool_t on = (fBView[1]->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (!on) return;
   if (fVolume->IsVisLeaves() == on) return;
   if (fVolume->IsRaytracing()) {
      fVolume->Raytrace(kFALSE);
      fBRaytrace->SetState(kButtonUp);
   }
   fVolume->SetVisLeaves(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for viewing volume only.

void TGeoVolumeEditor::DoViewOnly()
{
   Bool_t on = (fBView[2]->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (!on) return;
   if (fVolume->IsVisOnly() == on) return;
   if (fVolume->IsRaytracing()) {
      fVolume->Raytrace(kFALSE);
      fBRaytrace->SetState(kButtonUp);
   }
   fVolume->SetVisOnly(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for raytracing.

void TGeoVolumeEditor::DoRaytrace()
{
   Bool_t on = (fBRaytrace->GetState()==kButtonDown)?kTRUE:kFALSE;
   if (fVolume->IsRaytracing() == on) return;
   fVolume->Raytrace(on);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Modify division name.

void TGeoVolumeEditor::DoDivName()
{
   fApplyDiv->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Change division axis and preserve number of slices.

void TGeoVolumeEditor::DoDivSelAxis()
{
   Int_t iaxis = 1;
   for (Int_t i=0; i<3; i++) {
      if (fBDiv[i]->GetState()!=kButtonDown) continue;
      iaxis = i+1;
      break;
   }
   TGeoShape *shape = fVolume->GetShape();
   if (!shape) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   Double_t xlo, xhi;
   shape->GetAxisRange(iaxis, xlo, xhi);
   if (xhi <= xlo) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   fEDivFrom->SetNumber(xlo);
   fEDivStep->SetNumber(0);
   fApplyDiv->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle division range modification.

void TGeoVolumeEditor::DoDivFromTo()
{
   Double_t min, max, xlo, xhi, step;
   Int_t iaxis = 1;
   Int_t ndiv;
   for (Int_t i=0; i<3; i++) {
      if (fBDiv[i]->GetState()!=kButtonDown) continue;
      iaxis = i+1;
      break;
   }
   TGeoShape *shape = fVolume->GetShape();
   if (!shape) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   shape->GetAxisRange(iaxis, xlo, xhi);
   if (xhi-xlo <= 0) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   min = fEDivFrom->GetNumber();
   step = fEDivStep->GetNumber();
   ndiv = fEDivN->GetIntNumber();
   if (min<xlo) {
      min = xlo;
      fEDivFrom->SetNumber(xlo);
   }
   max = min + ndiv*step;
   if (max>xhi) {
      max = xhi;
      step = (max-min)/ndiv;
      fEDivStep->SetNumber(step);
   }
   if (min>=max) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   fApplyDiv->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle division step modification.

void TGeoVolumeEditor::DoDivStep()
{
   Double_t min, max, xlo, xhi;
   Int_t iaxis = 1;
   for (Int_t i=0; i<3; i++) {
      if (fBDiv[i]->GetState()!=kButtonDown) continue;
      iaxis = i+1;
      break;
   }
   TGeoShape *shape = fVolume->GetShape();
   if (!shape) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   shape->GetAxisRange(iaxis, xlo, xhi);
   if (xhi-xlo <= 0) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   min = fEDivFrom->GetNumber();
   Double_t step = fEDivStep->GetNumber();
   Int_t ndiv = fEDivN->GetIntNumber();
   max = min + ndiv*step;

   // Check if ndiv*step < max-min
   if (max <= xhi) {
      fApplyDiv->SetEnabled(kTRUE);
      return;
   }
   // Step too big - set value to fit range
   max = xhi;
   step = (max-min)/ndiv;
   fEDivStep->SetNumber(step);
   if (step < 0) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   fApplyDiv->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle division N modification.

void TGeoVolumeEditor::DoDivN()
{
   Double_t min, max, xlo, xhi;
   Int_t iaxis = 1;
   for (Int_t i=0; i<3; i++) {
      if (fBDiv[i]->GetState()!=kButtonDown) continue;
      iaxis = i+1;
      break;
   }
   TGeoShape *shape = fVolume->GetShape();
   if (!shape) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   shape->GetAxisRange(iaxis, xlo, xhi);
   if (xhi-xlo <= 0) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   Double_t step = fEDivStep->GetNumber();
   // If step=0 it is discounted
   if (step==0) {
      fApplyDiv->SetEnabled(kTRUE);
      return;
   }
   Int_t ndiv = fEDivN->GetIntNumber();
   min = fEDivFrom->GetNumber();
   max = min + ndiv*step;
   // Check if ndiv*step < max-min
   if (max <= xhi) {
      fApplyDiv->SetEnabled(kTRUE);
      return;
   }
   max = xhi;
   ndiv = (Int_t)((max-min)/step);
   fEDivN->SetNumber(ndiv);
   fApplyDiv->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Apply current division settings

void TGeoVolumeEditor::DoApplyDiv()
{
   Double_t xlo, xhi, step;
   Int_t iaxis = 1;
   Int_t ndiv;
   for (Int_t i=0; i<3; i++) {
      if (fBDiv[i]->GetState()!=kButtonDown) continue;
      iaxis = i+1;
      break;
   }
   TGeoShape *shape = fVolume->GetShape();
   if (!shape) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   shape->GetAxisRange(iaxis, xlo, xhi);
   if (xhi-xlo <= 0) {
      fApplyDiv->SetEnabled(kFALSE);
      return;
   }
   xlo = fEDivFrom->GetNumber();
   step = fEDivStep->GetNumber();
   ndiv = fEDivN->GetIntNumber();
   TGeoPatternFinder *finder = fVolume->GetFinder();
   if (finder) {
   // we have to remove first the existing division
      TObjArray *nodes = fVolume->GetNodes();
      nodes->Delete();
      nodes->Clear();
      delete finder;
      fVolume->SetFinder(0);
   }
   fVolume->Divide(fDivName->GetText(), iaxis, ndiv, xlo, step);
   fApplyDiv->SetEnabled(kFALSE);
   fGeometry->SetTopVisible();
   Update();
//   fVolume->Draw();
}
