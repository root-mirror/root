// @(#)root/ged:$Id$
// Author: Carsten Hof   16/08/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TH1Editor                                                           //
//  Editor for changing TH1 histogram attributes, rebinning & fitting.  //
//  For all possible draw options (there are a few which are not imple- //
//  mentable in graphical user interface) see THistPainter::Paint in    //
//  root/histpainter/THistPainter.cxx                                   //
//
//Begin_Html
/*
<img src="gif/TH1Editor_1.gif">
*/
//End_Html
//Begin_Html
/*
<img src="gif/TH1Editor_2.gif">
*/
//End_Html
//
//  These changes can be made via the TH1Editor:                        //
//    Style Tab:                                                        //
//      'Line'     : change Line attributes (color, thickness)          //
//                   see TAttLineEditor                                 //
//      'Fill'     : change Fill attributes (color, pattern)            //
//                   see TAttFillEditor                                 //
//      'Title'    : TextEntry: set the title of the histogram          //
//      'Histogram': change the draw options of the histogram           //
//          'Plot' : Radiobutton: draw a 2D or 3D plot of the histogram //
//                   according to the Plot dimension there will be      //
//                   different drawing possibilities (ComboBoxes/       //
//                   CheckBoxes)                                        //
//    2d Plot:                                                          //
//      'Error'   : ComboBox: add different error bars to the histogram //
//                  (no errors, simple, ...,  see THistPainter::Paint   //
//      'Add'     : ComboBox: further things which can be added to the  //
//                  histogram (None, simple/smooth line, fill area      //
//      'Simple Drawing': CheckBox: draw a simple histogram without     //
//                  errors (= "HIST" drawoption). In combination with   //
//                  some other draw options an outer line is drawn on   //
//                  top of the histogram                                //
//      'Show markers': CheckBox: draw a marker on to of each bin (="P" //
//                  drawoption)                                         //
//      'Draw bar chart': CheckBox: draw a bar chart (="B" drawoption)  //
//                  change the Fill Color with Fill in the Style Tab    //
//                  => will show Bar menue in the Style Tab             //
//      'Bar option': CheckBox: draw a bar chart (="BAR" drawoption)    //
//                  => will show Bar menue in the Style Tab             //
//    3d Plot:                                                          //
//      'Type'    : ComboBox: set histogram type Lego-Plot or Surface   //
//                  draw(Lego, Lego1.2, Surf, Surf1..5)                 //
//                  see THistPainter::Paint                             //
//      'Coords'  : ComboBox: set the coordinate system (Cartesian, ..  //
//                  Spheric) see THistPainter::Paint                    //
//      'Error'   : see 2D plot                                         //
//      'Bar'     : change the bar attributes                           //
//            'W' : change Bar Width                                    //
//            'O' : change Bar Offset                                   //
//      'Percentage': specifies the percentage of the bar which is drawn//
//                    brighter and darker (10% == BAR1 drawoption)      //
//      'Horizontal Bar': draw a horizontal bar chart                   //
//                                                                      //
//      'Marker'   : change the Marker attributes (color, appearance,   //
//                   thickness) see TAttMarkerEditor                    //
//Begin_Html
/*
<img src="gif/TH1Editor1.gif">
*/
//End_Html
//      This Tab has two different layouts. One is for a histogram which//
//      is not drawn from an ntuple. The other one is available for a   //
//      histogram which is drawn from an ntuple. In this case the rebin //
//      algorithm can create a rebinned histogram from the original data//
//      i.e. the ntuple.                                                //
//      To see te differences do:                                       //
//         TFile f("hsimple.root");                                     //
//         hpx->Draw("BAR1");        // non ntuple histogram            //
//         ntuple->Draw("px");       // ntuple histogram                //
//    Non ntuple histogram:                                             //
//       'Rebin': with the Slider the number of bins (shown in the field//
//                below the Slider) can be changed to any number which  //
//                divides the number of bins of the original histogram. //
//                Pushing 'Apply' will delete the origin histogram and  //
//                replace it by the rebinned one on the screen          //
//                Pushing 'Ignore' the origin histogram will be restored//
//    Histogram drawn from an ntuple:                                   //
//       'Rebin'  with the slider the number of bins can be enlarged by //
//                a factor of 2,3,4,5 (moving to the right) or reduced  //
//                by a factor of 1/2, 1/3, 1/4, 1/5                     //
//       'BinOffset': with the BinOffset slider the origin of the       //
//                histogram can be changed within one binwidth          //
//                Using this slider the effect of binning the data into //
//                bins can be made visible => statistical fluctuations  //
//       'Axis Range': with the DoubleSlider it is possible to zoom into//
//                the specified axis range. It is also possible to set  //
//                the upper and lower limit in fields below the slider  //
//       'Delayed drawing': all the Binning sliders can set to delay    //
//                draw mode. Then the changes on the histogram are only //
//                updated, when the Slider is released. This should be  //
//                activated if the redrawing of the histogram is too    //
//                time consuming.                                       //
//////////////////////////////////////////////////////////////////////////
//
//Begin_Html
/*
<img src="gif/TH1Editor1_1.gif">
*/
//End_Html
//Begin_Html
/*
<img src="gif/TH1Editor1_2.gif">
*/
//End_Html


#include "TH1Editor.h"
#include "TH1.h"
#include "TGedEditor.h"
#include "TGComboBox.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TVirtualPad.h"
#include "TString.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"
#include "TG3DLine.h"
#include "TGDoubleSlider.h"
#include "TGSlider.h"
#include "TView.h"
#include "TCanvas.h"
#include "TTreePlayer.h"
#include "TSelectorDraw.h"
#include "TGMsgBox.h"
#include "TGTab.h"
#include "TROOT.h"
#include "TVirtualX.h"

#include <cstdlib>

ClassImp(TH1Editor);

enum ETH1Wid{
   kTH1_TITLE,
   kTYPE_HIST,  kTYPE_LEGO,  kTYPE_LEGO1, kTYPE_LEGO2,
   kTYPE_SURF,  kTYPE_SURF1, kTYPE_SURF2, kTYPE_SURF3, kTYPE_SURF4, kTYPE_SURF5,
   kCOORDS_CAR, kCOORDS_CYL, kCOORDS_POL, kCOORDS_PSR, kCOORDS_SPH,
   kERRORS_NO,  kERRORS_SIMPLE, kERRORS_EDGES,
   kERRORS_REC, kERRORS_FILL,   kERRORS_CONTOUR,
   kHIST_TYPE,  kCOORD_TYPE, kERROR_TYPE, kMARKER_ONOFF, kB_ONOFF,  kBAR_ONOFF,
   kADD_TYPE,   kADD_NONE,   kADD_SIMPLE, kADD_SMOOTH,   kADD_FILL,
   kADD_BAR,    kADD_LINE,
   kDIM_SIMPLE, kDIM_COMPLEX,
   kPERCENT_TYPE, kPER_0, kPER_10, kPER_20, kPER_30, kPER_40,
   kBAR_H,      kBAR_WIDTH, kBAR_OFFSET,
   kSLIDER_MAX, kSLIDER_MIN,
   kDELAYED_DRAWING,
   kBINSLIDER, kBINSLIDER1, kBINOFFSET
};


////////////////////////////////////////////////////////////////////////////////
/// Constructor of histogram attribute GUI.

TH1Editor::TH1Editor(const TGWindow *p,  Int_t width,
                     Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back),
     fHist(0),
     fSameOpt(kFALSE),
     fBin(0),
     fBinHist(0)
{
   // TextEntry for changing the title of the histogram
   MakeTitle("Title");
   fTitlePrec = 2;
   fTitle = new TGTextEntry(this, new TGTextBuffer(50), kTH1_TITLE);
   fTitle->Resize(135, fTitle->GetDefaultHeight());
   fTitle->SetToolTipText("Enter the histogram title string");
   AddFrame(fTitle, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   // Histogram draw options
   TGCompositeFrame *fHistLbl = new TGCompositeFrame(this, 145, 10,
                                                           kHorizontalFrame |
                                                           kLHintsExpandX   |
                                                           kFixedWidth      |
                                                           kOwnBackground);
   fHistLbl->AddFrame(new TGLabel(fHistLbl,"Histogram"),
                      new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   fHistLbl->AddFrame(new TGHorizontal3DLine(fHistLbl),
                      new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 0));
   AddFrame(fHistLbl, new TGLayoutHints(kLHintsTop,0,0,2,0));

   // TGButtonGroup to change: 2D plot <-> 3D plot
   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fDimGroup = new TGHButtonGroup(f2,"Plot");
   fDimGroup->SetRadioButtonExclusive();
   fDim = new TGRadioButton(fDimGroup,"2-D",kDIM_SIMPLE);
   fDim->SetToolTipText("A 2-d plot of the histogram is dawn");
   fDim0 = new TGRadioButton(fDimGroup,"3-D",kDIM_COMPLEX);
   fDim0->SetToolTipText("A 3-d plot of the histogram is dawn");
   fDimGroup->SetLayoutHints(fDimlh=new TGLayoutHints(kLHintsLeft ,-2,3,3,-7),fDim);
   fDimGroup->SetLayoutHints(fDim0lh=new TGLayoutHints(kLHintsLeft ,16,-1,3,-7),fDim0);
   fDimGroup->Show();
   fDimGroup->ChangeOptions(kFitWidth | kChildFrame | kHorizontalFrame);
   f2->AddFrame(fDimGroup, new TGLayoutHints(kLHintsTop, 4, 1, 0, 0));
   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 2, 8));

   // Set the type of histogram (Lego0..2, Surf0..5) for 3D plot
   f3 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f3, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   TGCompositeFrame *f3a = new TGCompositeFrame(f3, 40, 20);
   f3->AddFrame(f3a, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   TGLabel *fType = new TGLabel(f3a, "Add: ");
   f3a->AddFrame(fType, new TGLayoutHints(kLHintsLeft, 6, 1, 4, 4));
   TGLabel *fCoords = new TGLabel(f3a, "Coords:");
   f3a->AddFrame(fCoords, new TGLayoutHints(kLHintsLeft, 6, 1, 4, 1));

   TGCompositeFrame *f3b = new TGCompositeFrame(f3, 40, 20);
   f3->AddFrame(f3b, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   fTypeCombo = BuildHistTypeComboBox(f3b, kHIST_TYPE);
   f3b->AddFrame(fTypeCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 1));
   fTypeCombo->Resize(80, 20);
   fTypeCombo->Associate(this);
   //Set the coordinate system (Cartesian, Spheric, ...)
   fCoordsCombo = BuildHistCoordsComboBox(f3b, kCOORD_TYPE);
   f3b->AddFrame(fCoordsCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 1));
   fCoordsCombo->Resize(80, 20);
   fCoordsCombo->Associate(this);

   // Set the Error (No error, error1..5)
   TGCompositeFrame *f5 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f5, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   TGCompositeFrame *f5a = new TGCompositeFrame(f5, 40, 20);
   f5->AddFrame(f5a, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   TGLabel *fError = new TGLabel(f5a, "Error:");
   f5a->AddFrame(fError, new TGLayoutHints(kLHintsLeft, 6, 2, 4, 1));

   TGCompositeFrame *f5b = new TGCompositeFrame(f5, 40, 20);
   f5->AddFrame(f5b, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   fErrorCombo = BuildHistErrorComboBox(f5b, kERROR_TYPE);
   f5b->AddFrame(fErrorCombo, new TGLayoutHints(kLHintsLeft, 15, 1, 2, 1));
   fErrorCombo->Resize(80, 20);
   fErrorCombo->Associate(this);

   // Further draw options: Smooth/Simple Line, Fill Area for 2D plot
   f6 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f6, new TGLayoutHints(kLHintsTop, 1, 1, 0, 3));

   TGCompositeFrame *f6a = new TGCompositeFrame(f6, 40, 20);
   f6->AddFrame(f6a, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   TGLabel *fAddLabel = new TGLabel(f6a, "Style:");
   f6a->AddFrame(fAddLabel, new TGLayoutHints(kLHintsLeft, 6, 2, 4, 1));

   TGCompositeFrame *f6b = new TGCompositeFrame(f6, 40, 20);
   f6->AddFrame(f6b, new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   fAddCombo = BuildHistAddComboBox(f6b, kADD_TYPE);
   f6b->AddFrame(fAddCombo, new TGLayoutHints(kLHintsLeft, 15, 1, 2, 1));
   fAddCombo->Resize(80, 20);
   fAddCombo->Associate(this);

   // option related to HIST: some changes needed here!
   // because of inconsistencies
   f15 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   fAddSimple = new TGCheckButton(f15, "Simple Drawing", kADD_LINE);
   fAddSimple ->SetToolTipText("A simple histogram without errors is drawn (draw option: Hist)");
   f15->AddFrame(fAddSimple, new TGLayoutHints(kLHintsLeft, 6, 1, 1, 0));
   AddFrame(f15, new TGLayoutHints(kLHintsTop, 1, 1, 0, -1));

   // Show Marker Checkbox: draw marker (or not)
   f7 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   fAddMarker = new TGCheckButton(f7, "Show markers", kMARKER_ONOFF);
   fAddMarker ->SetToolTipText("Make marker visible/invisible");
   f7->AddFrame(fAddMarker, new TGLayoutHints(kLHintsLeft, 6, 1, 1, 0));
   AddFrame(f7, new TGLayoutHints(kLHintsTop, 1, 1, 2, 0));

   // Bar Chart Checkbox: draw with option B
   f8 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   fAddB = new TGCheckButton(f8, "Draw bar chart", kB_ONOFF);
   fAddB ->SetToolTipText("Draw a bar chart");
   f8->AddFrame(fAddB, new TGLayoutHints(kLHintsLeft, 6, 1, 1, 0));
   AddFrame(f8, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   // Bar CheckBox: draw with option BAR +option selected by
   // fPercentCombo (0..4) e.g. BAR2
   f9 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   fAddBar = new TGCheckButton(f9, "Bar option", kBAR_ONOFF);
   fAddBar ->SetToolTipText("Draw bar chart with bar-option");
   f9->AddFrame(fAddBar, new TGLayoutHints(kLHintsLeft, 6, 1, 1, 0));
   AddFrame(f9, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   // Bar Menu => appears when the BAR checkbox is set
   f10 = new TGCompositeFrame(this, 145, 10, kHorizontalFrame |
                                             kLHintsExpandX   |
                                             kFixedWidth      |
                                             kOwnBackground);
   f10->AddFrame(new TGLabel(f10,"Bar"),
                 new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f10->AddFrame(new TGHorizontal3DLine(f10),
                 new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   AddFrame(f10, new TGLayoutHints(kLHintsTop,0,0,6,4));

   // NumberEntry to change the Bar Width
   f11 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   TGLabel *fWidthLbl = new TGLabel(f11, "W:");
   f11->AddFrame(fWidthLbl, new TGLayoutHints(kLHintsLeft, 1, 3, 4, 1));
   fBarWidth = new TGNumberEntry(f11, 1.00, 6, kBAR_WIDTH,
                                 TGNumberFormat::kNESRealTwo,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0.01, 1.);
   fBarWidth->GetNumberEntry()->SetToolTipText("Set bin bar width");
   fBarWidth->Resize(45,20);
   f11->AddFrame(fBarWidth, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 1));

   // NumberEntry to change the Bar OFfset
   TGLabel *foffsetLbl = new TGLabel(f11, "O:");
   f11->AddFrame(foffsetLbl, new TGLayoutHints(kLHintsLeft, 6,3, 4, 1));
   fBarOffset = new TGNumberEntry(f11, 0.00, 5, kBAR_OFFSET,
                                  TGNumberFormat::kNESRealTwo,
                                  TGNumberFormat::kNEAAnyNumber,
                                  TGNumberFormat::kNELLimitMinMax, -1., 1.);
   fBarOffset->GetNumberEntry()->SetToolTipText("Set bin bar offset");
   fBarOffset->Resize(50,20);
   f11->AddFrame(fBarOffset, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 1));
   AddFrame(f11, new TGLayoutHints(kLHintsTop, 1, 1, 0, 4));

   // ComboBox which specifies the width of the Bar which should be drawn
   // in another color i.e. specifies the number in BAR option e.g. BAR2
   f12 = new TGCompositeFrame(this, 80, 20, kVerticalFrame);
   TGCompositeFrame *f13 = new TGCompositeFrame(f12, 80, 20, kHorizontalFrame);
   TGLabel *percentLabel = new TGLabel(f13, "Percentage:");
   f13->AddFrame(percentLabel, new TGLayoutHints(kLHintsLeft, 6, 1, 3, 1));
   fPercentCombo = BuildPercentComboBox(f13, kPERCENT_TYPE);
   fPercentCombo->Resize(51, 20);
   fPercentCombo->Associate(f13);
   f13->AddFrame(fPercentCombo, new TGLayoutHints(kLHintsLeft, 14, 1, 2, 1));
   f12->AddFrame(f13,new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

   // CHeckBox for horizontal drawing of the Histogram
   fMakeHBar = new TGCheckButton(f12, "Horizontal Bar", kBAR_H);
   fMakeHBar ->SetToolTipText("Draw a horizontal bar chart with hBar-Option");
   f12->AddFrame(fMakeHBar, new TGLayoutHints(kLHintsLeft, 6, 1, 3, 0));
   AddFrame(f12, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   CreateBinTab();

   // add itself in the least of cleanups to be notified when attached histogram is deleted
   gROOT->GetListOfCleanups()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create binning tab.

void TH1Editor::CreateBinTab()
{
   fBin = CreateEditorTabSubFrame("Binning");

   TGCompositeFrame *title1 = new TGCompositeFrame(fBin, 145, 10,
                                                         kHorizontalFrame |
                                                         kLHintsExpandX   |
                                                         kFixedWidth      |
                                                         kOwnBackground);
   title1->AddFrame(new TGLabel(title1, "Rebin"),
                    new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   title1->AddFrame(new TGHorizontal3DLine(title1),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   fBin->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

   // Widgets for rebinning a histogram which does NOT derive from a ntuple
   fBinCont = new TGCompositeFrame(fBin, 80, 20, kVerticalFrame);
   TGCompositeFrame *f18 = new TGCompositeFrame(fBinCont, 80, 20,
                                                          kHorizontalFrame);
   fBinSlider  = new TGHSlider(f18, 100, kSlider1 | kScaleBoth);
   fBinSlider->Resize(107,20);
   f18->AddFrame(fBinSlider, new TGLayoutHints(kLHintsLeft, 3,0,0,3));
   fBinCont->AddFrame(f18, new TGLayoutHints(kLHintsTop, 15, 7, 3, 5));

   TGCompositeFrame *f20 = new TGCompositeFrame(fBinCont, 80, 20,
                                                          kHorizontalFrame);
   TGLabel *binLabel1 = new TGLabel(f20, "# of Bins:");
   f20->AddFrame(binLabel1, new TGLayoutHints(kLHintsLeft, 7, 1, 2, 1));
   fBinNumberEntry = new TGNumberEntryField(f20, kBINSLIDER, 0.0,
                                            TGNumberFormat::kNESInteger);
   ((TGTextEntry*)fBinNumberEntry)->SetToolTipText("Set the number of bins in the rebinned histogram");
   fBinNumberEntry->Resize(57,20);
   f20->AddFrame(fBinNumberEntry, new TGLayoutHints(kLHintsRight, 21, 0, 0, 0));
   fBinCont->AddFrame(f20, new TGLayoutHints(kLHintsTop, 0, 7, 3, 4));

   // Text buttons to Apply or Delete the rebinned histogram
   TGCompositeFrame *f23 = new TGCompositeFrame(fBinCont, 118, 20,
                                                          kHorizontalFrame |
                                                          kFixedWidth);
   fApply = new TGTextButton(f23, " &Apply ");
   f23->AddFrame(fApply,
                 new TGLayoutHints(kLHintsExpandX | kLHintsLeft , 0, 3, 4, 4));
   fCancel = new TGTextButton(f23, " &Ignore ");
   f23->AddFrame(fCancel,
                 new TGLayoutHints(kLHintsExpandX | kLHintsLeft, 3, 0, 4, 4));
   fBinCont->AddFrame(f23, new TGLayoutHints(kLHintsTop, 20, 3, 3, 4));
   fBin->AddFrame(fBinCont,new TGLayoutHints(kLHintsTop| kLHintsExpandX));

   // Widgets for rebinning a histogram which derives from a ntuple
   fBinCont1 = new TGCompositeFrame(fBin, 80, 20, kVerticalFrame);
   TGCompositeFrame *f21 = new TGCompositeFrame(fBinCont1, 80, 20,
                                                           kHorizontalFrame);
   fBinSlider1  = new TGHSlider(f21, 100, kSlider1 | kScaleBoth);
   fBinSlider1->Resize(107,20);
   fBinSlider1->SetRange(1,9);
   fBinSlider1->SetScale(12);
   fBinSlider1->SetPosition(5);
   f21->AddFrame(fBinSlider1, new TGLayoutHints(kLHintsLeft, 3,0,0,3));
   fBinCont1->AddFrame(f21, new TGLayoutHints(kLHintsTop, 15, 7, 5, 0));

   //  Lettering of the Rebin Slider
   TGCompositeFrame *f24 = new TGCompositeFrame(fBinCont1, 80, 20,
                                                           kHorizontalFrame);
   TGLabel *l1 = new TGLabel(f24, "-5");
   f24->AddFrame(l1, new TGLayoutHints(kLHintsLeft, 18, 1, -1, 0));
   TGLabel *l2 = new TGLabel(f24, "-2");
   f24->AddFrame(l2, new TGLayoutHints(kLHintsLeft, 26, 2, -1, 0));
   TGLabel *l3 = new TGLabel(f24, "2");
   f24->AddFrame(l3, new TGLayoutHints(kLHintsLeft, 17, 2, -1, 0));
   TGLabel *l4 = new TGLabel(f24, "5");
   f24->AddFrame(l4, new TGLayoutHints(kLHintsLeft, 32, 3, -1, 0));
   fBinCont1->AddFrame(f24, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));

   TGCompositeFrame *f22 = new TGCompositeFrame(fBinCont1, 140, 20,
                                                kHorizontalFrame);
   TGLabel *binLabel2 = new TGLabel(f22, "# of Bins:");
   f22->AddFrame(binLabel2, new TGLayoutHints(kLHintsLeft, 7, 1, 4, 1));

   fBinNumberEntry1 = new TGNumberEntryField(f22, kBINSLIDER1, 0.0,
                                             TGNumberFormat::kNESInteger);
   ((TGTextEntry*)fBinNumberEntry1)->SetToolTipText("Set the number of bins in the rebinned histogram");
   fBinNumberEntry1->Resize(57,20);
   f22->AddFrame(fBinNumberEntry1, new TGLayoutHints(kLHintsLeft, 21, 0, 2, 0));
   fBinCont1->AddFrame(f22, new TGLayoutHints(kLHintsTop, 0, 7, 2, 4));

   TGCompositeFrame *f26 = new TGCompositeFrame(fBinCont1, 80, 20,
                                                kHorizontalFrame);
   TGLabel *offsetLbl = new TGLabel(f26, "BinOffset:");
   f26->AddFrame(offsetLbl, new TGLayoutHints(kLHintsLeft, 6, 1, 2, 1));
   fOffsetNumberEntry = new TGNumberEntryField(f26, kBINOFFSET, 0.0,
                                               TGNumberFormat::kNESRealFour,
                                               TGNumberFormat::kNEAAnyNumber,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0., 1.);
   ((TGTextEntry*)fOffsetNumberEntry)->SetToolTipText("Add an offset to the origin of the histogram");
   fOffsetNumberEntry->Resize(57,20);
   f26->AddFrame(fOffsetNumberEntry,
                 new TGLayoutHints(kLHintsRight, 21, 0, 0, 0));
   fBinCont1->AddFrame(f26, new TGLayoutHints(kLHintsTop, 0, 7, 3, 1));

   TGCompositeFrame *f25 = new TGCompositeFrame(fBinCont1, 80, 20,
                                                           kHorizontalFrame);
   fBinOffsetSld  = new TGHSlider(f25, 100, kSlider1 | kScaleBoth);
   fBinOffsetSld->Resize(107,20);
   f25->AddFrame(fBinOffsetSld, new TGLayoutHints(kLHintsLeft, 15,0,0,2));
   fBinCont1->AddFrame(f25, new TGLayoutHints(kLHintsTop, 3, 7, 3, 3));
   fBin->AddFrame(fBinCont1, new TGLayoutHints(kLHintsTop));

   // Sliders for axis range
   TGCompositeFrame *sldCont = new TGCompositeFrame(fBin, 80, 20,
                                                    kVerticalFrame);
   TGCompositeFrame *title2 = new TGCompositeFrame(sldCont, 145, 10,
                                                            kHorizontalFrame |
                                                            kLHintsExpandX   |
                                                            kFixedWidth      |
                                                            kOwnBackground);
   title2->AddFrame(new TGLabel(title2, "Axis Range"),
                    new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   title2->AddFrame(new TGHorizontal3DLine(title2),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   sldCont->AddFrame(title2, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

   TGCompositeFrame *f14 = new TGCompositeFrame(sldCont, 80, 20,
                                                         kHorizontalFrame);
   TGLabel *fSliderLbl = new TGLabel(f14,"x:");
   f14->AddFrame(fSliderLbl,
                 new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 4,4, 4, 1));
   fSlider = new TGDoubleHSlider(f14, 1, 2);
   fSlider->Resize(118,20);
   f14->AddFrame(fSlider, new TGLayoutHints(kLHintsLeft));
   sldCont->AddFrame(f14, new TGLayoutHints(kLHintsTop, 3, 7, 4, 1));

   TGCompositeFrame *f16 = new TGCompositeFrame(sldCont, 80, 20,
                                                         kHorizontalFrame);
   fSldMin = new TGNumberEntryField(f16, kSLIDER_MIN, 0.0,
                                    TGNumberFormat::kNESRealTwo,
                                    TGNumberFormat::kNEAAnyNumber);
   ((TGTextEntry*)fSldMin)->SetToolTipText("Set the minimum value of the x-axis");
   fSldMin->Resize(57,20);
   f16->AddFrame(fSldMin, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
   fSldMax = new TGNumberEntryField(f16, kSLIDER_MAX, 0.0,
                                    TGNumberFormat::kNESRealTwo,
                                    TGNumberFormat::kNEAAnyNumber);
   ((TGTextEntry*)fSldMax)->SetToolTipText("Set the maximum value of the x-axis");
   fSldMax->Resize(57,20);
   f16->AddFrame(fSldMax, new TGLayoutHints(kLHintsLeft, 4, 0, 0, 0));
   sldCont->AddFrame(f16, new TGLayoutHints(kLHintsTop, 20, 3, 5, 0));

   TGCompositeFrame *f17 = new TGCompositeFrame(sldCont, 80, 20, kVerticalFrame);
   fDelaydraw = new TGCheckButton(f17, "Delayed drawing", kDELAYED_DRAWING);
   fDelaydraw ->SetToolTipText("Draw the new histogram only when any Slider is released");
   f17->AddFrame(fDelaydraw, new TGLayoutHints(kLHintsLeft, 6, 1, 2, 0));
   sldCont->AddFrame(f17, new TGLayoutHints(kLHintsTop, 1, 1, 5, 0));
   fBin->AddFrame(sldCont, new TGLayoutHints(kLHintsTop));

   // to avoid jumping from DoAddBar to DoAddB and vice versa
   fMakeB=kTRUE;
   // to avoid calling SetDrawoption after every change
   fMake=kTRUE;

   fBinHist = 0; // used to save a copy of the histogram

   // (when not drawn from an ntuple)
   fBinOffsetSld->SetRange(0,100);
   fBinOffsetSld->SetPosition(0);
   fOffsetNumberEntry->SetNumber(0.0000);
   fCancel->SetState(kButtonDisabled);
   fApply->SetState(kButtonDisabled);

}  // end bin tab

////////////////////////////////////////////////////////////////////////////////
/// Destructor of TH1 editor.

TH1Editor::~TH1Editor()
{
   // remove itselef from the list of cleanups
   gROOT->GetListOfCleanups()->Remove(this);

   // children of TGButonGroup are not deleted
   delete fDim;
   delete fDim0;
   delete fDimlh;
   delete fDim0lh;

   if (fBinHist) delete fBinHist;
   fBinHist = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TH1Editor::ConnectSignals2Slots()
{
   //widgets for draw options
   fAddB->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddB(Bool_t)");
   fAddBar->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddBar(Bool_t)");
   fTitle->Connect("TextChanged(const char *)", "TH1Editor", this, "DoTitle(const char *)");
   fTypeCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
   fCoordsCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
   fErrorCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
   fAddCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
   fAddMarker->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddMarker(Bool_t)");
   fAddSimple->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddSimple(Bool_t)");

   //change 2D <-> 3D plot
   fDimGroup->Connect("Clicked(Int_t)","TH1Editor",this,"DoHistView()");

   // change Bar Width/Offset, the second connection is needed to have the ability to confirm the value also with enter
   fBarWidth->Connect("ValueSet(Long_t)", "TH1Editor", this, "DoBarWidth()");
   (fBarWidth->GetNumberEntry())->Connect("ReturnPressed()", "TH1Editor", this, "DoBarWidth()");
   fBarOffset->Connect("ValueSet(Long_t)", "TH1Editor", this, "DoBarOffset()");
   (fBarOffset->GetNumberEntry())->Connect("ReturnPressed()", "TH1Editor", this, "DoBarOffset()");
   fPercentCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoPercent()");
   fMakeHBar-> Connect("Toggled(Bool_t)","TH1Editor",this,"DoHBar(Bool_t))");

   // Connections for rebinning are created - i.e. slider is
   // connected to the slots that perform the rebinning in the
   // case of a histogram not derived from an ntuple.
   fBinSlider->Connect("PositionChanged(Int_t)","TH1Editor",this, "DoBinMoved(Int_t)");
   fBinSlider->Connect("Released()","TH1Editor",this, "DoBinReleased()");
   fBinSlider->Connect("Pressed()","TH1Editor",this, "DoBinPressed()");
   // numberEntry which shows/sets the actual number of bins
   fBinNumberEntry->Connect("ReturnPressed()", "TH1Editor", this, "DoBinLabel()");
   // Buttons to accept/reject the rebinned histogram
   fApply->Connect("Clicked()", "TH1Editor", this, "DoApply()");
   fCancel->Connect("Pressed()", "TH1Editor", this, "DoCancel()");
   // in case of a histogram which is derived from an ntuple these slots are used
   fBinSlider1->Connect("Released()","TH1Editor",this, "DoBinReleased1()");
   fBinSlider1->Connect("PositionChanged(Int_t)","TH1Editor",this, "DoBinMoved1()");
   fBinNumberEntry1->Connect("ReturnPressed()", "TH1Editor", this, "DoBinLabel1()");
   // slider/slots to change the offset of the histogram
   fBinOffsetSld->Connect("PositionChanged(Int_t)", "TH1Editor", this,"DoOffsetMoved(Int_t)");
   fBinOffsetSld->Connect("Released()", "TH1Editor", this, "DoOffsetReleased()");
   fBinOffsetSld->Connect("Pressed()", "TH1Editor", this, "DoOffsetPressed()");
   fOffsetNumberEntry->Connect("ReturnPressed()", "TH1Editor", this, "DoBinOffset()");
   // slider/slots to set the visible axisrange
   fSlider->Connect("PositionChanged()","TH1Editor", this,"DoSliderMoved()");
   fSlider->Connect("Pressed()","TH1Editor", this, "DoSliderPressed()");
   fSlider->Connect("Released()","TH1Editor", this, "DoSliderReleased()");
   fSldMin->Connect("ReturnPressed()", "TH1Editor", this, "DoAxisRange()");
   fSldMax->Connect("ReturnPressed()", "TH1Editor", this, "DoAxisRange()");
   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if object is able to configure with this editor.

Bool_t TH1Editor::AcceptModel(TObject* obj)
{
   if (obj == 0 || !obj->InheritsFrom(TH1::Class()) ||
       ((TH1*)obj)->GetDimension()!=1 ||
       ((TH1*)obj)->GetEntries() == 0
       /*|| obj->InheritsFrom("TH2")  || obj->InheritsFrom("TProfile")*/) {
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up current values of histogram attributes.

void TH1Editor::SetModel(TObject* obj)
{

   if (fBinHist && (obj != fHist)) {
      //we have probably moved to a different pad.
      //let's restore the original histogram
      if (fHist) {
         fHist->Reset();
         fHist->SetBins(fBinHist->GetXaxis()->GetNbins(),
                        fBinHist->GetXaxis()->GetXmin(),
                        fBinHist->GetXaxis()->GetXmax());
         fHist->Add(fBinHist);
      }
      // delete in anycase fBinHist also when fHist is zero (i.e when it has been deleted)
      delete fBinHist; fBinHist = 0;
   }

   fHist = (TH1*)obj;
   fAvoidSignal = kTRUE;

   const char *text = fHist->GetTitle();
   fTitle->SetText(text);

   fMake=kFALSE;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   Bool_t errorset = kFALSE;
   // if no draw option is specified: (default options)
   if (str.IsNull() || str=="" ) {
      fDimGroup->SetButton(kDIM_SIMPLE, kTRUE);
      fDimGroup->SetButton(kDIM_COMPLEX, kFALSE);
      HideFrame(f3);  // Hiding the histogram type combo box
      ShowFrame(f6);
      ShowFrame(f7);
      ShowFrame(f8);
      ShowFrame(f9);
      HideFrame(f10);
      HideFrame(f11);
      HideFrame(f12);
      ShowFrame(f15);
      fCoordsCombo->Select(kCOORDS_CAR);
      fErrorCombo->Select(kERRORS_NO);
      errorset=kTRUE;
      fAddCombo->Select(kADD_NONE);
      fAddMarker->SetState(kButtonUp);
      fAddB->SetState(kButtonUp);
      fAddBar->SetState(kButtonUp);
      fAddSimple->SetState(kButtonDisabled);
      ChangeErrorCombo(1);
   // in case of a 2D plot:
   } else if (!str.Contains("LEGO") && !str.Contains("SURF")){
      fDimGroup->SetButton(kDIM_SIMPLE,kTRUE);
      fDimGroup->SetButton(kDIM_COMPLEX,kFALSE);
      HideFrame(f3);  // Hiding the histogram type combo box
      ShowFrame(f7);
      ShowFrame(f8);
      ShowFrame(f9);
      ShowFrame(f15);
      fCoordsCombo->Select(kCOORDS_CAR);
      // initialising fAddCombo
      if (str.Contains("C")) {
         if (str.Contains("CYL")) {
            TString dum = str;
            dum.Remove(strstr(dum.Data(),"CYL")-dum.Data(),3);
            if (dum.Contains("C")) fAddCombo->Select(kADD_SMOOTH);
         } else fAddCombo->Select(kADD_SMOOTH);
      }
      else if (str.Contains("LF2")) fAddCombo->Select(kADD_FILL);
      else if (str.Contains("L")){
         TString dum = str;
         if (str.Contains("CYL")) {
            dum.Remove(strstr(dum.Data(),"CYL")-dum.Data(),3);
            if (dum.Contains("L")) fAddCombo->Select(kADD_SIMPLE);
         }
         if (str.Contains("POL")) {
            dum.Remove(strstr(dum.Data(),"POL")-dum.Data(),3);
            if (dum.Contains("L")) fAddCombo->Select(kADD_SIMPLE);
         } else fAddCombo->Select(kADD_SIMPLE);
      } else fAddCombo->Select(kADD_NONE);

      if (fAddCombo->GetSelected()!=kADD_NONE)
         fAddSimple->SetState(kButtonDisabled);
      else if (str.Contains("HIST")) {
         if (str=="HIST") fAddSimple->SetState(kButtonDisabled);
         else fAddSimple->SetState(kButtonDown);
      } else fAddSimple->SetState(kButtonUp);

      if (str.Contains("B")) {
         TString dum = str;
         if (str.Contains("BAR")) {
            fAddBar->SetState(kButtonDown);
            fAddB->SetState(kButtonDisabled);
            ShowFrame(f10);
            ShowFrame(f11);
            ShowFrame(f12);
         } else {
            fAddB->SetState(kButtonDown);
            fAddBar->SetState(kButtonDisabled);
            fAddSimple->SetState(kButtonDisabled);
            ShowFrame(f10);
            ShowFrame(f11);
            HideFrame(f12);
         }
      } else {
         fAddB->SetState(kButtonUp);
         fAddBar->SetState(kButtonUp);
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
      }
      if (str.Contains("P") ) {
         fAddMarker->SetState(kButtonDown);
         fAddSimple->SetState(kButtonDisabled);
      } else if (!str.Contains("BAR")) fAddMarker->SetState(kButtonUp);
      ChangeErrorCombo(1);

   // in case of a 3D plot
   } else if (str.Contains("LEGO") || str.Contains("SURF")){
      fDimGroup->SetButton(kDIM_COMPLEX,kTRUE);
      fDimGroup->SetButton(kDIM_SIMPLE,kFALSE);
      TGListBox* lb;
      ChangeErrorCombo(0);
      // set Coordinate ComboBox
      if (str.Contains("SURF")){
         // surf cannot be combined with spheric and cartesian coordinates
         // i.e. remove them from the combobox
         fCoordsCombo->RemoveEntry(kCOORDS_SPH);
         fCoordsCombo->RemoveEntry(kCOORDS_CAR);
         lb = fCoordsCombo->GetListBox();
         lb->Resize(lb->GetWidth(), 49);
      } else {
         // surf cannot be combined with spheric and cartesian coordinates
         // if surf was selected before here the removed items were added the combobox again
         if (((TGLBContainer*)((TGListBox*)fCoordsCombo->GetListBox())->GetContainer())->GetPos(kCOORDS_SPH)==-1)
            fCoordsCombo->AddEntry("Spheric", kCOORDS_SPH);
         if (((TGLBContainer*)((TGListBox*)fCoordsCombo->GetListBox())->GetContainer())->GetPos(kCOORDS_CAR)==-1) {
            fCoordsCombo->AddEntry("Cartesian", kCOORDS_CAR);
            lb = fCoordsCombo->GetListBox();
            lb->Resize(lb->GetWidth(), 83);
         }
      }
      // initialising the Type Combobox
      if (str.Contains("LEGO2")) fTypeCombo->Select(kTYPE_LEGO2);
      else if (str.Contains("LEGO1")) fTypeCombo->Select(kTYPE_LEGO1);
      else if (str.Contains("LEGO")) fTypeCombo->Select(kTYPE_LEGO);
      else if (str.Contains("SURF5")) fTypeCombo->Select(kTYPE_SURF5);
      else if (str.Contains("SURF4")) fTypeCombo->Select(kTYPE_SURF4);
      else if (str.Contains("SURF3")) fTypeCombo->Select(kTYPE_SURF3);
      else if (str.Contains("SURF2")) fTypeCombo->Select(kTYPE_SURF2);
      else if (str.Contains("SURF1")) fTypeCombo->Select(kTYPE_SURF1);
      else if (str.Contains("SURF")) fTypeCombo->Select(kTYPE_SURF);

      if (str.Contains("CYL")) fCoordsCombo->Select(kCOORDS_CYL);
      else if (str.Contains("POL")) fCoordsCombo->Select(kCOORDS_POL);
      else if (str.Contains("SPH")) fCoordsCombo->Select(kCOORDS_SPH);
      else if (str.Contains("PSR")) fCoordsCombo->Select(kCOORDS_PSR);
      else fCoordsCombo->Select(kCOORDS_CAR); //default

      HideFrame(f6);
      HideFrame(f7);
      HideFrame(f8);
      HideFrame(f9);
      HideFrame(f15);
      if (str.Contains("LEGO")) {
         ShowFrame(f10);
         ShowFrame(f11);
         HideFrame(f12);
      } else {
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
      }
      fAddMarker->SetState(kButtonDisabled);
      fAddB->SetState(kButtonDisabled);
   }

   if (!errorset) {
      if (str.Contains("E1")) fErrorCombo->Select(kERRORS_EDGES);
      else if (str.Contains("E2")) fErrorCombo->Select(kERRORS_REC);
      else if (str.Contains("E3")) fErrorCombo->Select(kERRORS_FILL);
      else if (str.Contains("E4")) fErrorCombo->Select(kERRORS_CONTOUR);
      else if (str.Contains("E")) {
         if (str.Contains("LEGO")) {
            TString dum=str;
            dum.Remove(strstr(dum.Data(),"LEGO")-dum.Data(),4);
            if (dum.Contains("E")) fErrorCombo->Select(kERRORS_SIMPLE);
         } else fErrorCombo->Select(kERRORS_SIMPLE);
      } else fErrorCombo->Select(kERRORS_NO); //default
   }

   if (fErrorCombo->GetSelected() != kERRORS_NO){
      HideFrame(f7);
      HideFrame(f8);
   }
   if (str.Contains("BAR") || ((fAddBar->GetState()==kButtonDown) &&
       (fDim->GetState()==kButtonDown))) {
      ShowFrame(f10);
      ShowFrame(f11);
      ShowFrame(f12);
      fBarWidth->SetNumber(fHist->GetBarWidth());
      fBarOffset->SetNumber(fHist->GetBarOffset());
      if (str.Contains("HBAR")) fMakeHBar->SetState(kButtonDown);
      else fMakeHBar->SetState(kButtonUp);

      if (str.Contains("BAR4")) fPercentCombo->Select(kPER_40);
      else if (str.Contains("BAR3")) fPercentCombo->Select(kPER_30);
      else if (str.Contains("BAR2")) fPercentCombo->Select(kPER_20);
      else if (str.Contains("BAR1")) fPercentCombo->Select(kPER_10);
      else fPercentCombo->Select(kPER_0);
   }

   Int_t nx = fHist -> GetXaxis() -> GetNbins();
   Int_t nxbinmin = fHist -> GetXaxis() -> GetFirst();
   Int_t nxbinmax = fHist -> GetXaxis() -> GetLast();

   if (fDelaydraw->GetState()!=kButtonDown) fDelaydraw->SetState(kButtonUp);

   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();

   // Check if histogram is from ntupla/tree or not.
   // If it is a standard histogram or a ntupla based histogram
   // show a different frame in case of rebinning (fBinCont) with sliders and bin number entries
   // connected to different methods.
   // For example  the entry field fBinNumberEntry is connected to
   // the method DoBinLabel in case of non-ntupla histograms which just call Th1::Rebin
   // In csae of a tree based histogram the entry field fBinNumberEntry1 is used which is connected to
   // TH1Editor::DoBinLabel1 which is re-filling the histograms with the cached values from the TTreePlayer.
   // Since the actual number of histogram entry can be larger than the cache size of the TTreePlayer
   // (see JIRA ROOT-5900 or http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=17107 )
   // the GUI frame based on a non-tupla histogram is used when the number of entries of the histogram is
   // not the same as the number of filled entries in the TTreePlayer object.

   if (!player || player->GetHistogram()!=fHist ||
       fHist->GetEntries() != player->GetNfill()) {

      Int_t n = 0;
      if (fBinHist) n = fBinHist->GetXaxis()->GetNbins();
      else n = nx;
      if (n < 1) n = 1;
      fBin->HideFrame(fBinCont1);
      fBin->ShowFrame(fBinCont);
      Int_t* div = Dividers(n);
      Int_t up = 0;
      if (div[0]-1 <= 1) up = 2;
      else up = div[0]-1;
      fBinSlider->SetRange(1,up);
      Int_t i = 1;
      if (fBinSlider->GetMaxPosition()==2 && fBinSlider->GetPosition()==2)
         fBinSlider->SetPosition(2);
      else {
         while ( div[i] != nx) i ++;
         fBinSlider->SetPosition(div[0] - i + 1);
      }
      fBinNumberEntry->SetLimits(TGNumberFormat::kNELLimitMinMax , 2, n);
      fBinNumberEntry->SetIntNumber(nx);
      delete [] div;
   }
   else if (player && fHist==player->GetHistogram() && fHist->GetEntries() == player->GetNfill()) {
      // in case of a ntupla/tree based histogram with number of entries not exceeding the TTreePlayer cache
      fBin->HideFrame(fBinCont);
      fBin->ShowFrame(fBinCont1);
      fBinSlider->SetRange(0,1);
      fBinSlider->SetPosition(0);
      fBinSlider1->SetPosition(5);
      fBinNumberEntry1->SetLimits(TGNumberFormat::kNELLimitMinMax , 2, 10000);
      fBinNumberEntry1->SetIntNumber(nxbinmax-nxbinmin+1);
   }

   fSlider->SetRange(1,nx);
   fSlider->SetPosition((Double_t)nxbinmin,(Double_t)nxbinmax);

   fSldMin->SetNumber(fHist->GetXaxis()->GetBinLowEdge(nxbinmin));
   fSldMax->SetNumber(fHist->GetXaxis()->GetBinUpEdge(nxbinmax));

   fOffsetNumberEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, 0,
                                 fHist->GetXaxis()->GetBinWidth(1));

   if (fInit) ConnectSignals2Slots();
   fMake=kTRUE;
   fGedEditor->GetTab()->SetEnabled(1, kTRUE);
   fAvoidSignal = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the histogram title setting.

void TH1Editor::DoTitle(const char *text)
{
   if (fAvoidSignal) return;
   fHist->SetTitle(text);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the show markers check box.

void TH1Editor::DoAddMarker(Bool_t on)
{
   if (fAvoidSignal) return;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   TString dum = str;

   if (dum.Contains("POL")) dum.Remove(strstr(dum.Data(),"POL")-dum.Data(),3);
   if (dum.Contains("SPH")) dum.Remove(strstr(dum.Data(),"SPH")-dum.Data(),3);
   if (dum.Contains("PSR")) dum.Remove(strstr(dum.Data(),"PSR")-dum.Data(),3);
   if (on) {
      if (!dum.Contains("P")) str += "P";
      fAddSimple->SetState(kButtonDisabled);
      if (str.Contains("HIST"))
         str.Remove(strstr(str.Data(),"HIST")-str.Data(),4);
   } else if (fAddMarker->GetState()==kButtonUp) {
      if (str.Contains("POL") || str.Contains("SPH")) {
         while (dum.Contains("P"))
            dum.Remove(strstr(dum.Data(),"P")-dum.Data(),1);
         if (str.Contains("POL")) str = dum + "POL";
         if (str.Contains("SPH")) str = dum + "SPH";
         if (str.Contains("PSR")) str = dum + "PSR";
      } else if (str.Contains("P")) str.Remove(str.First("P"),1);
      if ((str=="HIST") || (str=="") ||
          (fAddB->GetState()==kButtonDown) ||
          fAddCombo->GetSelected() != kADD_NONE)
         fAddSimple->SetState(kButtonDisabled);
      else if (str.Contains("HIST"))
         fAddSimple->SetState(kButtonDown);
      else
         fAddSimple->SetState(kButtonUp);
   }
   if (fMake) {
      if (fSameOpt) str += "SAME";
      SetDrawOption(str);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the bar Add check box.

void TH1Editor::DoAddB(Bool_t on)
{
   if (fAvoidSignal) return;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   if (fMakeB) {
      fMakeB=kFALSE;
      if (on) {
         if (!str.Contains("B")) str += "B";
         ShowFrame(f10);
         ShowFrame(f11);
         HideFrame(f12);
         fAddBar->SetState(kButtonDisabled);
         fAddSimple->SetState(kButtonDisabled);
         fBarOffset->SetNumber(fHist->GetBarOffset());
         fBarWidth->SetNumber(fHist->GetBarWidth());
      } else if (fAddB->GetState()==kButtonUp) {
         while (str.Contains("B"))
            str.Remove(str.First("B"),1);
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
         fAddBar->SetState(kButtonUp);
         if (fAddMarker->GetState()!=kButtonDown &&
             !(str=="" || str=="HIST" ||
             fAddCombo->GetSelected()!=kADD_NONE))
            fAddSimple->SetState(kButtonUp);
      }
      if (fSameOpt) str += "SAME";
      if (fMake) SetDrawOption(str);
      Update();

      fMakeB=kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the bar Add check box.

void TH1Editor::DoAddBar(Bool_t on)
{
   if (fAvoidSignal) return;
   Disconnect(fAddMarker);
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   if (fMakeB) {
      fMakeB=kFALSE;
      Int_t o = 0;
      if (str.Contains("HBAR")) o=1;
      if (str.Contains("BAR4"))
         str.Remove(strstr(str.Data(),"BAR4")-str.Data()-o,4+o);
      else if (str.Contains("BAR3"))
         str.Remove(strstr(str.Data(),"BAR3")-str.Data()-o,4+o);
      else if (str.Contains("BAR2"))
         str.Remove(strstr(str.Data(),"BAR2")-str.Data()-o,4+o);
      else if (str.Contains("BAR1"))
         str.Remove(strstr(str.Data(),"BAR1")-str.Data()-o,4+o);
      else if (str.Contains("BAR0"))
         str.Remove(strstr(str.Data(),"BAR0")-str.Data()-o,4+o);
      else if (str.Contains("BAR"))
         str.Remove(strstr(str.Data(),"BAR")-str.Data()-o,3+o);
      if (on) {
         if ((fAddMarker->GetState()==kButtonDown) &&
             (fErrorCombo->GetSelected()==kERRORS_NO) &&
             (fAddSimple->GetState()!=kButtonDisabled))
            fAddSimple->SetState(kButtonDisabled);
         else if ((fAddMarker->GetState()!=kButtonDown) &&
                  (fAddSimple->GetState()==kButtonDisabled)) {
            if (str.Contains("HIST"))
               fAddSimple->SetState(kButtonDown);
            else if (fAddCombo->GetSelected()!=kADD_NONE)
               fAddSimple->SetState(kButtonDisabled);
            else
               fAddSimple->SetState(kButtonUp);
         }
         switch (fPercentCombo->GetSelected()){
            case(-1): {
               str += "BAR";
               fPercentCombo->Select(kPER_0);
               break;
            }
            case(kPER_0): {
               str += "BAR";
               break;
            }
            case(kPER_10): {
               str += "BAR1";
               break;
            }
            case(kPER_20): {
               str += "BAR2";
               break;
            }
            case(kPER_30): {
               str += "BAR3";
               break;
            }
            case(kPER_40): {
               str += "BAR4";
               break;
            }
         }
         ShowFrame(f10);
         ShowFrame(f11);
         ShowFrame(f12);
         if (fMakeHBar->GetState()==kButtonDown)
            str.Insert(strstr(str.Data(),"BAR")-str.Data(),"H");
         fBarOffset->SetNumber(fHist->GetBarOffset());
         fBarWidth->SetNumber(fHist->GetBarWidth());
         fAddB->SetState(kButtonDisabled);
      } else if (fAddBar->GetState()==kButtonUp) {
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
         fAddB->SetState(kButtonUp);
         if (fAddMarker->GetState()==kButtonDisabled)
            fAddMarker->SetState(kButtonUp);
         if (str=="" || str=="HIST" || fAddCombo->GetSelected() != kADD_NONE ||
             ((fAddMarker->GetState() == kButtonDown) &&
               fErrorCombo->GetSelected() == kERRORS_NO) )
            fAddSimple->SetState(kButtonDisabled);
      }
      if (fSameOpt) str += "SAME";
      if (fMake) SetDrawOption(str);
      Update();
      ((TGMainFrame*)GetMainFrame())->Layout();
      fMakeB=kTRUE;
   }
   fAddMarker->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddMarker(Bool_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to fAddSimple check box for drawing a simple histogram
/// without errors (== HIST draw option) in combination with some other
/// draw options. It draws an additional line on the top of the bins.

void TH1Editor::DoAddSimple(Bool_t on)
{
   if (fAvoidSignal) return;
   Disconnect(fAddMarker);
   //   Bool_t make=kFALSE;
   fMake = kFALSE;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   if (on) {
      if (!str.Contains("HIST")) {
         str += "HIST";
         fAddMarker->SetState(kButtonDisabled);
         fMake=kTRUE;
      }
   } else if (fAddSimple->GetState()==kButtonUp) {
      if (str.Contains("HIST")) {
         str.Remove(strstr(str.Data(),"HIST")-str.Data(),4);
         fAddMarker->SetState(kButtonUp);
         fMake=kTRUE;
      }
   }
   if (fSameOpt) str += "SAME";
   if (fMake) SetDrawOption(str);
   fAddMarker->Connect("Toggled(Bool_t)", "TH1Editor", this, "DoAddMarker(Bool_t)");
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the 'Plot' button group.

void TH1Editor::DoHistView()
{
   if (gPad && gPad->GetVirtCanvas())
      gPad->GetVirtCanvas()->SetCursor(kWatch);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kWatch));

   if (fDim->GetState() == kButtonDown)
      DoHistSimple();
   else
      DoHistComplex();

   if (gPad && gPad->GetVirtCanvas())
      gPad->GetVirtCanvas()->SetCursor(kPointer);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kPointer));
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the 2D radio button.

void TH1Editor::DoHistSimple()
{
   if (fAvoidSignal) return;
   if (fDim->GetState()==kButtonDown){
      TString str ="";
      fMake=kFALSE;
      TGListBox* lb;
      HideFrame(f3);
      ShowFrame(f6);
      ShowFrame(f9);
      ShowFrame(f15);
      ChangeErrorCombo(1);
      if ((fAddBar->GetState() != kButtonDown ||
           fAddMarker->GetState()==kButtonDown ) &&
          (fErrorCombo->GetSelected()==kERRORS_NO))
         fAddSimple->SetState(kButtonDisabled);
      else if ((fAddSimple->GetState()==kButtonDisabled) &&
               (fAddMarker->GetState()!=kButtonDown))
         fAddSimple->SetState(kButtonUp);
      else if (fAddSimple->GetState()!=kButtonUp)
         fAddSimple->SetState(kButtonDown);
      if (fAddMarker->GetState()==kButtonDisabled &&
          fAddSimple->GetState()!=kButtonDown)
         fAddMarker->SetState(kButtonUp);

      if (fErrorCombo->GetSelected()==kERRORS_NO) {
         ShowFrame(f7);
         ShowFrame(f8);
      } else {
         HideFrame(f7);
         HideFrame(f8);
         if (fAddBar->GetState()==kButtonDisabled)
            fAddBar->SetState(kButtonUp);
      }

      if ((fAddB->GetState() == kButtonDisabled)) {
         if (fAddBar->GetState()==kButtonDown) {
            ShowFrame(f10);
            ShowFrame(f11);
            ShowFrame(f12);
         } else {
            HideFrame(f10);
            HideFrame(f11);
            HideFrame(f12);
         }
      }
      if (fAddBar->GetState() == kButtonDisabled){
         ShowFrame(f10);
         ShowFrame(f11);
         HideFrame(f12);
      }
      if ((fAddBar->GetState() == kButtonUp) &&
          (fAddB->GetState() == kButtonUp)) {
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
      }
      if (fAddCombo->GetSelected()== -1 )fAddCombo->Select(kADD_NONE);
      if (fErrorCombo->GetSelected()!=kERRORS_NO) {
         fAddCombo->RemoveEntries(kADD_SIMPLE,kADD_FILL);
         lb = fAddCombo->GetListBox();
         lb->Resize(lb->GetWidth(),19);
         Disconnect(fAddCombo);
         fAddCombo->Select(kADD_NONE);
         fAddCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
      } else {
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_SIMPLE)==-1)
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Simple Line", kADD_SIMPLE);
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_SMOOTH)==-1)
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Smooth Line", kADD_SMOOTH);
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_FILL)==-1) {
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Fill Area",kADD_FILL);
            lb = fAddCombo->GetListBox();
            lb->Resize(lb->GetWidth(),76);
         }
      }
      if (fAddSimple->GetState()==kButtonDown) str+="HIST";
      str += GetHistErrorLabel()+GetHistAddLabel();
      if (fSameOpt) str += "SAME";
      SetDrawOption(str);
      Update();
      //fGedEditor->GetTab()->Layout();
      ((TGMainFrame*)GetMainFrame())->Layout();
      fMake=kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the 3D radio button.

void TH1Editor::DoHistComplex()
{
   if (fAvoidSignal) return;
   if (fDim0->GetState()==kButtonDown) {
      TString str ="";
      fMake=kFALSE;
      ShowFrame(f3);
      HideFrame(f6);
      HideFrame(f7);
      HideFrame(f8);
      HideFrame(f9);
      HideFrame(f15);
      ChangeErrorCombo(0);
      if (fTypeCombo->GetSelected()==-1 && fCoordsCombo->GetSelected()==-1) {
         str = "LEGO"+GetHistErrorLabel();
         fTypeCombo->Select(kTYPE_LEGO);
         fCoordsCombo->Select(kCOORDS_CAR);
      } else if (fTypeCombo->GetSelected()==-1){
         str = "LEGO"+GetHistErrorLabel();
         fTypeCombo->Select(kTYPE_LEGO);
      } else if (fCoordsCombo->GetSelected()==-1) {
         str = GetHistTypeLabel()+GetHistErrorLabel();
         fCoordsCombo->Select(kCOORDS_CAR);
      } else {
         str = GetHistTypeLabel()+GetHistCoordsLabel()+GetHistErrorLabel();
      }
      if (str.Contains("LEGO")) {
         ShowFrame(f10);
         ShowFrame(f11);
         HideFrame(f12);
      } else {
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
      }
      if (fSameOpt) str += "SAME";
      SetDrawOption(str);
      Update();
      ((TGMainFrame*)GetMainFrame())->Layout();
      fGedEditor->GetTab()->Layout();
      fMake=kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the histogram type, the coordinate type, the error type
/// and the Add combo box.

void TH1Editor::DoHistChanges()
{
   if (fAvoidSignal) return;
   fMakeB= kFALSE;
   TGListBox* lb;
   if (GetHistTypeLabel().Contains("SURF")) {
      if (fCoordsCombo->GetSelected()==kCOORDS_CAR ||
          fCoordsCombo->GetSelected()==kCOORDS_SPH)
         fCoordsCombo->Select(kCOORDS_POL);
      fCoordsCombo->RemoveEntry(kCOORDS_SPH);
      fCoordsCombo->RemoveEntry(kCOORDS_CAR);
      lb = fCoordsCombo->GetListBox();
      lb->Resize(lb->GetWidth(), 49);
   } else {
      if (((TGLBContainer*)((TGListBox*)fCoordsCombo->GetListBox())->GetContainer())->GetPos(kCOORDS_SPH)==-1)
         ((TGListBox*)fCoordsCombo->GetListBox())->AddEntrySort("Spheric", kCOORDS_SPH);
      if (((TGLBContainer*)((TGListBox*)fCoordsCombo->GetListBox())->GetContainer())->GetPos(kCOORDS_CAR)==-1) {
         ((TGListBox*)fCoordsCombo->GetListBox())->AddEntrySort("Cartesian", kCOORDS_CAR);
         lb = fCoordsCombo->GetListBox();
         lb->Resize(lb->GetWidth(), 83);
      }
   }
   if (fDim->GetState()!=kButtonUp){
      if (fErrorCombo->GetSelected() != kERRORS_NO){
         HideFrame(f7);
         HideFrame(f8);
         ShowFrame(f9);
         fAddMarker->SetState(kButtonDisabled);
         fAddB->SetState(kButtonDisabled);
         if (fAddBar->GetState()==kButtonDisabled)
            fAddBar->SetState(kButtonUp);
         if (fAddSimple->GetState()==kButtonDisabled)
            fAddSimple->SetState(kButtonUp);
         fAddCombo->RemoveEntries(kADD_SIMPLE,kADD_FILL);
         lb = fAddCombo->GetListBox();
         lb->Resize(lb->GetWidth(),19);
         Disconnect(fAddCombo);
         fAddCombo->Select(kADD_NONE);
         fAddCombo->Connect("Selected(Int_t)", "TH1Editor", this, "DoHistChanges()");
         if (fAddBar->GetState()==kButtonDown) {
            ShowFrame(f10);
            ShowFrame(f11);
            ShowFrame(f12);
         } else {
            HideFrame(f10);
            HideFrame(f11);
            HideFrame(f12);
         }
      } else {
         Bool_t on = fMake;
         fMake=kFALSE;
         ShowFrame(f7);
         ShowFrame(f8);
         ShowFrame(f9);
         if (fAddMarker->GetState()==kButtonDisabled)
            fAddMarker->SetState(kButtonUp);
         if (fAddBar->GetState() != kButtonDown &&
             fAddB->GetState()==kButtonDisabled)
            fAddB->SetState(kButtonUp);
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_SIMPLE)==-1)
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Simple Line", kADD_SIMPLE);
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_SMOOTH)==-1)
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Smooth Line", kADD_SMOOTH);
         if (((TGLBContainer*)((TGListBox*)fAddCombo->GetListBox())->GetContainer())->GetPos(kADD_FILL)==-1) {
            ((TGListBox*)fAddCombo->GetListBox())->AddEntry("Fill Area",kADD_FILL);
            lb = fAddCombo->GetListBox();
            lb->Resize(lb->GetWidth(),76);
         }
         fMake=on;
      }
      if (fAddCombo->GetSelected()!=kADD_NONE) {
         fAddSimple->SetState(kButtonDisabled);
      } else {
         if (fAddMarker->GetState()==kButtonDown)
            fAddSimple->SetState(kButtonDisabled);
         else if (fAddSimple->GetState()==kButtonDisabled)
            fAddSimple->SetState(kButtonUp);
      }
   } else if (fDim0->GetState()==kButtonDown) {
      if (GetHistTypeLabel().Contains("LEGO")) {
         ShowFrame(f10);
         ShowFrame(f11);
         HideFrame(f12);
      } else {
         HideFrame(f10);
         HideFrame(f11);
         HideFrame(f12);
      }
   }
   if (fMake) {
      TString str = "";
      if (fDim->GetState()==kButtonDown)
         str = GetHistErrorLabel()+GetHistAddLabel();
      else if (fDim0->GetState()==kButtonDown)
         str = GetHistTypeLabel()+GetHistCoordsLabel()+GetHistErrorLabel();
      if (fAddSimple->GetState()==kButtonDown)
         str += "HIST";
      if (fSameOpt)
         str += "SAME";
      SetDrawOption(str);
      if (str=="" || str=="HIST") fAddSimple->SetState(kButtonDisabled);
      Update();
   }
   ((TGMainFrame*)GetMainFrame())->Layout();
   //   fGedEditor->GetTab()->Layout();
   fMakeB=kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Bar Width of the Bar Charts.

void TH1Editor::DoBarWidth()
{
   if (fAvoidSignal) return;
   fHist->SetBarWidth(fBarWidth->GetNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Bar Offset of the Bar Charts.

void TH1Editor::DoBarOffset()
{
   if (fAvoidSignal) return;
   Float_t f = fBarOffset->GetNumber();
   fHist->SetBarOffset(f);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the bar percentage settings.

void TH1Editor::DoPercent()
{
   if (fAvoidSignal) return;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   Int_t o = 0;
   if (str.Contains("HBAR")) o=1;
   if (str.Contains("BAR4"))
      str.Remove(strstr(str.Data(),"BAR4")-str.Data()-1,4+o);
   else if (str.Contains("BAR3"))
      str.Remove(strstr(str.Data(),"BAR3")-str.Data()-o,4+o);
   else if (str.Contains("BAR2"))
      str.Remove(strstr(str.Data(),"BAR2")-str.Data()-o,4+o);
   else if (str.Contains("BAR1"))
      str.Remove(strstr(str.Data(),"BAR1")-str.Data()-o,4+o);
   else if (str.Contains("BAR0"))
      str.Remove(strstr(str.Data(),"BAR0")-str.Data()-o,4+o);
   else if (str.Contains("BAR"))
      str.Remove(strstr(str.Data(),"BAR")-str.Data()-o,3+o);

   if (fMakeHBar->GetState()==kButtonDown) str+="H";
   switch (fPercentCombo->GetSelected()){
      case (kPER_0) :{ str += "BAR"; break;}
      case (kPER_10):{ str += "BAR1"; break;}
      case (kPER_20):{ str += "BAR2"; break;}
      case (kPER_30):{ str += "BAR3"; break;}
      case (kPER_40):{ str += "BAR4"; break;}
   }
   if (fSameOpt) str += "SAME";
   if (fMake) SetDrawOption(str);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Horizontal Bar check button.

void TH1Editor::DoHBar(Bool_t on)
{
   if (fAvoidSignal) return;
   TString str = GetDrawOption();
   str.ToUpper();
   if (str.Contains("SAME"))
      fSameOpt = kTRUE;
   else
      fSameOpt = kFALSE;
   if (on) {
      if (!str.Contains("HBAR"))
         str.Insert(strstr(str.Data(),"BAR")-str.Data(),"H");
   }
   else if (fMakeHBar->GetState()==kButtonUp) {
      if (str.Contains("HBAR"))
         str.Remove(strstr(str.Data(),"BAR")-str.Data()-1,1);
   }
   if (fSameOpt) str += "SAME";
   if (fMake) SetDrawOption(str);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-Slider for redrawing of the histogram
/// according to the new Slider range.

void TH1Editor::DoSliderMoved()
{
   if (fAvoidSignal) return;
   if (fGedEditor->GetPad()->GetCanvas())
      fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
   fGedEditor->GetPad()->cd();
   if (fDelaydraw->GetState()==kButtonDown && fDim->GetState()==kButtonDown) {
      static Int_t px1,py1,px2,py2;
      static Float_t ymin,ymax,xleft,xright;
      xleft = fHist->GetXaxis()->GetBinLowEdge((Int_t)((fSlider->GetMinPosition())+0.5));
      xright =  fHist->GetXaxis()->GetBinUpEdge((Int_t)((fSlider->GetMaxPosition())+0.5));
      ymin  = fGedEditor->GetPad()->GetUymin();
      ymax  = fGedEditor->GetPad()->GetUymax();
      px1   = fGedEditor->GetPad()->XtoAbsPixel(xleft);
      py1   = fGedEditor->GetPad()->YtoAbsPixel(ymin);
      px2   = fGedEditor->GetPad()->XtoAbsPixel(xright);
      py2   = fGedEditor->GetPad()->YtoAbsPixel(ymax);
      if (fGedEditor->GetPad()->GetCanvas())
         fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
      fGedEditor->GetPad()->SetLineWidth(1);
      fGedEditor->GetPad()->SetLineColor(2);
      fGedEditor->GetPad()->SetLineWidth(1);
      fGedEditor->GetPad()->SetLineColor(2);
      fGedEditor->GetPad()->cd();
      gVirtualX->DrawBox(fPx1old, fPy1old, fPx2old, fPy2old, TVirtualX::kHollow);
      gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
      fPx1old = px1;
      fPy1old = py1;
      fPx2old = px2 ;
      fPy2old = py2;
      gVirtualX->Update(0);
      fSldMin->SetNumber(xleft);
      fSldMax->SetNumber(xright);
   }  else  if (fDelaydraw->GetState() == kButtonDown &&
                fDim0->GetState() == kButtonDown &&
                fCoordsCombo->GetSelected() == kCOORDS_CAR) {
      static Float_t p1[3], p2[3], p3[3], p4[3], p5[3], p6[3], p7[3], p8[3];
      TView *fView = fGedEditor->GetPad()->GetView();
      if (!fView) return;
      Double_t *rmin = fView->GetRmin();
      if (!rmin) return;
      Double_t *rmax = fView->GetRmax();
      if (!rmax) return;
      p1[0] = p4[0] = p5[0] = p8[0] =
            fHist->GetXaxis()->GetBinLowEdge((Int_t)((fSlider->GetMinPosition())+0.5));
      p2[0] = p3[0] = p6[0] = p7[0] =
            fHist->GetXaxis()->GetBinUpEdge((Int_t)((fSlider->GetMaxPosition())+0.5));
      p1[1] = p2[1] = p3[1] = p4[1] = rmin[1];
      p5[1] = p6[1] = p7[1] = p8[1] = rmax[1];
      p1[2] = p2[2] = p5[2] = p6[2] = rmin[2];
      p3[2] = p4[2] = p7[2] = p8[2] = rmax[2];
      fGedEditor->GetPad()->SetLineWidth(1);
      fGedEditor->GetPad()->SetLineColor(2);
      PaintBox3D(fP2old, fP3old, fP7old, fP6old);
      PaintBox3D(fP1old, fP4old, fP8old, fP5old);
      PaintBox3D(p2, p3, p7, p6);
      PaintBox3D(p1, p4, p8, p5);
      for (Int_t i = 0; i<3; i++){
         fP1old[i] = p1[i];
         fP2old[i] = p2[i];
         fP3old[i] = p3[i];
         fP4old[i] = p4[i];
         fP5old[i] = p5[i];
         fP6old[i] = p6[i];
         fP7old[i] = p7[i];
         fP8old[i] = p8[i];
      }
      fSldMin->SetNumber(p1[0]);
      fSldMax->SetNumber(p2[0]);
   } else  if (fDelaydraw->GetState() == kButtonDown &&
               fDim0->GetState() == kButtonDown) {
      fSldMin->SetNumber(fHist->GetXaxis()->GetBinLowEdge((Int_t)((fSlider->GetMinPosition())+0.5)));
      fSldMax->SetNumber(fHist->GetXaxis()->GetBinUpEdge((Int_t)((fSlider->GetMaxPosition())+0.5)));
   } else {
      fHist->GetXaxis()->SetRange((Int_t)((fSlider->GetMinPosition())+0.5),
                                  (Int_t)((fSlider->GetMaxPosition())+0.5));
      fSldMin->SetNumber(fHist->GetXaxis()->GetBinLowEdge(fHist->GetXaxis()->GetFirst()));
      fSldMax->SetNumber(fHist->GetXaxis()->GetBinUpEdge(fHist->GetXaxis()->GetLast()));
      fClient->NeedRedraw(fSlider,kTRUE);
      Update();
   }
   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
   if (player && player->GetHistogram() == fHist) {
      Int_t last = fHist->GetXaxis()->GetLast();
      Int_t first = fHist->GetXaxis()->GetFirst();
      fBinNumberEntry1->SetIntNumber(last-first+1);
      // How to redraw the NumberEntry without calling Update??
      // Update kills the "virtual" painted box in Delayed draw mode
      fClient->NeedRedraw(fBinNumberEntry1,kTRUE);
      //      fGedEditor->GetTab()->Layout();
   }
   fClient->NeedRedraw(fSldMin,kTRUE);
   fClient->NeedRedraw(fSldMax,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-axis Range slider for initialising the
/// values of the slider movement.

void TH1Editor::DoSliderPressed()
{
   if (fAvoidSignal) return;
   if (fGedEditor->GetPad()->GetCanvas())
      fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
   fGedEditor->GetPad()->cd();
   static Float_t ymin,ymax,xleft,xright;
   Int_t sldmin = (Int_t)((fSlider->GetMinPosition())+0.5);
   Int_t sldmax = (Int_t)((fSlider->GetMaxPosition())+0.5);
   if (fDelaydraw->GetState() == kButtonDown &&
       fDim->GetState()==kButtonDown) {
      if (fGedEditor->GetPad()->GetCanvas())
         fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
      fGedEditor->GetPad()->SetLineWidth(1);
      fGedEditor->GetPad()->SetLineColor(2);
      xleft = fHist->GetXaxis()->GetBinLowEdge(sldmin);
      xright =  fHist->GetXaxis()->GetBinUpEdge(sldmax);
      ymin  = fGedEditor->GetPad()->GetUymin();
      ymax  = fGedEditor->GetPad()->GetUymax();
      fPx1old   = fGedEditor->GetPad()->XtoAbsPixel(xleft);
      fPy1old   = fGedEditor->GetPad()->YtoAbsPixel(ymin);
      fPx2old   = fGedEditor->GetPad()->XtoAbsPixel(xright);
      fPy2old   = fGedEditor->GetPad()->YtoAbsPixel(ymax);
      gVirtualX->DrawBox(fPx1old, fPy1old, fPx2old, fPy2old, TVirtualX::kHollow);
   } else if (fDelaydraw->GetState() == kButtonDown &&
              fDim0->GetState() == kButtonDown &&
              fCoordsCombo->GetSelected() == kCOORDS_CAR) {
      TView *fView = fGedEditor->GetPad()->GetView();
      if (!fView) return;
      Double_t *rmin = fView->GetRmin();
      if (!rmin) return;
      Double_t *rmax = fView->GetRmax();
      if (!rmax) return;
      fP1old[0] = fP4old[0] = fP5old[0] = fP8old[0] =
                  fHist->GetXaxis()->GetBinLowEdge(sldmin);
      fP2old[0] = fP3old[0] = fP6old[0] = fP7old[0] =
                  fHist->GetXaxis()->GetBinUpEdge(sldmax);
      fP1old[1] = fP2old[1] = fP3old[1] = fP4old[1] = rmin[1];
      fP5old[1] = fP6old[1] = fP7old[1] = fP8old[1] = rmax[1];
      fP1old[2] = fP2old[2] = fP5old[2] = fP6old[2] = rmin[2];
      fP3old[2] = fP4old[2] = fP7old[2] = fP8old[2] = rmax[2];
      if (fGedEditor->GetPad()->GetCanvas())
         fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
      fGedEditor->GetPad()->SetLineWidth(1);
      fGedEditor->GetPad()->SetLineColor(2);
      PaintBox3D(fP2old, fP3old, fP7old, fP6old);
      PaintBox3D(fP1old, fP4old, fP8old, fP5old);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the x-axis Range slider for finalizing the
/// values of the slider movement.

void TH1Editor::DoSliderReleased()
{
   if (fAvoidSignal) return;
   if (fDelaydraw->GetState()==kButtonDown) {
      fHist->GetXaxis()->SetRange((Int_t)((fSlider->GetMinPosition())+0.5),
                                  (Int_t)((fSlider->GetMaxPosition())+0.5));
      fSldMin->SetNumber(fHist->GetXaxis()->GetBinLowEdge(fHist->GetXaxis()->GetFirst()));
      fSldMax->SetNumber(fHist->GetXaxis()->GetBinUpEdge(fHist->GetXaxis()->GetLast()));
      Update();
   }
   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
   if (player) if (player->GetHistogram() == fHist) {
      Int_t last = fHist->GetXaxis()->GetLast();
      Int_t first = fHist->GetXaxis()->GetFirst();
      fBinNumberEntry1->SetIntNumber(last-first+1);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the number entry fields containing the Max/Min
/// value of the x-axis.

void TH1Editor::DoAxisRange()
{
   if (fAvoidSignal) return;
   Int_t nx = fHist->GetXaxis()->GetNbins();
   Double_t width = fHist->GetXaxis()->GetBinWidth(1);
   Double_t lowLimit = fHist->GetXaxis()->GetBinLowEdge(1);
   Double_t upLimit = fHist->GetXaxis()->GetBinUpEdge(nx);
   if ((fSldMin->GetNumber()+width/2) < (lowLimit))
      fSldMin->SetNumber(lowLimit);
   if ((fSldMax->GetNumber()-width/2) > (upLimit))
      fSldMax->SetNumber(upLimit);
// Set the histogram range and the axis range slider
   fHist->GetXaxis()->SetRangeUser(fSldMin->GetNumber()+width/2,
                                   fSldMax->GetNumber()-width/2);
   Int_t nxbinmin = fHist->GetXaxis()->GetFirst();
   Int_t nxbinmax = fHist->GetXaxis()->GetLast();
   fSlider->SetPosition((Double_t)(nxbinmin),(Double_t)(nxbinmax));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the rebin slider in case of a not ntuple histogram
/// Updates some other widgets which are related to the rebin slider.

void TH1Editor::DoBinReleased()
{
   // draw the rebinned histogram in case of delay draw mode
   if (fAvoidSignal) return;
   if (fDelaydraw->GetState()==kButtonDown){
      if (!fBinHist) {
         fBinHist = (TH1*)fHist->Clone("BinHist");
         // we will manage this histogram
         fBinHist->SetDirectory(0);
      }
      Int_t nx = fBinHist->GetXaxis()->GetNbins();
      Int_t numx = fBinSlider->GetPosition();
      Int_t* divx = Dividers(nx);
      if (divx[0]==2) fBinSlider->SetPosition(2);
      if (divx[0]==2) {
         delete [] divx;
         return;
      }
      // delete the histogram which is on the screen
      fGedEditor->GetPad()->cd();
      fHist->Reset();
      fHist->SetBins(nx,fBinHist->GetXaxis()->GetXmin(),
                     fBinHist->GetXaxis()->GetXmax());
      fHist->Add(fBinHist);
      fHist->SetCanExtend(TH1::kNoAxis);
      fHist->Rebin(divx[numx]);
      // fModel=fHist;
      if (divx[0]!=2) {
         TAxis* xaxis = fHist->GetXaxis();
         Double_t xBinWidth = xaxis->GetBinWidth(1);
         xaxis->SetRangeUser(fSldMin->GetNumber()+xBinWidth/2,
                             fSldMax->GetNumber()-xBinWidth/2);
         fSlider->SetRange(1,(Int_t)nx/divx[numx]);
         fSlider->SetPosition(xaxis->FindBin(fSldMin->GetNumber()+xBinWidth/2),
                              xaxis->FindBin(fSldMax->GetNumber()-xBinWidth/2));
         // the x-axis range could be changed a little bit by Rebin algorithm
         fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
         fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
      }
      if (fCancel->GetState()==kButtonDisabled)
         fCancel->SetState(kButtonUp);
      if (fApply->GetState()==kButtonDisabled)
         fApply->SetState(kButtonUp);
      Update();
      delete [] divx;
   }
//   fGedEditor->GetPad()->GetCanvas()->Selected(fGedEditor->GetPad(), fHist,  0);
   //  fModel = fHist;
      Refresh(fHist);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the rebin slider in case of a not ntuple histogram
/// (does the Rebinning of the histogram).

void TH1Editor::DoBinMoved(Int_t numx)
{
   // create a clone in the background, when the slider is moved for
   // the first time
   if (fAvoidSignal) return;
   if (!fBinHist /*&& fDelaydraw->GetState()!=kButtonDown*/) {
      Int_t* divx = Dividers(fHist->GetXaxis()->GetNbins());
      if (divx[0]==2) {
         delete [] divx;
         return;
      }
      fBinHist = (TH1*)fHist->Clone("BinHist");
      // the TH1Editor class  manage this histogram
      fBinHist->SetDirectory(0);
      delete [] divx;
   }
   // if the slider already has been moved and the clone is saved
   Int_t nx = fBinHist->GetXaxis()->GetNbins();
   Int_t* divx = Dividers(nx);
   if (divx[0]==2) {
      fBinSlider->SetPosition(2);
      numx=1;
      delete [] divx;
      return;
   }
   Int_t maxx = (Int_t)nx/divx[numx];
   if (maxx==1) maxx=2;
   if (fDelaydraw->GetState() == kButtonUp) {
      fGedEditor->GetPad()->cd();
      fHist->Reset();
      fHist->SetBins(nx,fBinHist->GetXaxis()->GetXmin(),
                     fBinHist->GetXaxis()->GetXmax());
      fHist->Add(fBinHist);
      fHist->SetCanExtend(TH1::kNoAxis);
      fHist->Rebin(divx[numx]);
      //fModel=fHist;
      TAxis* xaxis = fHist->GetXaxis();
      Double_t xBinWidth = xaxis->GetBinWidth(1);
      xaxis->SetRangeUser(fSldMin->GetNumber()+xBinWidth/2,
                          fSldMax->GetNumber()-xBinWidth/2);
      fSlider->SetRange(1,maxx);
      fSlider->SetPosition(xaxis->FindBin(fSldMin->GetNumber()+xBinWidth/2),
                           xaxis->FindBin(fSldMax->GetNumber()-xBinWidth/2));
      // the axis range could be changed a little bit by the Rebin algorithm
      fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
      fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
      fClient->NeedRedraw(fBinSlider,kTRUE);
      Update();
   }
   if (fCancel->GetState()==kButtonDisabled)
      fCancel->SetState(kButtonUp);
   if (fApply->GetState()==kButtonDisabled)
      fApply->SetState(kButtonUp);
   if (fBinNumberEntry->GetNumber()!=maxx)
      fBinNumberEntry->SetNumber(maxx);
   delete [] divx;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the rebin slider in case of a not ntuple histogram.

void TH1Editor::DoBinPressed()
{
   if (fAvoidSignal) return;
   Int_t* d = Dividers(fHist->GetXaxis()->GetNbins());
   if (d[0]==2 && !fBinHist) {
      new TGMsgBox(fClient->GetDefaultRoot(), this->GetMainFrame(),
                   "TH1 Editor", "It is not possible to rebin the histogram",
                   kMBIconExclamation, kMBOk, 0, kVerticalFrame);
      gVirtualX->GrabPointer(fBinSlider->GetId(),0,0,0);
   }
   delete [] d;
   // calling the MessageBox again does NOT work!*/
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the BinNumber Slider in case of a ntuple histogram
/// (does the Rebinning of the histogram).

void TH1Editor::DoBinReleased1()
{
   if (fAvoidSignal) return;
   Double_t oldOffset = fOffsetNumberEntry->GetNumber();
   Int_t number = fBinSlider1->GetPosition();
   if (number==5) return;
   Int_t fact = 0;
   Int_t binNumber = 0;
   TAxis* xaxis = fHist->GetXaxis();
   // "compute" the scaling factor:
   if (number > 5) fact = number - 4;
   else fact = number - 6;
   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
   if (!player) return;
   Int_t first = xaxis->GetFirst();
   Int_t last = xaxis->GetLast();
   Int_t nx = xaxis->GetNbins();
   Double_t min = xaxis->GetBinLowEdge(1);       // overall min in user coords
   Double_t max = xaxis->GetBinUpEdge(nx);       // overall max in user coords
   Double_t rmin = xaxis->GetBinLowEdge(first);  // recent min in user coords
   Double_t rmax = xaxis->GetBinUpEdge(last);    // recent max in user coords

   ((TH1*)player->GetHistogram())->SetCanExtend(TH1::kNoAxis);
   ((TH1*)player->GetHistogram())->Reset();

   // get new Number of bins
   if (fact > 0) binNumber = fact*nx;
   if (fact < 0) binNumber = (Int_t) ((-1)*nx/fact+0.5);
   if (binNumber < 1) binNumber = 1;
   if (binNumber > 10000) binNumber= 10000;
   Double_t newOffset = 1.*fBinOffsetSld->GetPosition()/100*((max-min)/binNumber);
   // create new histogram - the main job is done by sel->TakeAction()
   ((TH1*)player->GetHistogram())->SetBins(binNumber,
                                           min-oldOffset+newOffset,
                                           max-oldOffset+newOffset);
   TSelectorDraw *sel = (TSelectorDraw*)player->GetSelector();
   if (!sel) return;
   sel->TakeAction();

   // restore and set all the attributes which were changed by TakeAction()
   fHist = (TH1*)((TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer())->GetHistogram();
   fSlider->SetRange(1,binNumber);
   Double_t binWidth = fHist->GetXaxis()->GetBinWidth(1);
   fSlider->SetPosition(xaxis->FindBin(rmin), xaxis->FindBin(rmax));
   Double_t offset = 1.*fBinOffsetSld->GetPosition()/100*binWidth;
   xaxis->SetRange(xaxis->FindBin(rmin+binWidth/2),
                   xaxis->FindBin(rmax-binWidth/2));   // SetRange in binNumbers!
   fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
   fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
   fBinNumberEntry1->SetNumber(xaxis->GetLast() - xaxis->GetFirst() + 1);
   fBinSlider1->SetPosition(5);
   fOffsetNumberEntry->SetNumber(offset);
   fOffsetNumberEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, 0,
                                 xaxis->GetBinWidth(1));
   fClient->NeedRedraw(fBinSlider1);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the rebin slider in case of an ntuple histogram.
/// It updates the BinNumberEntryField during the BinSlider movement.

void TH1Editor::DoBinMoved1()
{
   if (fAvoidSignal) return;
   TAxis* xaxis = fHist->GetXaxis();
   Int_t first = xaxis->GetFirst();
   Int_t last = xaxis->GetLast();
   Int_t number = fBinSlider1->GetPosition();
   Int_t n = last -first+1;
   Int_t fact = 0;
   Int_t binNumber = 0;
   if (number >= 5) fact = number - 4;
   else fact = number - 6;
   if (fact > 0) binNumber = fact*n;
   if (fact < 0) binNumber = (Int_t) ((-1)*n/fact+0.5);
   if (binNumber < 1) binNumber = 1;
   if (binNumber > 10000) binNumber= 10000;
   fBinNumberEntry1->SetIntNumber(binNumber);
//   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Bin number entry of the Rebinning tab.

void TH1Editor::DoBinLabel()
{
   if (fAvoidSignal) return;
   Int_t num = (Int_t)(fBinNumberEntry->GetNumber());
   Int_t nx = 0;
   if (fBinHist) nx = fBinHist->GetXaxis()->GetNbins();
   else nx = fHist->GetXaxis()->GetNbins();
   if (nx < 2) return;
   Int_t *div = Dividers(nx);
   Int_t diff = TMath::Abs(num - div[1]);
   Int_t c = 1;
   for (Int_t i = 2; i <= div[0]; i++) {
      if ((TMath::Abs(num - div[i])) < diff) {
         c = i;
         diff = TMath::Abs(num - div[i]);
      }
   }
   fBinNumberEntry->SetNumber(div[c]);
   fBinSlider->SetPosition(div[0] - c +1);
   if (fDelaydraw->GetState()==kButtonUp) DoBinMoved(div[0] - c +1);
   else DoBinReleased();
//   fGedEditor->GetPad()->GetCanvas()->Selected(fGedEditor->GetPad(), fHist,  0);
   // fModel = fHist;
   Refresh(fHist);
   delete [] div;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Bin number entry of the Rebinning tab.

void TH1Editor::DoBinLabel1()
{
   if (fAvoidSignal) return;
   Double_t oldOffset = fOffsetNumberEntry->GetNumber();
   Int_t num = (Int_t)fBinNumberEntry1->GetNumber();
   TAxis* xaxis = fHist->GetXaxis();
   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
   if (!player) return;
   Int_t first = xaxis->GetFirst();
   Int_t last = xaxis->GetLast();
   Int_t nx = xaxis->GetNbins();
   Double_t min = xaxis->GetBinLowEdge(1);        // overall min in user coords
   Double_t max = xaxis->GetBinUpEdge(nx);        // overall max in user coords
   Double_t rmin = xaxis->GetBinLowEdge(first);   // recent min in user coords
   Double_t rmax = xaxis->GetBinUpEdge(last);     // recent max in user coords

   ((TH1*)player->GetHistogram())->SetCanExtend(TH1::kNoAxis);
   ((TH1*)player->GetHistogram())->Reset();

// Calculate the new number of bins in the complete range
   Int_t binNumber = (Int_t) ((max-min)/(rmax - rmin)*num + 0.5);
   if (binNumber < 1) binNumber = 1;
   if (binNumber > 10000) binNumber = 10000;
   Double_t offset = 1.*(fBinOffsetSld->GetPosition())/100*(max-min)/binNumber;
// create new histogram - the main job is done by sel->TakeAction()
   ((TH1*)player->GetHistogram())->SetBins(binNumber,
                                           min-oldOffset+offset,
                                           max-oldOffset+offset);
   TSelectorDraw *sel = (TSelectorDraw*)player->GetSelector();
   if (!sel) return;
   sel->TakeAction();

// Restore and set all the attributes which were changed by TakeAction()
   fHist = (TH1*)((TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer())->GetHistogram();
   fSlider->SetRange(1,binNumber);
   Double_t binWidth = xaxis->GetBinWidth(1);
   fSlider->SetPosition(xaxis->FindBin(rmin), xaxis->FindBin(rmax));
   offset = 1.*fBinOffsetSld->GetPosition()/100*binWidth;
   xaxis->SetRange(xaxis->FindBin(rmin+binWidth/2),
                   xaxis->FindBin(rmax-binWidth/2));   // SetRange in binNumbers!
   fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
   fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
   fOffsetNumberEntry->SetNumber(offset);
   fOffsetNumberEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, 0, binWidth);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the OffSetSlider that saves the OldBinOffset
/// (nessesary for delay draw mode).

void TH1Editor::DoOffsetPressed()
{
   if (fAvoidSignal) return;
   fOldOffset = fOffsetNumberEntry->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the OffSetSlider.
/// It changes the origin of the histogram inbetween a binwidth and
/// rebin the histogram with the new Offset given by the Slider.

void TH1Editor::DoOffsetReleased()
{
   // !!problem: histogram with variable binwidth??
   // computes the new histogram in "delay draw" mode

   if (fAvoidSignal) return;
   if (fDelaydraw->GetState()==kButtonDown) {
      Int_t num = (Int_t) fBinOffsetSld->GetPosition();
      TAxis* xaxis = fHist->GetXaxis();
      Double_t binWidth = xaxis->GetBinWidth(1);
      Double_t offset =  1.*num/100*binWidth;
      Double_t oldOffset = fOldOffset;
      Int_t nx = xaxis->GetNbins();
      TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
      if (!player) return;
      Int_t first = xaxis->GetFirst();
      Int_t last = xaxis->GetLast();
      Double_t min = xaxis->GetBinLowEdge(1);          // overall min in user coords
      Double_t max = xaxis->GetBinUpEdge(nx);          // overall max in user coords
      Double_t rmin = xaxis->GetBinLowEdge(first);     // recent min in user coords
      Double_t rmax = xaxis->GetBinUpEdge(last);       // recent max in user coords

      ((TH1*)player->GetHistogram())->SetCanExtend(TH1::kNoAxis);
      ((TH1*)player->GetHistogram())->Reset();

      ((TH1*)player->GetHistogram())->SetBins(nx,
                                              min+offset-oldOffset,
                                              max+offset-oldOffset);
      TSelectorDraw *sel = (TSelectorDraw*)player->GetSelector();
      if (!sel) return;
      sel->TakeAction();

      // Restore all the attributes which were changed by TakeAction()
      fHist = (TH1*)((TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer())->GetHistogram();
      xaxis->SetRange(xaxis->FindBin(rmin+offset-oldOffset+binWidth/2),
                      xaxis->FindBin(rmax+offset-oldOffset-binWidth/2)); // in binNumbers!
      fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
      fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
      fOffsetNumberEntry->SetNumber(offset);
      Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the OffSetSlider.
/// It changes the origin of the histogram inbetween a binwidth and
/// rebin the histogram with the new offset given by the Slider.

void TH1Editor::DoOffsetMoved(Int_t num)
{
   // !!histogram with variable binwidth??
   // !!only works for histograms with fixed binwidth

   if (fAvoidSignal) return;
   TAxis* xaxis = fHist->GetXaxis();
   Double_t binWidth = xaxis->GetBinWidth(1);
   Double_t offset =  1.*num/100*binWidth;
   if (fDelaydraw->GetState()==kButtonUp) {
      Double_t oldOffset = fOffsetNumberEntry->GetNumber();
      Int_t nx = xaxis->GetNbins();
      TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
      if (!player) return;
      Int_t first = xaxis->GetFirst();
      Int_t last = xaxis->GetLast();
      Double_t min = xaxis->GetBinLowEdge(1);          // overall min in user coords
      Double_t max = xaxis->GetBinUpEdge(nx);          // overall max in user coords
      Double_t rmin = xaxis->GetBinLowEdge(first);     // recent min in user coords
      Double_t rmax = xaxis->GetBinUpEdge(last);       // recent max in user coords

      ((TH1*)player->GetHistogram())->SetCanExtend(TH1::kNoAxis);
      ((TH1*)player->GetHistogram())->Reset();

      ((TH1*)player->GetHistogram())->SetBins(nx,
                                              min+offset-oldOffset,
                                              max+offset-oldOffset);
      TSelectorDraw *sel = (TSelectorDraw*)player->GetSelector();
      if (!sel) return;
      sel->TakeAction();

   // Restore all the attributes which were changed by TakeAction()
      fHist = (TH1*)((TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer())->GetHistogram();
      xaxis->SetRange(xaxis->FindBin(rmin+offset-oldOffset+binWidth/2),
                      xaxis->FindBin(rmax+offset-oldOffset-binWidth/2)); // in binNumbers!
      fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
      fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
      fClient->NeedRedraw(fBinOffsetSld,kTRUE);
   }
   fOffsetNumberEntry->SetNumber(offset);
   fClient->NeedRedraw(fOffsetNumberEntry,kTRUE);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the OffSetNumberEntry which is related to the
/// OffSetSlider changes the origin of the histogram inbetween a binwidth.

void TH1Editor::DoBinOffset()
{
   if (fAvoidSignal) return;
   TAxis* xaxis = fHist->GetXaxis();
   Double_t binWidth = xaxis->GetBinWidth(1);
   Double_t offset =  fOffsetNumberEntry->GetNumber();
   Double_t oldOffset = 1.*fBinOffsetSld->GetPosition()/100*binWidth;
   Int_t nx = xaxis->GetNbins();
   TTreePlayer *player = (TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer();
   if (!player) return;
   Int_t first = xaxis->GetFirst();
   Int_t last = xaxis->GetLast();
   Double_t min = xaxis->GetBinLowEdge(1);          // overall min in user coords
   Double_t max = xaxis->GetBinUpEdge(nx);          // overall max in user coords
   Double_t rmin = xaxis->GetBinLowEdge(first);     // recent min in user coords
   Double_t rmax = xaxis->GetBinUpEdge(last);       // recent max in user coords

   ((TH1*)player->GetHistogram())->SetCanExtend(TH1::kNoAxis);
   ((TH1*)player->GetHistogram())->Reset();

   ((TH1*)player->GetHistogram())->SetBins(nx,
                                           min+offset-oldOffset,
                                           max+offset-oldOffset);
   TSelectorDraw *sel = (TSelectorDraw*)player->GetSelector();
   if (!sel) return;
   sel->TakeAction();

    // Restore all the attributes which were changed by TakeAction()
   fHist = (TH1*)((TTreePlayer*)TVirtualTreePlayer::GetCurrentPlayer())->GetHistogram();
   xaxis->SetRange(xaxis->FindBin(rmin+offset-oldOffset+binWidth/2),
                   xaxis->FindBin(rmax+offset-oldOffset-binWidth/2)); // in binNumbers!
   fSldMin->SetNumber(xaxis->GetBinLowEdge(xaxis->GetFirst()));
   fSldMax->SetNumber(xaxis->GetBinUpEdge(xaxis->GetLast()));
   fBinOffsetSld->SetPosition((Int_t)(offset/binWidth*100));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Apply button of the Binning tab.

void TH1Editor::DoApply()
{
   Int_t ret = 0;
   new TGMsgBox(fClient->GetDefaultRoot(), this->GetMainFrame(),
                "TH1 Editor", "Replace origin histogram with rebinned one?",
                kMBIconQuestion, kMBYes | kMBNo, &ret, kVerticalFrame);
   if (ret==1) {
      if (fBinHist) {
         delete fBinHist;
         fBinHist = 0;
      }
      Int_t nx = fHist->GetXaxis()->GetNbins();
      Int_t *div = Dividers(nx);
      Int_t up = 0;
      if (div[0]-1 <= 1) up = 2;
      else up = div[0]-1;
      fBinSlider->SetRange(1,up);
      if (fBinSlider->GetMaxPosition()==2 && div[0]==2 )
         fBinSlider->SetPosition(2);
      else
         fBinSlider->SetPosition(1);
      fCancel->SetState(kButtonDisabled);
      fApply->SetState(kButtonDisabled);
      Update();
      delete [] div;
   } else if (ret==2) DoCancel();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the Cancel button of the Binning tab.

void TH1Editor::DoCancel()
{
   if (fBinHist) {
      fGedEditor->GetPad()->cd();
      fHist->Reset();
      fHist->SetBins(fBinHist->GetXaxis()->GetNbins(),
                     fBinHist->GetXaxis()->GetXmin(),
                     fBinHist->GetXaxis()->GetXmax());
      fHist->Add(fBinHist);
      fHist->GetXaxis()->SetRange(fBinHist->GetXaxis()->GetFirst(),
                                  fBinHist->GetXaxis()->GetLast());
      delete fBinHist;
      fBinHist = 0;
      fCancel->SetState(kButtonDisabled);
      fApply->SetState(kButtonDisabled);
      Int_t* divx = Dividers(fHist->GetXaxis()->GetNbins());
      if (divx[0]!=2) fBinSlider->SetPosition(1);
      // Consigning the new Histogram to all other Editors
//      fGedEditor->GetPad()->GetCanvas()->Selected(fGedEditor->GetPad(), fHist,  0);
      Update();
      //fModel = fHist;
      Refresh(fHist);
      delete [] divx;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the selected histogram type (HIST, LEGO1-2, SURF1-5).

TString TH1Editor::GetHistTypeLabel()
{
   TString s="";
   switch (fTypeCombo->GetSelected()){
      case (-1)         : {s = "LEGO"; break;}
      case (kTYPE_LEGO ): {s = "LEGO"; break;}
      case (kTYPE_LEGO1): {s = "LEGO1"; break;}
      case (kTYPE_LEGO2): {s = "LEGO2"; break;}
      case (kTYPE_SURF ): {s = "SURF"; break;}
      case (kTYPE_SURF1): {s = "SURF1"; break;}
      case (kTYPE_SURF2): {s = "SURF2"; break;}
      case (kTYPE_SURF3): {s = "SURF3"; break;}
      case (kTYPE_SURF4): {s = "SURF4"; break;}
      case (kTYPE_SURF5): {s = "SURF5"; break;}
      default:  break;
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the selected coordinate system of the histogram (POL,CYL,SPH,PSR).

TString TH1Editor::GetHistCoordsLabel()
{
   TString s="";
   if (fDim->GetState()!=kButtonDown) {
      switch (fCoordsCombo->GetSelected()){
         case (-1)         : {s = "POL"; break;}
         case (kCOORDS_CAR): {s = ""; break;}
         case (kCOORDS_POL): {s = "POL"; break;}
         case (kCOORDS_CYL): {s = "CYL"; break;}
         case (kCOORDS_SPH): {s = "SPH"; break;}
         case (kCOORDS_PSR): {s = "PSR"; break;}
         default:  break;
      }
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the selected error type (E,E1-5).

TString TH1Editor::GetHistErrorLabel()
{
   TString s="";
   switch (fErrorCombo->GetSelected()){
      case (-1)             : {s = ""; break;}
      case (kERRORS_NO)     : {s = ""; break;}
      case (kERRORS_SIMPLE) : {s = "E"; break;}
      case (kERRORS_EDGES)  : {s = "E1"; break;}
      case (kERRORS_REC)    : {s = "E2"; break;}
      case (kERRORS_FILL)   : {s = "E3"; break;}
      case (kERRORS_CONTOUR): {s = "E4"; break;}
      default:  break;
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the selected shape of the histogram (C, L, LF2).

TString TH1Editor::GetHistAddLabel()
{
   TString s="";
   switch (fAddCombo->GetSelected()){
      case (-1)         : {s = "" ; break;}
      case (kADD_NONE)  : {s = "" ; break;}
      case (kADD_SMOOTH): {s = "C"; break;}
      case (kADD_SIMPLE): {s = "L"; break;}
      case (kADD_FILL)  : {s = "LF2"; break;}
      default           :  break;
   }
   if (fAddMarker->GetState()==kButtonDown) s += "P";
   if (fAddB->GetState()==kButtonDown) s += "B";
   if (fAddBar->GetState()==kButtonDown){
      if (fMakeHBar->GetState()==kButtonDown) s+="H";
      switch (fPercentCombo->GetSelected()){
         case (kPER_0) : { s += "BAR" ; break;}
         case (kPER_10): { s += "BAR1"; break;}
         case (kPER_20): { s += "BAR2"; break;}
         case (kPER_30): { s += "BAR3"; break;}
         case (kPER_40): { s += "BAR4"; break;}
      }
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Create histogram type combo box.

TGComboBox* TH1Editor::BuildHistTypeComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("Lego" , kTYPE_LEGO);
   c->AddEntry("Lego1", kTYPE_LEGO1);
   c->AddEntry("Lego2", kTYPE_LEGO2);
   c->AddEntry("Surf" , kTYPE_SURF);
   c->AddEntry("Surf1", kTYPE_SURF1);
   c->AddEntry("Surf2", kTYPE_SURF2);
   c->AddEntry("Surf3", kTYPE_SURF3);
   c->AddEntry("Surf4", kTYPE_SURF4);
   c->AddEntry("Surf5", kTYPE_SURF5);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create coordinate system type combo box.

TGComboBox* TH1Editor::BuildHistCoordsComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("Cartesian", kCOORDS_CAR);
   c->AddEntry("Cylindric", kCOORDS_CYL);
   c->AddEntry("Polar", kCOORDS_POL);
   c->AddEntry("Rapidity", kCOORDS_PSR);
   c->AddEntry("Spheric", kCOORDS_SPH);
   TGListBox* lb = c->GetListBox();
   lb->Resize(lb->GetWidth(), 83);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create error type combo box.

TGComboBox* TH1Editor::BuildHistErrorComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("No Errors", kERRORS_NO);
   c->AddEntry("Simple", kERRORS_SIMPLE);
   c->AddEntry("Edges", kERRORS_EDGES);
   c->AddEntry("Rectangles",kERRORS_REC);
   c->AddEntry("Fill", kERRORS_FILL);
   c->AddEntry("Contour", kERRORS_CONTOUR);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create Line/Bar combo box.

TGComboBox* TH1Editor::BuildHistAddComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("No Line", kADD_NONE);
   c->AddEntry("Simple Line", kADD_SIMPLE);
   c->AddEntry("Smooth Line", kADD_SMOOTH);
   c->AddEntry("Fill Area",kADD_FILL);
   TGListBox* lb = c->GetListBox();
   lb->Resize(lb->GetWidth(), 76);
   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create Percentage combo box for bar option.

TGComboBox* TH1Editor::BuildPercentComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry(" 0 %", kPER_0);
   c->AddEntry("10 %", kPER_10);
   c->AddEntry("20 %", kPER_20);
   c->AddEntry("30 %", kPER_30);
   c->AddEntry("40 %", kPER_40);
   TGListBox* lb = c->GetListBox();
   lb->Resize(lb->GetWidth(), 83);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Change the error combo box entry.

void TH1Editor::ChangeErrorCombo(Int_t i)
{
   switch (i){
      case 0: {
         if (((TGLBContainer*)((TGListBox*)fErrorCombo->GetListBox())->GetContainer())->GetPos(kERRORS_EDGES)!=-1)
            fErrorCombo->RemoveEntries(kERRORS_EDGES,kERRORS_CONTOUR);
         if (!((fErrorCombo->GetSelected()== kERRORS_NO) || (fErrorCombo->GetSelected()== kERRORS_SIMPLE)))
            fErrorCombo->Select(kERRORS_NO);
         TGListBox* lb = fErrorCombo->GetListBox();
         lb->Resize(lb->GetWidth(),36);
         break;
      }
      case 1: {
         if (((TGLBContainer*)((TGListBox*)fErrorCombo->GetListBox())->GetContainer())->GetPos(kERRORS_EDGES)==-1) {
            fErrorCombo->AddEntry("Edges", kERRORS_EDGES);
            fErrorCombo->AddEntry("Rectangles",kERRORS_REC);
            fErrorCombo->AddEntry("Fill", kERRORS_FILL);
            fErrorCombo->AddEntry("Contour", kERRORS_CONTOUR);
            TGListBox* lb = fErrorCombo->GetListBox();
            lb->Resize(lb->GetWidth(),100);
         }
         break;
      }
   }
}

 ///////////////////////////////////////////////////////////////////////////////
 /// Paint a 3D box.

void TH1Editor::PaintBox3D(Float_t *p1, Float_t *p2,Float_t *p3, Float_t *p4)
{
   if (fGedEditor->GetPad()->GetCanvas())
      fGedEditor->GetPad()->GetCanvas()->FeedbackMode(kTRUE);
   fGedEditor->GetPad()->SetLineWidth(1);
   fGedEditor->GetPad()->SetLineColor(2);
   fGedEditor->GetPad()->cd();
   fGedEditor->GetPad()->PaintLine3D(p1, p2);
   fGedEditor->GetPad()->PaintLine3D(p2, p3);
   fGedEditor->GetPad()->PaintLine3D(p3, p4);
   fGedEditor->GetPad()->PaintLine3D(p4, p1);
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of dividers of n (without the trivial divider n).
/// The number of dividers is saved in the first entry.

Int_t* TH1Editor::Dividers(Int_t n)
{
   Int_t* div;
   if (n <= 0) {
      div = new Int_t[1];
      div[0]=0;
   } else if (n == 1) {
      div = new Int_t[2];
      div[0]=div[1]=1;
   } else {
      div = new Int_t[(Int_t) n/2+2];
      div[0]=0;
      div[1]=1;

      Int_t num = 1;
      for (Int_t i=2; i <= n/2; i++) {
         if (n % i == 0) {
            num++;
            div[num] = i;
         }
      }
      num++;
      div[num]=n;
      div[0] = num;
//   for (Int_t a=0; a <= div[0]; a++) printf("div[%d] = %d\n", a , div[a]);
   }
   return div;
}

////////////////////////////////////////////////////////////////////////////////
/// If the contained histogram obj is deleted we must set its pointer to zero

void TH1Editor::RecursiveRemove(TObject* obj)
{
   if (obj == fHist) {
      fHist = 0;
   }
}
