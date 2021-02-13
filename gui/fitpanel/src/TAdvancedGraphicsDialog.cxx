// @(#)root/fitpanel:$Id$
// Author: David Gonzalez Maline 11/12/2008

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAdvancedGraphicsDialog                                              //
//                                                                      //
// Allows to create advanced graphics from the last fit made in the     //
// fitpanel. This includes the scan graphics, the contour and the       //
// confidence levels.                                                   //
//////////////////////////////////////////////////////////////////////////

#include "TAdvancedGraphicsDialog.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TVirtualPad.h"
#include "TColor.h"

#include "Fit/BinData.h"
#include "Math/IParamFunction.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"

#include <sstream>

#include <string>

TAdvancedGraphicsDialog::TAdvancedGraphicsDialog(const TGWindow *p, const TGWindow *main):
   TGTransientFrame(p, main, 10, 10, kVerticalFrame),
   fFitter((TBackCompFitter *) TVirtualFitter::GetFitter())
{
   // Creates the Advanced Graphics Dialog.


   if (!p && !main) {
      MakeZombie();
      return;
   }
   SetCleanup(kDeepCleanup);

   fMainFrame = new TGVerticalFrame(this);

   fTab = new TGTab(fMainFrame, 10, 10);
   fMainFrame->AddFrame(fTab, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX, 5,5,5,0));
   fTab->SetCleanup(kDeepCleanup);
   fTab->Associate(this);

   // Add the first method to the dialog (Contour)
   CreateContourFrame();
   fTab->AddTab("Contour", fContourFrame);

   // Add the second method to the dialog (Scan)
   CreateScanFrame();
   fTab->AddTab("Scan", fScanFrame);

   CreateConfFrame();
   fTab->AddTab("Conf Intervals", fConfFrame);

   TGCompositeFrame * frame = new TGHorizontalFrame(fMainFrame);

   fDraw = new TGTextButton(frame, "&Draw", kAGD_BDRAW);
   fDraw->Associate(this);
   frame->AddFrame(fDraw, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0));

   fClose = new TGTextButton(frame, "&Close", kAGD_BCLOSE);
   fClose->Associate(this);
   frame->AddFrame(fClose, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0));

   UInt_t width = 0, height = 0;
   height = fClose->GetDefaultHeight();
   width  = TMath::Max(width, fClose->GetDefaultWidth());
   frame->Resize((width + 20) * 2, height);

   fMainFrame->AddFrame(frame, new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 5, 0));

   this->AddFrame(fMainFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX,0,0,5,5));

   ConnectSlots();

   SetWindowName("Advanced Drawing Tools");

   // map all widgets and calculate size of dialog
   MapSubwindows();

   width  = GetDefaultWidth();
   height = GetDefaultHeight();

   Resize(width, height);
   MapWindow();

   // position relative to the parent's window
   CenterOnParent();

   // make the message box non-resizable
   SetWMSize(width, height);
   SetWMSizeHints(width, height, width, height, 0, 0);

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
               kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
               kMWMFuncMinimize,
               kMWMInputModeless);

   // popup dialog and wait till user replies
   gClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create the frame that contains all the necessary information for
/// the Contour method.

void TAdvancedGraphicsDialog::CreateContourFrame()
{
   fContourFrame = new TGVerticalFrame(fTab);
   TGHorizontalFrame* frame = new TGHorizontalFrame(fContourFrame);

   TGLabel* label = new TGLabel(frame, "Number of Points: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fContourPoints = new TGNumberEntry(frame, 40,
                                5, kAGD_SCANMIN,
                                TGNumberFormat::kNESInteger,
                                TGNumberFormat::kNEAPositive,
                                TGNumberFormat::kNELNoLimits);
   fContourPoints->Resize(130, 20);
   fContourPoints->GetNumberEntry()->SetToolTipText("Sets the number of points used for the contour");
   frame->AddFrame(fContourPoints, new TGLayoutHints(kLHintsNormal, 8, 0, 5, 0));
   fContourFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fContourFrame);
   label = new TGLabel(frame, "Parameter 1: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fContourPar1 = new TGComboBox(frame, kAGD_CONTPAR1);
   AddParameters(fContourPar1);
   fContourPar1->Resize(130, 20);
   fContourPar1->Associate(this);
   TGListBox *lb = fContourPar1->GetListBox();
   lb->Resize(lb->GetWidth(), 200);
   frame->AddFrame(fContourPar1, new TGLayoutHints(kLHintsNormal, 37, 0, 5, 0));
   fContourFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fContourFrame);

   label = new TGLabel(frame, "Parameter 2: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fContourPar2 = new TGComboBox(frame, kAGD_CONTPAR2);
   AddParameters(fContourPar2);
   fContourPar2->Select(kAGD_PARCOUNTER+1, kFALSE);
   fContourPar2->Resize(130, 20);
   fContourPar2->Associate(this);
   lb = fContourPar2->GetListBox();
   lb->Resize(lb->GetWidth(), 200);
   frame->AddFrame(fContourPar2, new TGLayoutHints(kLHintsNormal, 37, 0, 5, 0));

   fContourFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fContourFrame);

   label = new TGLabel(frame, "Confidence Level: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fContourError = new TGNumberEntry(frame, 0.683, 5, kAGD_CONTERR,
                                     TGNumberFormat::kNESRealThree,
                                     TGNumberFormat::kNEANonNegative,
                                     TGNumberFormat::kNELNoLimits);
   fContourError->Resize(130, 20);
   fContourError->GetNumberEntry()->SetToolTipText("Sets the contour confidence level");
   frame->AddFrame(fContourError, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));

   fContourFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 5));

   frame = new TGHorizontalFrame(fContourFrame);

   label = new TGLabel(frame, "Fill Colour: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fContourColor = new TGColorSelect(frame, TColor::Number2Pixel(kYellow - 10), kAGD_CONTCOLOR);
   frame->AddFrame(fContourColor, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));

   fContourOver = new TGCheckButton(frame, "Superimpose", kAGD_CONTOVER);
   fContourOver->SetToolTipText("If checked, the new contour will overlap the previous one");
   frame->AddFrame(fContourOver, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));

   fContourFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 5));
}

////////////////////////////////////////////////////////////////////////////////
/// Create the frame that contains all the necessary information for
/// the Scan method.

void TAdvancedGraphicsDialog::CreateScanFrame()
{
   fScanFrame = new TGVerticalFrame(fTab);
   TGHorizontalFrame* frame = new TGHorizontalFrame(fScanFrame);

   TGLabel* label = new TGLabel(frame, "Number of Points: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fScanPoints = new TGNumberEntry(frame, 40,
                                5, kAGD_SCANMIN,
                                TGNumberFormat::kNESInteger,
                                TGNumberFormat::kNEAPositive,
                                TGNumberFormat::kNELNoLimits);
   fScanPoints->Resize(140, 20);
   fScanPoints->GetNumberEntry()->SetToolTipText("Sets the number of points used in the scan");
   frame->AddFrame(fScanPoints, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
   fScanFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fScanFrame);

   label = new TGLabel(frame, "Parameter: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fScanPar = new TGComboBox(frame, kAGD_SCANPAR);
   AddParameters(fScanPar);
   fScanPar->Resize(140, 20);
   fScanPar->Associate(this);
   TGListBox *lb = fScanPar->GetListBox();
   lb->Resize(lb->GetWidth(), 200);
   frame->AddFrame(fScanPar, new TGLayoutHints(kLHintsNormal, 39, 0, 5, 0));
   fScanFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fScanFrame);
   label = new TGLabel(frame, "Min: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   double val = fFitter->GetParameter( fScanPar->GetSelected() - kAGD_PARCOUNTER );
   double err = fFitter->GetParError( fScanPar->GetSelected() - kAGD_PARCOUNTER );
   fScanMin = new TGNumberEntry(frame, val - 2.*err ,
                                5, kAGD_SCANMIN,
                                TGNumberFormat::kNESRealFour,
                                TGNumberFormat::kNEAAnyNumber,
                                TGNumberFormat::kNELNoLimits);
   fScanMin->Resize(70, 20);
   fScanMin->GetNumberEntry()->SetToolTipText("Sets the minimum parameter value");
   frame->AddFrame(fScanMin, new TGLayoutHints(kLHintsNormal, 2, 0, 5, 0));

   label = new TGLabel(frame, "Max: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 35, 5, 5, 0));
   fScanMax = new TGNumberEntry(frame, val + 2.*err,
                                5, kAGD_SCANMAX,
                                TGNumberFormat::kNESRealFour,
                                TGNumberFormat::kNEAAnyNumber,
                                TGNumberFormat::kNELNoLimits);
   fScanMax->Resize(70, 20);
   fScanMax->GetNumberEntry()->SetToolTipText("Sets the maximum parameter value");
   frame->AddFrame(fScanMax, new TGLayoutHints(kLHintsNormal, 2, 0, 5, 0));
   fScanFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

}

////////////////////////////////////////////////////////////////////////////////
/// Create the frame that contains all the necessary information for
/// the Confidence Level  method.

void TAdvancedGraphicsDialog::CreateConfFrame()
{
   fConfFrame = new TGVerticalFrame(fTab);
   TGHorizontalFrame* frame = new TGHorizontalFrame(fConfFrame);

   TGLabel* label = new TGLabel(frame, "Confidence Level: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fConfLevel = new TGNumberEntry(frame, 0.95,
                                  5, kAGD_SCANMIN,
                                  TGNumberFormat::kNESRealTwo,
                                  TGNumberFormat::kNEAPositive,
                                  TGNumberFormat::kNELLimitMinMax,
                                  0, 0.9999);
   fConfLevel->Resize(140, 20);
   fConfLevel->GetNumberEntry()->SetToolTipText("Sets the value of the confidence level");
   frame->AddFrame(fConfLevel, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));
   fConfFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   frame = new TGHorizontalFrame(fConfFrame);

   label = new TGLabel(frame, "Fill Colour: ");
   frame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   fConfColor = new TGColorSelect(frame, TColor::Number2Pixel(kYellow - 10), kAGD_CONTCOLOR);
   frame->AddFrame(fConfColor, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));

   fConfFrame->AddFrame(frame, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 5));
}

void TAdvancedGraphicsDialog::AddParameters(TGComboBox* comboBox)
{
   // Add all the parameters of the VirtualFitter into a comboBox
   // (helper method)

   for ( Int_t i = 0; i < fFitter->GetNumberTotalParameters(); ++i ) {
      comboBox->AddEntry(fFitter->GetParName(i), kAGD_PARCOUNTER + i);
   }
   comboBox->Select(kAGD_PARCOUNTER, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Connect the slots (buttons mainly + specific methods)

void TAdvancedGraphicsDialog::ConnectSlots()
{
   // Buttons
   fClose->Connect("Clicked()", "TAdvancedGraphicsDialog", this, "CloseWindow()");
   fDraw->Connect("Clicked()", "TAdvancedGraphicsDialog", this, "DoDraw()");

   // Slots for the Scan method
   fScanPar->Connect("Selected(Int_t)", "TAdvancedGraphicsDialog", this, "DoChangedScanPar(Int_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the Min and Max default values of the scan method,
/// depending on the selected parameter.

void TAdvancedGraphicsDialog::DoChangedScanPar(Int_t selected)
{
   double val = fFitter->GetParameter( selected - kAGD_PARCOUNTER );
   double err = fFitter->GetParError(  selected - kAGD_PARCOUNTER );
   fScanMin->SetNumber( val -2 * err );
   fScanMax->SetNumber( val +2 * err );
}

////////////////////////////////////////////////////////////////////////////////
/// Calls the correspoding method, depending on the selected tab.

void TAdvancedGraphicsDialog::DoDraw()
{
   if ( fTab->GetCurrent() == 0 ) {
      DrawContour();
   } else if ( fTab->GetCurrent() == 1 ) {
      DrawScan();
   } else if ( fTab->GetCurrent() == 2 ) {
      DrawConfidenceLevels();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generates all necessary data for the Contour method from its
/// tab. Then it call Virtual Fitter to perform it.

void TAdvancedGraphicsDialog::DrawContour()
{
   static TGraph * graph = 0;
   std::string options;
   if ( ! (fContourOver->GetState() == kButtonDown) ) {
      if ( graph )
         delete graph;
      options = "ALF";
   } else
      options = "LF";
   graph = new TGraph( static_cast<int>(fContourPoints->GetNumber()) );
   Int_t par1 = fContourPar1->GetSelected() - kAGD_PARCOUNTER;
   Int_t par2 = fContourPar2->GetSelected() - kAGD_PARCOUNTER;
   if ( par1 == par2 ) {
      Error("TAdvancedGraphicsDialog::DrawContour", "Parameters cannot be the same");
      return;
   }
   // contour error is actually the desired confidence level
   Double_t cl = fContourError->GetNumber();
   fFitter->Contour( par1, par2, graph, cl);
   graph->SetFillColor( TColor::GetColor( fContourColor->GetColor() ) );
   graph->GetXaxis()->SetTitle( fFitter->GetParName(par1) );
   graph->GetYaxis()->SetTitle( fFitter->GetParName(par2) );
   graph->Draw( options.c_str() );
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Generates all necessary data for the Scan method from its
/// tab. Then it call Virtual Fitter to perform it.

void TAdvancedGraphicsDialog::DrawScan()
{
   static TGraph * graph = 0;
   if ( graph )
      delete graph;
   graph = new TGraph( static_cast<int>(fScanPoints->GetNumber()) );
   Int_t par = fScanPar->GetSelected() - kAGD_PARCOUNTER;
   fFitter->Scan( par, graph,
                  fScanMin->GetNumber(),
                  fScanMax->GetNumber() );
   graph->SetLineColor(kBlue);
   graph->SetLineWidth(2);
   graph->GetXaxis()->SetTitle(fFitter->GetParName(par) );
   graph->GetYaxis()->SetTitle("FCN" );
   graph->Draw("APL");
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Generates all necessary data for the Scan method from its
/// tab. Then it call Virtual Fitter to perform it.

void TAdvancedGraphicsDialog::DrawConfidenceLevels()
{
   const ROOT::Fit::FitResult& result = fFitter->GetFitResult();
   const ROOT::Fit::FitResult::IModelFunction* function = result.FittedFunction();
   const ROOT::Fit::BinData* data = dynamic_cast<const ROOT::Fit::BinData*>(&(fFitter->GetFitData()));
   if ( !data )
   {
      Error("DrawConfidenceLevels","Unbinned data set cannot draw confidence levels.");
      return;
   }

   if ( !function )
   {
      Error("DrawConfidenceLevels","Fit Function does not exist!");
      return;
   }

   std::vector<Double_t> ci(data->Size());
   result.GetConfidenceIntervals(*data, &ci[0], fConfLevel->GetNumber());

   if ( data->NDim() == 1 )
   {
      TGraphErrors* g = new TGraphErrors(ci.size());
      for (unsigned int i = 0; i < ci.size(); ++i)
      {
         const Double_t *x = data->Coords(i);
         const Double_t y = (*function)(x);
         g->SetPoint(i, *x, y);
         g->SetPointError(i, 0, ci[i]);
      }
      std::ostringstream os;
      os << "Confidence Intervals with " << fConfLevel->GetNumber()
         << " conf. band.";
      g->SetTitle(os.str().c_str());
      g->SetLineColor( TColor::GetColor( fConfColor->GetColor() ));
      g->SetFillColor( TColor::GetColor( fConfColor->GetColor() ));
      g->SetFillStyle(3001);
      g->Draw("C3same");
   } else if ( data->NDim() == 2 )
   {
      TGraph2DErrors* g = new TGraph2DErrors(ci.size());
      for (unsigned int i = 0; i < ci.size(); ++i)
      {
         const Double_t *x = data->Coords(i);
         const Double_t y = (*function)(x);
         g->SetPoint(i, x[0], x[1], y);
         g->SetPointError(i, 0, 0, ci[i]);
      }
      std::ostringstream os;
      os << "Confidence Intervals with " << fConfLevel->GetNumber()
         << " conf. band.";
      g->SetTitle(os.str().c_str());
      g->SetLineColor( TColor::GetColor( fConfColor->GetColor() ));
      g->SetFillColor( TColor::GetColor( fConfColor->GetColor() ));
      g->SetFillStyle(3001);
      g->Draw("C3same");
   }
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup dialog.

TAdvancedGraphicsDialog::~TAdvancedGraphicsDialog()
{
   Cleanup();
}
