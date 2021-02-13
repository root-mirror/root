// @(#)root/fitpanel:$Id: ed8d59036b6a51c67cd739c2c75aa7780b847bf8 $
// Author: Ilka Antcheva, Lorenzo Moneta 10/08/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFitEditor                                                           //
//                                                                      //
// Allows to perform, explore and compare various fits.                 //
//                                                                      //
// To display the new Fit panel interface right click on a histogram    //
// or a graph to pop up the context menu and then select the menu       //
// entry 'Fit Panel'.                                                   //
//                                                                      //
// "General" Tab                                                        //
//                                                                      //
// The first set of GUI elements is related to the function choice      //
// and settings. The status bar on the bottom provides information      //
// about the current minimization settings using the following          //
// abbreviations:                                                       //
// LIB - shows the current choice between Minuit/Minuit2/Fumili         //
// MIGRAD or FUMILI points to the current minimization method in use.   //
// Itr: - shows the maximum number of iterations nnnn set for the fit.  //
// Prn: - can be DEF/VER/QT and shows the current print option in use.  //
//                                                                      //
// "Predefined" combo box - contains a list of predefined functions     //
// in ROOT. The default one is Gaussian.                                //
//                                                                      //
// "Operation" radio button group defines selected operational mode     //
// between functions: NOP - no operation (default); ADD - addition      //
// CONV - convolution (will be implemented in the future).              //
//                                                                      //
// Users can enter the function expression in a text entry field.       //
// The entered string is checked after Enter key was pressed. An        //
// error message shows up if the string is not accepted. The current    //
// prototype is limited and users have no freedom to enter file/user    //
// function names in this field.                                        //
//                                                                      //
// "Set Parameters" button opens a dialog for parameters settings.      //
//                                                                      //
// "Fit Settings" provides user interface elements related to the       //
// fitter. Currently there are two method choices: Chi-square and       //
// Binned Likelihood.                                                   //
//                                                                      //
// "Linear Fit" check button sets the use of Linear fitter is it is     //
// selected. Otherwise the option 'F' is applied if polN is selected.   //
// "Robust" number entry sets the robust value when fitting graphs.     //
// "No Chi-square" check button sets ON/OFF option 'C' - do not         //
// calculate Chi-square (for Linear fitter).                            //
//                                                                      //
// Fit options:                                                         //
// "Integral" check button switch ON/OFF option 'I' - use integral      //
// of function instead of value in bin center.                          //
// "Best Errors" sets ON/OFF option 'E' - better errors estimation      //
// using Minos technique.                                               //
// "All weights = 1" sets ON/OFF option 'W' - all weights set to 1,     //
// excluding empty bins and ignoring error bars.                        //
// "Empty bins, weights=1" sets ON/OFF option 'WW' -  all weights       //
// equal to 1, including  empty bins, error bars ignored.               //
// "Use range" sets ON/OFF option 'R' - fit only data within the        //
// specified function range with the slider.                            //
// "Improve fit results" sets ON/OFF option 'M' - after minimum is      //
// found, search for a new one.                                         //
// "Add to list" sets On/Off option '+'- add function to the list       //
// without deleting the previous.                                       //
//                                                                      //
// Draw options:                                                        //
// "SAME" sets On/Off function drawing on the same pad.                 //
// "No drawing" sets On/Off option '0'- do not draw function graphics.  //
// "Do not store/draw" sets On/Off option 'N'- do not store the         //
// function, do not draw it.                                            //
//                                                                      //
// Sliders settings are used if option 'R' - use range is active.       //
// Users can change min/max values by pressing the left mouse button    //
// near to the left/right slider edges. It is possible o change both    //
// values simultaneously by pressing the left mouse button near to its  //
// center and moving it to a new desire position.                       //
//                                                                      //
// "Minimization" Tab                                                   //
//                                                                      //
// "Library" group allows you to use Minuit, Minuit2 or Fumili          //
// minimization packages for your fit.                                  //
//  "Minuit" - the popular Minuit minimization package.                 //
//  "Minuit2" - a new object-oriented implementation of Minuit in C++.  //
//  "Fumili" - the popular Fumili minimization package.                 //
//                                                                      //
// "Method" group has currently restricted functionality.               //
//  "MIGRAD" method is available for Minuit and Minuit2                 //
//  "FUMILI" method is available for Fumili and Minuit2                 //
//  "SIMPLEX" method is disabled (will come with the new fitter design) //
//                                                                      //
// "Minimization Settings' group allows users to set values for:        //
//  "Error definition" - between 0.0 and 100.0  (default is 1.0).       //
//  "Maximum tolerance" - the fit relative precision in use.            //
//  "Maximum number of iterations" - default is 5000.                   //
//                                                                      //
// Print options:                                                       //
//  "Default" - between Verbose and Quiet.                              //
//  "Verbose" - prints results after each iteration.                    //
//  "Quiet" - no fit information is printed.                            //
//                                                                      //
// Fit button - performs a fit.                                         //
// Reset - resets all GUI elements and related fit settings to the      //
// default ones.                                                        //
// Close - closes this window.                                          //
//                                                                      //
// Begin_Html                                                           //
/*
<img src="gif/TFitEditor.gif">
*/
//End_Html
//////////////////////////////////////////////////////////////////////////

#include "TFitEditor.h"
#include "TROOT.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TGTab.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGComboBox.h"
#include "TGTextEntry.h"
#include "TGGC.h"
#include "TGButtonGroup.h"
#include "TGNumberEntry.h"
#include "TGDoubleSlider.h"
#include "TGStatusBar.h"
#include "TFitParametersDialog.h"
#include "TGMsgBox.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "HFitInterface.h"
#include "TF1.h"
#include "TF1NormSum.h"
#include "TF1Convolution.h"
#include "TF2.h"
#include "TF3.h"
#include "THStack.h"
#include "Fit/UnBinData.h"
#include "Fit/BinData.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TVirtualTreePlayer.h"
#include "TSelectorDraw.h"
#include "TTreeInput.h"
#include "TAdvancedGraphicsDialog.h"
#include "TVirtualX.h"
#include "strlcpy.h"

#include "RConfigure.h"
#include "TPluginManager.h"

#include <vector>
#include <queue>
using std::vector;
using std::pair;

#include "CommonDefs.h"

// #include <iostream>
// using std::cout;
// using std::endl;

void SearchCanvases(TSeqCollection* canvases, std::vector<TObject*>& objects);

typedef std::multimap<TObject*, TF1*> FitFuncMap_t;

////////////////////////////////////////////////////////////////////////////////
/// This method looks among the functions stored by the fitpanel, the
/// one that is currently selected in the fFuncList

TF1* TFitEditor::FindFunction()
{
   // Get the title/name of the function from fFuncList
   TGTextLBEntry *te = (TGTextLBEntry *)fFuncList->GetSelectedEntry();
   if ( !te ) return 0;
   TString name(te->GetTitle());

   // Look for a system function if it's USER DEFINED function
   if ( fTypeFit->GetSelected() == kFP_UFUNC ) {
      for (auto f : fSystemFuncs) {
         if ( strcmp( f->GetName(), name ) == 0 )
            // If found, return it.
            return f;
      }
   // If we are looking for previously fitted functions, look in the
   // fPrevFit data structure.
   } else if ( fTypeFit->GetSelected() == kFP_PREVFIT ) {
      std::pair<fPrevFitIter, fPrevFitIter> look = fPrevFit.equal_range(fFitObject);
      for ( fPrevFitIter it = look.first; it != look.second; ++it ) {
         TF1* f = it->second;
         if ( strcmp( f->GetName(), name ) == 0 )
            // If found, return it
            return f;
      }
   }

   // Return a pointer to null if the function does not exist. This
   // will eventually create a segmentation fault, but the line should
   // never be executed.
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///Copies f into a new TF1 to be stored in the fitpanel with it's
///own ownership. This is taken from Fit::StoreAndDrawFitFunction in
///HFitImpl.cxx

TF1* copyTF1(TF1 *f)
{
   double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;

   // no need to use kNotGlobal bit. TF1::Copy does not add in the list by default
   if ( dynamic_cast<TF3 *>(f) != 0 ) {
      TF3* fnew = (TF3 *)f->IsA()->New();
      f->Copy(*fnew);
      f->GetRange(xmin,ymin,zmin,xmax,ymax,zmax);
      fnew->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
      fnew->SetParent( nullptr );
      fnew->AddToGlobalList(false);
      return fnew;
   } else if ( dynamic_cast<TF2 *>(f) != 0 ) {
      TF2* fnew = (TF2 *)f->IsA()->New();
      f->Copy(*fnew);
      f->GetRange(xmin,ymin,xmax,ymax);
      fnew->SetRange(xmin,ymin,xmax,ymax);
      fnew->Save(xmin,xmax,ymin,ymax,0,0);
      fnew->SetParent( nullptr );
      fnew->AddToGlobalList(false);
      return fnew;
   } else {
      TF1* fnew = (TF1 *)f->IsA()->New();
      f->Copy(*fnew);
      f->GetRange(xmin,xmax);
      fnew->SetRange(xmin,xmax);
      // This next line is added, as fnew-Save fails with gausND! As
      // the number of dimensions is unknown...
      if ( '\0' != fnew->GetExpFormula()[0] )
         fnew->Save(xmin,xmax,0,0,0,0);
      fnew->SetParent( nullptr );
      fnew->AddToGlobalList(false);
      return fnew;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stores the parameters of the given function into pars

void GetParameters(TFitEditor::FuncParams_t & pars, TF1* func)
{
   int npar = func->GetNpar();
   if (npar != (int) pars.size() ) pars.resize(npar);
   for ( Int_t i = 0; i < npar; ++i )
   {
      Double_t par_min, par_max;
      pars[i][PAR_VAL] = func->GetParameter(i);
      func->GetParLimits(i, par_min, par_max);
      pars[i][PAR_MIN] = par_min;
      pars[i][PAR_MAX] = par_max;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Restore the parameters from pars into the function

void SetParameters(TFitEditor::FuncParams_t & pars, TF1* func)
{
   int npar = func->GetNpar();
   if (npar > (int) pars.size() ) pars.resize(npar);
   for ( Int_t i = 0; i < npar; ++i )
   {
      func->SetParameter(i, pars[i][PAR_VAL]);
      func->SetParLimits(i, pars[i][PAR_MIN], pars[i][PAR_MAX]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Parameter initialization for the function

template<class FitObject>
void InitParameters(TF1* func, FitObject * fitobj)
{
   const int special = func->GetNumber();
   if (100 == special || 400 == special) {
      ROOT::Fit::BinData data;
      ROOT::Fit::FillData(data,fitobj,func);
      ROOT::Fit::InitGaus(data, func);
      // case gaussian or Landau
   } else if ( 110 == special || 410 == special ) {
      ROOT::Fit::BinData data;
      ROOT::Fit::FillData(data,fitobj,func);
      ROOT::Fit::Init2DGaus(data,func);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Splits the entry in fDataSet to get the selected variables and cuts
/// from the text.

void GetTreeVarsAndCuts(TGComboBox* dataSet, TString& variablesStr, TString& cutsStr)
{
   // Get the entry
   TGTextLBEntry* textEntry =
      static_cast<TGTextLBEntry*>( dataSet->GetListBox()->GetEntry( dataSet->GetSelected() ) );
   if (!textEntry) return;
   // Get the name of the tree
   TString nameStr ( textEntry->GetText()->GetString() );
   // Get the variables selected
   variablesStr = nameStr(nameStr.First('(') + 2, nameStr.First(',') - nameStr.First('(') - 3);
   // Get the cuts selected
   cutsStr = nameStr( nameStr.First(',') + 3, nameStr.First(')') - nameStr.First(',') - 4 );
}


ClassImp(TFitEditor);

TFitEditor *TFitEditor::fgFitDialog = 0;

////////////////////////////////////////////////////////////////////////////////
/// Static method - opens the fit panel.

TFitEditor * TFitEditor::GetInstance(TVirtualPad* pad, TObject *obj)
{
   // Get the default pad if not provided.
   if (!pad)
   {
      if (!gPad)
         gROOT->MakeDefCanvas();
      pad = gPad;
   }

   if (!fgFitDialog) {
      fgFitDialog = new TFitEditor(pad, obj);
   } else {
      fgFitDialog->Show(pad, obj);
   }
   return fgFitDialog;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor of fit editor. 'obj' is the object to be fitted and
/// 'pad' where it is drawn.

TFitEditor::TFitEditor(TVirtualPad* pad, TObject *obj) :
   TGMainFrame(gClient->GetRoot(), 20, 20),
   fParentPad   (0),
   fFitObject   (0),
   fDim         (0),
   fXaxis       (0),
   fYaxis       (0),
   fZaxis       (0),
   fSumFunc     (0),
   fConvFunc    (0),
   fFuncPars    (0),
   fChangedParams (kFALSE)
{
   fType = kObjectHisto;
   SetCleanup(kDeepCleanup);

   TGCompositeFrame *tf = new TGCompositeFrame(this, 350, 26,
                                                kHorizontalFrame);
   TGLabel *label = new TGLabel(tf,"Data Set: ");
   tf->AddFrame(label, new TGLayoutHints(kLHintsNormal, 15, 0, 5, 0));

   fDataSet = new TGComboBox(tf, kFP_DATAS);
   FillDataSetList();
   fDataSet->Resize(264, 20);

   tf->AddFrame(fDataSet, new TGLayoutHints(kLHintsNormal, 13, 0, 5, 0));
   fDataSet->Associate(this);

   this->AddFrame(tf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX,0,0,5,5));

   CreateFunctionGroup();

   fTab = new TGTab(this, 10, 10);
   AddFrame(fTab, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));
   fTab->SetCleanup(kDeepCleanup);
   fTab->Associate(this);

   TGHorizontalFrame *cf1 = new TGHorizontalFrame(this, 350, 20, kFixedWidth);
   cf1->SetCleanup(kDeepCleanup);
   fUpdateButton = new TGTextButton(cf1, "&Update", kFP_UPDATE);
   fUpdateButton->Associate(this);
   cf1->AddFrame(fUpdateButton, new TGLayoutHints(kLHintsTop |
                                                  kLHintsExpandX, 0, 20, 2, 2));


   fFitButton = new TGTextButton(cf1, "&Fit", kFP_FIT);
   fFitButton->Associate(this);
   cf1->AddFrame(fFitButton, new TGLayoutHints(kLHintsTop |
                                               kLHintsExpandX, 15, -6, 2, 2));

   fResetButton = new TGTextButton(cf1, "&Reset", kFP_RESET);
   fResetButton->Associate(this);
   cf1->AddFrame(fResetButton, new TGLayoutHints(kLHintsTop |
                                                 kLHintsExpandX, 11, -2, 2, 2));

   fCloseButton = new TGTextButton(cf1, "&Close", kFP_CLOSE);
   fCloseButton->Associate(this);
   cf1->AddFrame(fCloseButton, new TGLayoutHints(kLHintsTop |
                                                 kLHintsExpandX, 7, 2, 2, 2));
   AddFrame(cf1, new TGLayoutHints(kLHintsNormal |
                                   kLHintsRight, 0, 5, 5, 5));

   // Create status bar
   int parts[] = { 20, 20, 20, 20, 20 };
   fStatusBar = new TGStatusBar(this, 10, 10);
   fStatusBar->SetParts(parts, 5);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom |
                                          kLHintsLeft   |
                                          kLHintsExpandX));

   CreateGeneralTab();
   CreateMinimizationTab();

   gROOT->GetListOfCleanups()->Add(this);

   MapSubwindows();
   fGeneral->HideFrame(fSliderZParent);

   // do not allow resizing
   TGDimension size = GetDefaultSize();
   SetWindowName("Fit Panel");
   SetIconName("Fit Panel");
   SetClassHints("ROOT", "Fit Panel");

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);

   ConnectSlots();

   GetFunctionsFromSystem();

   if (!obj) {
      TList* l = new TList();
      l->Add(pad);
      std::vector<TObject*> v;
      SearchCanvases(l, v);
      if ( v.size() )
         obj = v[0];
      delete l;
   }

   SetFitObject(pad, obj, kButton1Down);

   // In case we want to make it without a default canvas. This will
   // be implemented after the 5.21/06 Release. Remember to take out
   // any reference to the pad/canvas when the fitpanel is shown
   // and/or built.

   //SetCanvas(0 /*pad->GetCanvas()*/);

   if ( pad ) {
      SetCanvas(pad->GetCanvas());
      if ( obj )
         pad->GetCanvas()->Selected(pad, obj, kButton1Down);
   }

   UInt_t dw = fClient->GetDisplayWidth();
   UInt_t cw = 0;
   UInt_t cx = 0;
   UInt_t cy = 0;
   if (pad && pad->GetCanvas() ) {
      cw = pad->GetCanvas()->GetWindowWidth();
      cx = (UInt_t)pad->GetCanvas()->GetWindowTopX();
      cy = (UInt_t)pad->GetCanvas()->GetWindowTopY();
   }

   Resize(size);
   MapWindow();

   if (cw + size.fWidth < dw) {
      Int_t gedx = 0, gedy = 0;
      gedx = cx+cw+4;
      gedy = (cy > 20) ? cy-20 : 0;
      MoveResize(gedx, gedy, size.fWidth, size.fHeight);
      SetWMPosition(gedx, gedy);
   }

   gVirtualX->RaiseWindow(GetId());

   ChangeOptions(GetOptions() | kFixedSize);
   SetWMSize(size.fWidth, size.fHeight);
   SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Fit editor destructor.

TFitEditor::~TFitEditor()
{
   DisconnectSlots();

   // Disconnect all the slot that were no disconnected in DisconnecSlots
   fCloseButton ->Disconnect("Clicked()");
   fDataSet     ->Disconnect("Selected(Int_t)");
   fUpdateButton->Disconnect("Clicked()");
   TQObject::Disconnect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)",
                         this, "SetFitObject(TVirtualPad *, TObject *, Int_t)");
   gROOT->GetListOfCleanups()->Remove(this);

   //Clean up the members that are not automatically cleaned.
   Cleanup();
   delete fLayoutNone;
   delete fLayoutAdd;
   delete fLayoutConv;

   if (fConvFunc) delete fConvFunc;
   if (fSumFunc) delete fSumFunc;

   // release memory used by stored functions of previous fits
   for (auto &entry : fPrevFit)
      delete entry.second;
   fPrevFit.clear();

   // release memory used by copies of system functions
   for (auto func : fSystemFuncs)
      delete func;
   fSystemFuncs.clear();

   // Set the singleton reference to null
   fgFitDialog = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the Frame that contains oll the information about the
/// function.

void TFitEditor::CreateFunctionGroup()
{
   TGGroupFrame     *gf1 = new TGGroupFrame(this, "Fit Function", kFitWidth);
   TGCompositeFrame *tf0 = new TGCompositeFrame(gf1, 350, 26, kHorizontalFrame);
   TGLabel *label1       = new TGLabel(tf0,"Type:");
   tf0 -> AddFrame(label1, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 0));

   fTypeFit =  new TGComboBox(tf0, kFP_TLIST);
   fTypeFit -> AddEntry("User Func", kFP_UFUNC);
   fTypeFit -> AddEntry("Predef-1D", kFP_PRED1D);
   fTypeFit -> Resize(90, 20);
   fTypeFit -> Select(kFP_PRED1D, kFALSE);

   TGListBox *lb = fTypeFit->GetListBox();
   lb->Resize(lb->GetWidth(), 200);
   tf0->AddFrame(fTypeFit, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));
   fTypeFit->Associate(this);

   fFuncList = new TGComboBox(tf0, kFP_FLIST);
   FillFunctionList();
   fFuncList->Resize(194, 20);
   fFuncList->Select(kFP_GAUS, kFALSE);

   lb = fFuncList->GetListBox();
   lb -> Resize(lb->GetWidth(), 500);
   tf0 -> AddFrame(fFuncList, new TGLayoutHints(kLHintsNormal, 5, 0, 5, 0));
   fFuncList->Associate(this);

   gf1->AddFrame(tf0, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));

   TGCompositeFrame *tf1 = new TGCompositeFrame(gf1, 350, 26,  kHorizontalFrame);
   TGHButtonGroup *bgr   = new TGHButtonGroup(tf1, "Operation");

   bgr      -> SetRadioButtonExclusive();
   fNone    = new TGRadioButton(bgr, "Nop", kFP_NONE);
   fAdd     = new TGRadioButton(bgr, "Add", kFP_ADD);
   fNormAdd = new TGRadioButton(bgr, "NormAdd", kFP_NORMADD);
   fConv    = new TGRadioButton(bgr, "Conv", kFP_CONV);

   fNone    -> SetToolTipText("No operation defined");
   fNone    -> SetState(kButtonDown, kFALSE);
   fAdd     -> SetToolTipText("Addition");
  // fAdd     -> SetState(kButtonDown, kFALSE);
   fNormAdd -> SetToolTipText("NormAddition");
   //fNormAdd -> SetState(kButtonDown, kFALSE);
   fConv    -> SetToolTipText("Convolution");
   //fConv    -> SetState(kButtonDown, kTRUE);

   fLayoutNone    = new TGLayoutHints(kLHintsLeft,0 ,5,3,-10);
   fLayoutAdd     = new TGLayoutHints(kLHintsLeft,10,5,3,-10);
   fLayoutNormAdd = new TGLayoutHints(kLHintsLeft,10,5,3,-10);
   fLayoutConv    = new TGLayoutHints(kLHintsLeft,10,5,3,-10);

   bgr -> SetLayoutHints(fLayoutNone,   fNone);
   bgr -> SetLayoutHints(fLayoutAdd,    fAdd);
   bgr -> SetLayoutHints(fLayoutNormAdd,fNormAdd);
   bgr -> SetLayoutHints(fLayoutConv,   fConv);
   bgr -> Show();
   bgr -> ChangeOptions(kFitWidth | kHorizontalFrame);

   tf1 -> AddFrame(bgr, new TGLayoutHints(kLHintsExpandX, 0, 0, 3, 0));
   gf1 -> AddFrame(tf1, new TGLayoutHints(kLHintsExpandX));

   TGCompositeFrame *tf2 = new TGCompositeFrame(gf1, 350, 26,
                                                kHorizontalFrame);
   fEnteredFunc = new TGTextEntry(tf2, new TGTextBuffer(0), kFP_FILE);
   //fEnteredFunc->SetMaxLength(4000);  // use default value (~4000)
   fEnteredFunc->SetAlignment(kTextLeft);
   TGTextLBEntry *te = (TGTextLBEntry *)fFuncList->GetSelectedEntry();
   assert(te);
   fEnteredFunc->SetText(te->GetTitle());
   fEnteredFunc->SetToolTipText("Enter file_name/function_name or a function expression");
   fEnteredFunc->Resize(250,fEnteredFunc->GetDefaultHeight());
   tf2->AddFrame(fEnteredFunc, new TGLayoutHints(kLHintsLeft    |
                                                 kLHintsCenterY |
                                                 kLHintsExpandX, 2, 2, 2, 2));
   gf1->AddFrame(tf2, new TGLayoutHints(kLHintsNormal |
                                        kLHintsExpandX, 0, 0, 2, 0));

   TGHorizontalFrame *s1 = new TGHorizontalFrame(gf1);
   TGLabel *label21 = new TGLabel(s1, "Selected: ");
   s1->AddFrame(label21, new TGLayoutHints(kLHintsNormal |
                                           kLHintsCenterY, 2, 2, 2, 0));
   TGHorizontal3DLine *hlines = new TGHorizontal3DLine(s1);
   s1->AddFrame(hlines, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
   gf1->AddFrame(s1, new TGLayoutHints(kLHintsExpandX));

   TGCompositeFrame *tf4 = new TGCompositeFrame(gf1, 350, 26,
                                                kHorizontalFrame);
   TGTextLBEntry *txt = (TGTextLBEntry *)fFuncList->GetSelectedEntry();
   TString s = txt->GetTitle();
   fSelLabel = new TGLabel(tf4, s.Sizeof()>30?s(0,30)+"...":s);
   tf4->AddFrame(fSelLabel, new TGLayoutHints(kLHintsNormal |
                                              kLHintsCenterY, 0, 6, 2, 0));
   Pixel_t color;
   gClient->GetColorByName("#336666", color);
   fSelLabel->SetTextColor(color, kFALSE);
   TGCompositeFrame *tf5 = new TGCompositeFrame(tf4, 120, 20,
                                                kHorizontalFrame | kFixedWidth);
   fSetParam = new TGTextButton(tf5, "Set Parameters...", kFP_PARS);
   tf5->AddFrame(fSetParam, new TGLayoutHints(kLHintsRight   |
                                              kLHintsCenterY |
                                              kLHintsExpandX));
   fSetParam->SetToolTipText("Open a dialog for parameter(s) settings");
   tf4->AddFrame(tf5, new TGLayoutHints(kLHintsRight |
                                        kLHintsTop, 5, 0, 2, 2));
   gf1->AddFrame(tf4, new TGLayoutHints(kLHintsNormal |
                                             kLHintsExpandX, 5, 0, 0, 0));

   this->AddFrame(gf1, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

}

////////////////////////////////////////////////////////////////////////////////
/// Create 'General' tab.

void TFitEditor::CreateGeneralTab()
{
   fTabContainer = fTab->AddTab("General");
   fGeneral = new TGCompositeFrame(fTabContainer, 10, 10, kVerticalFrame);
   fTabContainer->AddFrame(fGeneral, new TGLayoutHints(kLHintsTop |
                                                       kLHintsExpandX,
                                                       5, 5, 2, 2));

   // 'options' group frame
   TGGroupFrame *gf = new TGGroupFrame(fGeneral, "Fit Settings", kFitWidth);

   // 'method' sub-group
   TGHorizontalFrame *h1 = new TGHorizontalFrame(gf);
   TGLabel *label4 = new TGLabel(h1, "Method");
   h1->AddFrame(label4, new TGLayoutHints(kLHintsNormal |
                                          kLHintsCenterY, 2, 2, 0, 0));
   TGHorizontal3DLine *hline1 = new TGHorizontal3DLine(h1);
   h1->AddFrame(hline1, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
   gf->AddFrame(h1, new TGLayoutHints(kLHintsExpandX));

   TGHorizontalFrame *h2 = new TGHorizontalFrame(gf);
   TGVerticalFrame *v1 = new TGVerticalFrame(h2);
   fMethodList = BuildMethodList(v1, kFP_MLIST);
   fMethodList->Select(1, kFALSE);
   fMethodList->Resize(140, 20);
   TGListBox *lb = fMethodList->GetListBox();
   Int_t lbe = lb->GetNumberOfEntries();
   lb->Resize(lb->GetWidth(), lbe*16);
   v1->AddFrame(fMethodList, new TGLayoutHints(kLHintsLeft, 0, 0, 2, 5));

   fLinearFit = new TGCheckButton(v1, "Linear fit", kFP_MLINF);
   fLinearFit->Associate(this);
   fLinearFit->SetToolTipText("Perform Linear fitter if selected");
   v1->AddFrame(fLinearFit, new TGLayoutHints(kLHintsNormal, 0, 0, 8, 2));


   h2->AddFrame(v1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   TGVerticalFrame *v2 = new TGVerticalFrame(h2);
   TGCompositeFrame *v21 = new TGCompositeFrame(v2, 120, 20,
                                                kHorizontalFrame | kFixedWidth);
   fUserButton = new TGTextButton(v21, "User-Defined...", kFP_MUSR);
   v21->AddFrame(fUserButton, new TGLayoutHints(kLHintsRight   |
                                                kLHintsCenterY |
                                                kLHintsExpandX));
   fUserButton->SetToolTipText("Open a dialog for entering a user-defined method");
   fUserButton->SetState(kButtonDisabled);
   v2->AddFrame(v21, new TGLayoutHints(kLHintsRight | kLHintsTop));

   TGHorizontalFrame *v1h = new TGHorizontalFrame(v2);
   fEnableRobust = new TGCheckButton(v1h, "Robust:", -1);
   fEnableRobust->Associate(this); // needed ???
   fEnableRobust->SetToolTipText("Perform Linear Robust fitter if selected");
   v1h->AddFrame(fEnableRobust, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));
   fRobustValue = new TGNumberEntry(v1h, 0.95, 5, kFP_RBUST,
                                    TGNumberFormat::kNESRealTwo,
                                    TGNumberFormat::kNEAPositive,
                                    TGNumberFormat::kNELLimitMinMax,0.,0.99);
   v1h->AddFrame(fRobustValue, new TGLayoutHints(kLHintsLeft));
   v2->AddFrame(v1h, new TGLayoutHints(kLHintsNormal, 0, 0, 12, 2));
   fRobustValue->SetState(kFALSE);
   fRobustValue->GetNumberEntry()->SetToolTipText("Available only for graphs");

   fNoChi2 = 0;
   // fNoChi2 = new TGCheckButton(v2, "No Chi-square", kFP_NOCHI);
   // fNoChi2->Associate(this);
   // fNoChi2->SetToolTipText("'C'- do not calculate Chi-square (for Linear fitter)");
   // v2->AddFrame(fNoChi2, new TGLayoutHints(kLHintsNormal, 0, 0, 34, 2));

   h2->AddFrame(v2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 20, 0, 0, 0));
   gf->AddFrame(h2, new TGLayoutHints(kLHintsExpandX, 20, 0, 0, 0));

   // 'fit option' sub-group
   TGHorizontalFrame *h3 = new TGHorizontalFrame(gf);
   TGLabel *label5 = new TGLabel(h3, "Fit Options");
   h3->AddFrame(label5, new TGLayoutHints(kLHintsNormal |
                                          kLHintsCenterY, 2, 2, 0, 0));
   TGHorizontal3DLine *hline2 = new TGHorizontal3DLine(h3);
   h3->AddFrame(hline2, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
   gf->AddFrame(h3, new TGLayoutHints(kLHintsExpandX));

   TGHorizontalFrame *h = new TGHorizontalFrame(gf);
   TGVerticalFrame *v3 = new TGVerticalFrame(h);
   fIntegral = new TGCheckButton(v3, "Integral", kFP_INTEG);
   fIntegral->Associate(this);
   fIntegral->SetToolTipText("'I'- use integral of function instead of value in bin center");
   v3->AddFrame(fIntegral, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fBestErrors = new TGCheckButton(v3, "Best errors", kFP_IMERR);
   fBestErrors->Associate(this);
   fBestErrors->SetToolTipText("'E'- better errors estimation using Minos technique");
   v3->AddFrame(fBestErrors, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fAllWeights1 = new TGCheckButton(v3, "All weights = 1", kFP_ALLW1);
   fAllWeights1->Associate(this);
   fAllWeights1->SetToolTipText("'W'- all weights=1 for non empty bins; error bars ignored");
   v3->AddFrame(fAllWeights1, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fEmptyBinsWghts1 = new TGCheckButton(v3, "Empty bins, weights=1", kFP_EMPW1);
   fEmptyBinsWghts1->Associate(this);
   fEmptyBinsWghts1->SetToolTipText("'WW'- all weights=1 including empty bins; error bars ignored");
   v3->AddFrame(fEmptyBinsWghts1, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   h->AddFrame(v3, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   TGVerticalFrame *v4 = new TGVerticalFrame(h);
   fUseRange = new TGCheckButton(v4, "Use range", kFP_USERG);
   fUseRange->Associate(this);
   fUseRange->SetToolTipText("'R'- fit only data within the specified function range");
   v4->AddFrame(fUseRange, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fImproveResults = new TGCheckButton(v4, "Improve fit results", kFP_IFITR);
   fImproveResults->Associate(this);
   fImproveResults->SetToolTipText("'M'- after minimum is found, search for a new one");
   v4->AddFrame(fImproveResults, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fAdd2FuncList = new TGCheckButton(v4, "Add to list", kFP_ADDLS);
   fAdd2FuncList->Associate(this);
   fAdd2FuncList->SetToolTipText("'+'- add function to the list without deleting the previous");
   v4->AddFrame(fAdd2FuncList, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fUseGradient = new TGCheckButton(v4, "Use Gradient", kFP_ADDLS);
   fUseGradient->Associate(this);
   fUseGradient->SetToolTipText("'G'- Use the gradient as an aid for the fitting");
   v4->AddFrame(fUseGradient, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   h->AddFrame(v4, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 20, 0, 0, 0));
   gf->AddFrame(h, new TGLayoutHints(kLHintsExpandX, 20, 0, 0, 0));

   // 'draw option' sub-group
   TGHorizontalFrame *h5 = new TGHorizontalFrame(gf);
   TGLabel *label6 = new TGLabel(h5, "Draw Options");
   h5->AddFrame(label6, new TGLayoutHints(kLHintsNormal |
                                          kLHintsCenterY, 2, 2, 2, 2));
   TGHorizontal3DLine *hline3 = new TGHorizontal3DLine(h5);
   h5->AddFrame(hline3, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
   gf->AddFrame(h5, new TGLayoutHints(kLHintsExpandX));

   TGHorizontalFrame *h6 = new TGHorizontalFrame(gf);
   TGVerticalFrame *v5 = new TGVerticalFrame(h6);

   fDrawSame = new TGCheckButton(v5, "SAME", kFP_DSAME);
   fDrawSame->Associate(this);
   fDrawSame->SetToolTipText("Superimpose on previous picture in the same pad");
   v5->AddFrame(fDrawSame, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fNoDrawing = new TGCheckButton(v5, "No drawing", kFP_DNONE);
   fNoDrawing->Associate(this);
   fNoDrawing->SetToolTipText("'0'- do not draw function graphics");
   v5->AddFrame(fNoDrawing, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   fNoStoreDrawing = new TGCheckButton(v5, "Do not store/draw", kFP_DNOST);
   fNoStoreDrawing->Associate(this);
   fNoStoreDrawing->SetToolTipText("'N'- do not store the function, do not draw it");
   v5->AddFrame(fNoStoreDrawing, new TGLayoutHints(kLHintsNormal, 0, 0, 2, 2));

   h6->AddFrame(v5, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   TGVerticalFrame *v6 = new TGVerticalFrame(h6);
   TGCompositeFrame *v61 = new TGCompositeFrame(v6, 120, 20,
                                                kHorizontalFrame | kFixedWidth);
   fDrawAdvanced = new TGTextButton(v61, "&Advanced...", kFP_DADVB);
   v61->AddFrame(fDrawAdvanced, new TGLayoutHints(kLHintsRight   |
                                                  kLHintsCenterY |
                                                  kLHintsExpandX));
   fDrawAdvanced->SetToolTipText("Open a dialog for advanced draw options");
   fDrawAdvanced->SetState(kButtonDisabled);

   v6->AddFrame(v61, new TGLayoutHints(kLHintsRight | kLHintsTop,
                                       0, 0, (4+fDrawSame->GetHeight())*2, 0));

   h6->AddFrame(v6, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   gf->AddFrame(h6, new TGLayoutHints(kLHintsExpandX, 20, 0, 2, 0));

   fGeneral->AddFrame(gf, new TGLayoutHints(kLHintsExpandX |
                                            kLHintsExpandY, 5, 5, 0, 0));
   // sliderX
   fSliderXParent = new TGHorizontalFrame(fGeneral);
   TGLabel *label8 = new TGLabel(fSliderXParent, "X");
   fSliderXParent->AddFrame(label8, new TGLayoutHints(kLHintsLeft |
                                                      kLHintsCenterY, 0, 5, 0, 0));

   fSliderXMin = new TGNumberEntry(fSliderXParent, 0, 5, kFP_XMIN,
                                   TGNumberFormat::kNESRealTwo,
                                   TGNumberFormat::kNEAAnyNumber,
                                   TGNumberFormat::kNELLimitMinMax, -1,1);
   fSliderXParent->AddFrame(fSliderXMin, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fSliderX = new TGDoubleHSlider(fSliderXParent, 1, kDoubleScaleBoth);
   fSliderX->SetScale(5);
   fSliderXParent->AddFrame(fSliderX, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));


   fSliderXMax = new TGNumberEntry(fSliderXParent, 0, 5, kFP_XMIN,
                                   TGNumberFormat::kNESRealTwo,
                                   TGNumberFormat::kNEAAnyNumber,
                                   TGNumberFormat::kNELLimitMinMax, -1,1);
   fSliderXParent->AddFrame(fSliderXMax, new TGLayoutHints(kLHintsRight | kLHintsCenterY));
   fGeneral->AddFrame(fSliderXParent, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   // sliderY
   fSliderYParent = new TGHorizontalFrame(fGeneral);
   TGLabel *label9 = new TGLabel(fSliderYParent, "Y");
   fSliderYParent->AddFrame(label9, new TGLayoutHints(kLHintsLeft |
                                                      kLHintsCenterY, 0, 5, 0, 0));

   fSliderYMin = new TGNumberEntry(fSliderYParent, 0, 5, kFP_YMIN,
                                   TGNumberFormat::kNESRealTwo,
                                   TGNumberFormat::kNEAAnyNumber,
                                   TGNumberFormat::kNELLimitMinMax, -1,1);
   fSliderYParent->AddFrame(fSliderYMin, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fSliderY = new TGDoubleHSlider(fSliderYParent, 1, kDoubleScaleBoth);
   fSliderY->SetScale(5);
   fSliderYParent->AddFrame(fSliderY, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

   fSliderYMax = new TGNumberEntry(fSliderYParent, 0, 5, kFP_YMIN,
                                   TGNumberFormat::kNESRealTwo,
                                   TGNumberFormat::kNEAAnyNumber,
                                   TGNumberFormat::kNELLimitMinMax, -1,1);
   fSliderYParent->AddFrame(fSliderYMax, new TGLayoutHints(kLHintsRight | kLHintsCenterY));
   fGeneral->AddFrame(fSliderYParent, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   // sliderZ
   fSliderZParent = new TGHorizontalFrame(fGeneral);
   TGLabel *label10 = new TGLabel(fSliderZParent, "Z:");
   fSliderZParent->AddFrame(label10, new TGLayoutHints(kLHintsLeft |
                                                       kLHintsCenterY, 0, 5, 0, 0));
   fSliderZ = new TGDoubleHSlider(fSliderZParent, 1, kDoubleScaleBoth);
   fSliderZ->SetScale(5);
   fSliderZParent->AddFrame(fSliderZ, new TGLayoutHints(kLHintsExpandX |
                                                        kLHintsCenterY));
   fGeneral->AddFrame(fSliderZParent, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));
}


////////////////////////////////////////////////////////////////////////////////
/// Create 'Minimization' tab.

void TFitEditor::CreateMinimizationTab()
{
   fTabContainer = fTab->AddTab("Minimization");
   fMinimization = new TGCompositeFrame(fTabContainer, 10, 10, kVerticalFrame);
   fTabContainer->AddFrame(fMinimization, new TGLayoutHints(kLHintsTop |
                                                            kLHintsExpandX,
                                                            5, 5, 2, 2));
   MakeTitle(fMinimization, "Library");

   TGHorizontalFrame *hl = new TGHorizontalFrame(fMinimization);
   fLibMinuit = new TGRadioButton(hl, "Minuit", kFP_LMIN);
   fLibMinuit->Associate(this);
   fLibMinuit->SetToolTipText("Use minimization from libMinuit (default)");
   hl->AddFrame(fLibMinuit, new TGLayoutHints(kLHintsNormal, 40, 0, 0, 1));
   fStatusBar->SetText("LIB Minuit",1);

   fLibMinuit2 = new TGRadioButton(hl, "Minuit2", kFP_LMIN2);
   fLibMinuit2->Associate(this);
   fLibMinuit2->SetToolTipText("New C++ version of Minuit");
   hl->AddFrame(fLibMinuit2, new TGLayoutHints(kLHintsNormal, 35, 0, 0, 1));

   fLibFumili = new TGRadioButton(hl, "Fumili", kFP_LFUM);
   fLibFumili->Associate(this);
   fLibFumili->SetToolTipText("Use minimization from libFumili");
   hl->AddFrame(fLibFumili, new TGLayoutHints(kLHintsNormal, 30, 0, 0, 1));
   fMinimization->AddFrame(hl, new TGLayoutHints(kLHintsExpandX, 20, 0, 5, 1));

   TGHorizontalFrame *hl2 = new TGHorizontalFrame(fMinimization);

   fLibGSL = new TGRadioButton(hl2, "GSL", kFP_LGSL);
   #ifdef R__HAS_MATHMORE
   fLibGSL->Associate(this);
   fLibGSL->SetToolTipText("Use minimization from libGSL");
   #else
   fLibGSL->SetState(kButtonDisabled);
   fLibGSL->SetToolTipText("Needs GSL to be compiled");
   #endif
   hl2->AddFrame(fLibGSL, new TGLayoutHints(kLHintsNormal, 40, 0, 0, 1));

   fLibGenetics = new TGRadioButton(hl2, "Genetics", kFP_LGAS);
   if (gPluginMgr->FindHandler("ROOT::Math::Minimizer","Genetic") ||
       gPluginMgr->FindHandler("ROOT::Math::Minimizer","GAlibMin") )
   {
      fLibGenetics->Associate(this);
      fLibGenetics->SetToolTipText("Different GAs implementations");
   } else {
      fLibGenetics->SetState(kButtonDisabled);
      fLibGenetics->SetToolTipText("Needs any of the genetic"
                                   "minimizers to be compiled");
   }
   hl2->AddFrame(fLibGenetics, new TGLayoutHints(kLHintsNormal, 45, 0, 0, 1));

   fMinimization->AddFrame(hl2, new TGLayoutHints(kLHintsExpandX, 20, 0, 5, 1));

   MakeTitle(fMinimization, "Method");

   TGHorizontalFrame *hm0 = new TGHorizontalFrame(fMinimization);
   fMinMethodList = new TGComboBox(hm0, kFP_MINMETHOD);
   fMinMethodList->Resize(290, 20);
   fMinMethodList->Select(kFP_GAUS, kFALSE);

   TGListBox *lb = fMinMethodList->GetListBox();
   lb->Resize(lb->GetWidth(), 500);
   fMinMethodList->Associate(this);

   hm0->AddFrame(fMinMethodList, new TGLayoutHints(kLHintsNormal));
   fMinimization->AddFrame(hm0, new TGLayoutHints(kLHintsExpandX, 60, 0, 5, 1));

   // Set the status to the default minimization options!
   if ( ROOT::Math::MinimizerOptions::DefaultMinimizerType() == "Fumili" ) {
      fLibFumili->SetState(kButtonDown);
   } else if ( ROOT::Math::MinimizerOptions::DefaultMinimizerType() == "Minuit" ) {
      fLibMinuit->SetState(kButtonDown);
   } else {
      fLibMinuit2->SetState(kButtonDown);
   }
   FillMinMethodList();

   MakeTitle(fMinimization, "Settings");
   TGLabel *hslabel1 = new TGLabel(fMinimization,"Use ENTER key to validate a new value or click");
   fMinimization->AddFrame(hslabel1, new TGLayoutHints(kLHintsNormal, 61, 0, 5, 1));
   TGLabel *hslabel2 = new TGLabel(fMinimization,"on Reset button to set the defaults.");
   fMinimization->AddFrame(hslabel2, new TGLayoutHints(kLHintsNormal, 61, 0, 1, 10));

   TGHorizontalFrame *hs = new TGHorizontalFrame(fMinimization);

   TGVerticalFrame *hsv1 = new TGVerticalFrame(hs, 180, 10, kFixedWidth);
   TGLabel *errlabel = new TGLabel(hsv1,"Error definition (default = 1): ");
   hsv1->AddFrame(errlabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                              1, 1, 5, 7));
   TGLabel *tollabel = new TGLabel(hsv1,"Max tolerance (precision): ");
   hsv1->AddFrame(tollabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                              1, 1, 5, 7));
   TGLabel *itrlabel = new TGLabel(hsv1,"Max number of iterations: ");
   hsv1->AddFrame(itrlabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,
                                              1, 1, 5, 5));
   hs->AddFrame(hsv1, new TGLayoutHints(kLHintsNormal, 60, 0, 0, 0));

   TGVerticalFrame *hsv2 = new TGVerticalFrame(hs, 90,10, kFixedWidth);
   fErrorScale = new TGNumberEntryField(hsv2, kFP_MERR, ROOT::Math::MinimizerOptions::DefaultErrorDef(),
                                        TGNumberFormat::kNESRealTwo,
                                        TGNumberFormat::kNEAPositive,
                                        TGNumberFormat::kNELLimitMinMax,0.,100.);
   hsv2->AddFrame(fErrorScale, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,
                                                 1, 1, 0, 3));
   fTolerance = new TGNumberEntryField(hsv2, kFP_MTOL, ROOT::Math::MinimizerOptions::DefaultTolerance(),
                                       TGNumberFormat::kNESReal,
                                       TGNumberFormat::kNEAPositive,
                                       TGNumberFormat::kNELLimitMinMax, 0., 1.);
   fTolerance->SetNumber(ROOT::Math::MinimizerOptions::DefaultTolerance());
   hsv2->AddFrame(fTolerance, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,
                                                1, 1, 3, 3));
   fIterations = new TGNumberEntryField(hsv2, kFP_MITR, 5000,
                                   TGNumberFormat::kNESInteger,
                                   TGNumberFormat::kNEAPositive,
                                   TGNumberFormat::kNELNoLimits);
   fIterations->SetNumber(ROOT::Math::MinimizerOptions::DefaultMaxIterations());
   hsv2->AddFrame(fIterations, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,
                                                 1, 1, 3, 3));
   hs->AddFrame(hsv2, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   fMinimization->AddFrame(hs, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 1));
   fStatusBar->SetText(Form("Itr: %d",ROOT::Math::MinimizerOptions::DefaultMaxIterations()),3);

   MakeTitle(fMinimization, "Print Options");

   TGHorizontalFrame *h8 = new TGHorizontalFrame(fMinimization);
   fOptDefault = new TGRadioButton(h8, "Default", kFP_PDEF);
   fOptDefault->Associate(this);
   fOptDefault->SetToolTipText("Default is between Verbose and Quiet");
   h8->AddFrame(fOptDefault, new TGLayoutHints(kLHintsNormal, 40, 0, 0, 1));
   fOptDefault->SetState(kButtonDown);
   fStatusBar->SetText("Prn: DEF",4);

   fOptVerbose = new TGRadioButton(h8, "Verbose", kFP_PVER);
   fOptVerbose->Associate(this);
   fOptVerbose->SetToolTipText("'V'- print results after each iteration");
   h8->AddFrame(fOptVerbose, new TGLayoutHints(kLHintsNormal, 30, 0, 0, 1));

   fOptQuiet = new TGRadioButton(h8, "Quiet", kFP_PQET);
   fOptQuiet->Associate(this);
   fOptQuiet->SetToolTipText("'Q'- no print");
   h8->AddFrame(fOptQuiet, new TGLayoutHints(kLHintsNormal, 25, 0, 0, 1));

   fMinimization->AddFrame(h8, new TGLayoutHints(kLHintsExpandX, 20, 0, 5, 1));

}

////////////////////////////////////////////////////////////////////////////////
/// Connect GUI signals to fit panel slots.

void TFitEditor::ConnectSlots()
{
   // list of data sets to fit
   fDataSet          -> Connect("Selected(Int_t)", "TFitEditor", this, "DoDataSet(Int_t)");
   // list of predefined functions
   fTypeFit          -> Connect("Selected(Int_t)", "TFitEditor", this, "FillFunctionList(Int_t)");
   // list of predefined functions
   fFuncList         -> Connect("Selected(Int_t)", "TFitEditor", this, "DoFunction(Int_t)");
   // entered formula or function name
   fEnteredFunc      -> Connect("ReturnPressed()", "TFitEditor", this, "DoEnteredFunction()");
   // set parameters dialog
   fSetParam         -> Connect("Clicked()",       "TFitEditor", this, "DoSetParameters()");
   // allowed function operations
   fAdd              -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoAddition(Bool_t)");
   //fNormAdd          -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoNormAddition(Bool_t)");
   //fConv             -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoConvolution(Bool_t)");
   // fit options
   fAllWeights1      -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoAllWeights1()");
   fUseRange         -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoUseFuncRange()");
   fEmptyBinsWghts1  -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoEmptyBinsAllWeights1()");
   // linear fit
   fLinearFit        -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLinearFit()");
   fEnableRobust     -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoRobustFit()");
   //fNoChi2->Connect("Toggled(Bool_t)","TFitEditor",this,"DoNoChi2()");
   // draw options
   fNoStoreDrawing   -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoNoStoreDrawing()");
   // fit, reset, close buttons
   fUpdateButton     -> Connect("Clicked()",       "TFitEditor", this, "DoUpdate()");
   fFitButton        -> Connect("Clicked()",       "TFitEditor", this, "DoFit()");
   fResetButton      -> Connect("Clicked()",       "TFitEditor", this, "DoReset()");
   fCloseButton      -> Connect("Clicked()",       "TFitEditor", this, "DoClose()");
   // user method button
   fUserButton       -> Connect("Clicked()",       "TFitEditor", this, "DoUserDialog()");
   // advanced draw options
   fDrawAdvanced     -> Connect("Clicked()",       "TFitEditor", this, "DoAdvancedOptions()");

   if (fType != kObjectTree)
   {
      fSliderX       -> Connect("PositionChanged()","TFitEditor",this, "DoSliderXMoved()");
      fSliderXMax    -> Connect("ValueSet(Long_t)", "TFitEditor",this, "DoNumericSliderXChanged()");
      fSliderXMin    -> Connect("ValueSet(Long_t)", "TFitEditor",this, "DoNumericSliderXChanged()");
   }
   if (fDim > 1)
   {
      fSliderY       -> Connect("PositionChanged()","TFitEditor",this, "DoSliderYMoved()");
      fSliderYMax    -> Connect("ValueSet(Long_t)", "TFitEditor",this, "DoNumericSliderYChanged()");
      fSliderYMin    -> Connect("ValueSet(Long_t)", "TFitEditor",this, "DoNumericSliderYChanged()");
   }
   if (fDim > 2)
      fSliderZ       -> Connect("PositionChanged()","TFitEditor",this, "DoSliderZMoved()");

   if ( fParentPad )
      fParentPad     -> Connect("RangeAxisChanged()","TFitEditor",this, "UpdateGUI()");
   // 'Minimization' tab
   // library
   fLibMinuit        -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLibrary(Bool_t)");
   fLibMinuit2       -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLibrary(Bool_t)");
   fLibFumili        -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLibrary(Bool_t)");
   fLibGSL           -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLibrary(Bool_t)");
   fLibGenetics      -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoLibrary(Bool_t)");

   // minimization method
   fMinMethodList    -> Connect("Selected(Int_t)", "TFitEditor", this, "DoMinMethod(Int_t)");
   // fitter settings
   fIterations       -> Connect("ReturnPressed()", "TFitEditor", this, "DoMaxIterations()");
   // print options
   fOptDefault       -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoPrintOpt(Bool_t)");
   fOptVerbose       -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoPrintOpt(Bool_t)");
   fOptQuiet         -> Connect("Toggled(Bool_t)", "TFitEditor", this, "DoPrintOpt(Bool_t)");

}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect GUI signals from fit panel slots.

void TFitEditor::DisconnectSlots()
{
   Disconnect("CloseWindow()");

   fFuncList    -> Disconnect("Selected(Int_t)");
   fEnteredFunc -> Disconnect("ReturnPressed()");
   fSetParam    -> Disconnect("Clicked()");
   fAdd         -> Disconnect("Toggled(Bool_t)");
   // fNormAdd     -> Disconnect("Toggled(Bool_t)");
   // fConv        -> Disconnect("Toggled(Bool_t)");

   // fit options
   fAllWeights1      -> Disconnect("Toggled(Bool_t)");
   fEmptyBinsWghts1  -> Disconnect("Toggled(Bool_t)");

   // linear fit
   fLinearFit        -> Disconnect("Toggled(Bool_t)");
   fEnableRobust     -> Disconnect("Toggled(Bool_t)");
   //fNoChi2->Disconnect("Toggled(Bool_t)");

   // draw options
   fNoStoreDrawing -> Disconnect("Toggled(Bool_t)");

   // fit, reset, close buttons
   fFitButton     -> Disconnect("Clicked()");
   fResetButton   -> Disconnect("Clicked()");

   // other methods
   fUserButton    -> Disconnect("Clicked()");
   fDrawAdvanced  -> Disconnect("Clicked()");

   if (fType != kObjectTree)
   {
      fSliderX    -> Disconnect("PositionChanged()");
      fSliderXMax -> Disconnect("ValueChanged(Long_t)");
      fSliderXMin -> Disconnect("ValueChanged(Long_t)");
   }
   if (fDim > 1)
   {
      fSliderY    -> Disconnect("PositionChanged()");
      fSliderYMax -> Disconnect("ValueChanged(Long_t)");
      fSliderYMin -> Disconnect("ValueChanged(Long_t)");
   }
   if (fDim > 2)
      fSliderZ    -> Disconnect("PositionChanged()");
   // slots related to 'Minimization' tab
   fLibMinuit     -> Disconnect("Toggled(Bool_t)");
   fLibMinuit2    -> Disconnect("Toggled(Bool_t)");
   fLibFumili     -> Disconnect("Toggled(Bool_t)");
   fLibGSL        -> Disconnect("Toggled(Bool_t)");
   fLibGenetics   -> Disconnect("Toggled(Bool_t)");
   // minimization method
   fMinMethodList -> Disconnect("Selected(Int_t)");
   // fitter settings
   fIterations    -> Disconnect("ReturnPressed()");
   // print options
   fOptDefault    -> Disconnect("Toggled(Bool_t)");
   fOptVerbose    -> Disconnect("Toggled(Bool_t)");
   fOptQuiet      -> Disconnect("Toggled(Bool_t)");

}

////////////////////////////////////////////////////////////////////////////////
/// Connect to another canvas.

void TFitEditor::SetCanvas(TCanvas * /*newcan*/)
{
   // The next line is commented because it is stablishing a
   // connection with the particular canvas, while right the following
   // line will connect all the canvas in a general way.

   // It would also make the fitpanel crash if there is no object
   // defined to be fitted in the construction (as a side effect of
   // it).

//    newcan->Connect("Selected(TVirtualPad*,TObject*,Int_t)", "TFitEditor",
//                    this, "SetFitObject(TVirtualPad *, TObject *, Int_t)");

   TQObject::Connect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)",
                     "TFitEditor",this,
                     "SetFitObject(TVirtualPad *, TObject *, Int_t)");
   TQObject::Connect("TCanvas", "Closed()", "TFitEditor", this, "DoNoSelection()");
}

////////////////////////////////////////////////////////////////////////////////
/// Hide the fit panel and set it to non-active state.

void TFitEditor::Hide()
{
   if (fgFitDialog) {
      fgFitDialog->UnmapWindow();
   }
   if (fParentPad) {
      fParentPad->Disconnect("RangeAxisChanged()");
      DoReset();
      TQObject::Disconnect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)",
                           this, "SetFitObject(TVirtualPad *, TObject *, Int_t)");
   }
   fParentPad = 0;
   fFitObject = 0;
   gROOT->GetListOfCleanups()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Show the fit panel (possible only via context menu).

void TFitEditor::Show(TVirtualPad* pad, TObject *obj)
{
   if (!gROOT->GetListOfCleanups()->FindObject(this))
      gROOT->GetListOfCleanups()->Add(this);

   if (!fgFitDialog->IsMapped()) {
      fgFitDialog->MapWindow();
      gVirtualX->RaiseWindow(GetId());
   }
   fParentPad = static_cast<TPad*>(pad);
   SetCanvas(pad->GetCanvas());
   SetFitObject(pad, obj, kButton1Down);
}

////////////////////////////////////////////////////////////////////////////////
/// Close fit panel window.

void TFitEditor::CloseWindow()
{
   Hide();
}

//______________________________________________________________________________
// TFitEditor *&TFitEditor::GetFP()
// {
//    // Static: return main fit panel
//    return fgFitDialog;
// }

////////////////////////////////////////////////////////////////////////////////
///  Called to delete the fit panel.

void TFitEditor::Terminate()
{
   TQObject::Disconnect("TCanvas", "Closed()");
   delete fgFitDialog;
   fgFitDialog = 0;
}

////////////////////////////////////////////////////////////////////////////////
///  Set the fit panel GUI according to the selected object.

void TFitEditor::UpdateGUI()
{
   if (!fFitObject) return;

   DrawSelection(true);

   if ( fType == kObjectTree )
      // Don't do anything with the sliders, as they work with TAxis
      // that are not defined for the TTree
      return;

   // sliders
   if (fType != kObjectTree) { // This is as fDim > 0
      TH1* hist = 0;
      switch (fType) {
         case kObjectHisto:
            hist = (TH1*)fFitObject;
            break;

         case kObjectGraph:
            hist = ((TGraph*)fFitObject)->GetHistogram();
            break;

         case kObjectMultiGraph:
            hist = ((TMultiGraph*)fFitObject)->GetHistogram();
            break;

         case kObjectGraph2D:
            hist = ((TGraph2D*)fFitObject)->GetHistogram("empty");
            break;

         case kObjectHStack:
            hist = (TH1 *)((THStack *)fFitObject)->GetHists()->First();

         case kObjectTree:
         default:
            break;
      }


      if (!hist) {
         Error("UpdateGUI","No hist is present - this should not happen, please report."
               "The FitPanel might be in an inconsistent state");
         //assert(hist);
         return;
      }

      fSliderX->Disconnect("PositionChanged()");
      fSliderXMin->Disconnect("ValueChanged()");
      fSliderXMax->Disconnect("ValueChanged()");

      if (!fSliderXParent->IsMapped())
         fSliderXParent->MapWindow();

      fXaxis = hist->GetXaxis();
      fYaxis = hist->GetYaxis();
      fZaxis = hist->GetZaxis();
      Int_t ixrange = fXaxis->GetNbins();
      Int_t ixmin = fXaxis->GetFirst();
      Int_t ixmax = fXaxis->GetLast();

      if (ixmin > 1 || ixmax < ixrange) {
         fSliderX->SetRange(ixmin,ixmax);
         fSliderX->SetPosition(ixmin, ixmax);
      } else {
         fSliderX->SetRange(1,ixrange);
         fSliderX->SetPosition(ixmin,ixmax);
      }

      fSliderX->SetScale(5);

      fSliderXMin->SetLimits(TGNumberFormat::kNELLimitMinMax,
                             fXaxis->GetBinLowEdge( static_cast<Int_t>( fSliderX->GetMinPosition() ) ),
                             fXaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderX->GetMaxPosition() ) ));
      fSliderXMin->SetNumber( fXaxis->GetBinLowEdge( static_cast<Int_t>( fSliderX->GetMinPosition() ) ));
      fSliderXMax->SetLimits(TGNumberFormat::kNELLimitMinMax,
                             fXaxis->GetBinLowEdge( static_cast<Int_t>( fSliderX->GetMinPosition() ) ),
                             fXaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderX->GetMaxPosition() ) ));
      fSliderXMax->SetNumber( fXaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderX->GetMaxPosition() ) ));

      fSliderX->Connect("PositionChanged()","TFitEditor",this, "DoSliderXMoved()");
      fSliderXMax->Connect("ValueSet(Long_t)", "TFitEditor", this, "DoNumericSliderXChanged()");
      fSliderXMin->Connect("ValueSet(Long_t)", "TFitEditor", this, "DoNumericSliderXChanged()");
   }

   if (fDim > 1) {
      fSliderY->Disconnect("PositionChanged()");
      fSliderYMin->Disconnect("ValueChanged()");
      fSliderYMax->Disconnect("ValueChanged()");

      if (!fSliderYParent->IsMapped())
         fSliderYParent->MapWindow();
      if (fSliderZParent->IsMapped())
         fSliderZParent->UnmapWindow();

      Int_t iymin = 0, iymax = 0, iyrange = 0;
      switch (fType) {
         case kObjectHisto:
         case kObjectGraph2D:
         case kObjectHStack:
            iyrange = fYaxis->GetNbins();
            iymin = fYaxis->GetFirst();
            iymax = fYaxis->GetLast();
            break;

         case kObjectGraph:
         case kObjectMultiGraph:
         case kObjectTree:
         default:
            //not implemented
            break;
      }

      if (iymin > 1 || iymax < iyrange) {
         fSliderY->SetRange(iymin,iymax);
         fSliderY->SetPosition(iymin, iymax);
      } else {
         fSliderY->SetRange(1,iyrange);
         fSliderY->SetPosition(iymin,iymax);
      }

      fSliderY->SetScale(5);

      fSliderYMin->SetLimits(TGNumberFormat::kNELLimitMinMax,
                             fYaxis->GetBinLowEdge( static_cast<Int_t>( fSliderY->GetMinPosition() ) ),
                             fYaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderY->GetMaxPosition() ) ));
      fSliderYMin->SetNumber(fYaxis->GetBinLowEdge( static_cast<Int_t>( fSliderY->GetMinPosition() ) ));
      fSliderYMax->SetLimits(TGNumberFormat::kNELLimitMinMax,
                             fYaxis->GetBinLowEdge( static_cast<Int_t>( fSliderY->GetMinPosition() ) ),
                             fYaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderY->GetMaxPosition() ) ));
      fSliderYMax->SetNumber( fYaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderY->GetMaxPosition() ) ));

      fSliderY->Connect("PositionChanged()","TFitEditor",this, "DoSliderYMoved()");
      fSliderYMax->Connect("ValueSet(Long_t)", "TFitEditor", this, "DoNumericSliderYChanged()");
      fSliderYMin->Connect("ValueSet(Long_t)", "TFitEditor", this, "DoNumericSliderYChanged()");
   }


   if (fDim > 2) {
      fSliderZ->Disconnect("PositionChanged()");

      if (!fSliderZParent->IsMapped())
         fSliderZParent->MapWindow();

      Int_t izmin = 0, izmax = 0, izrange = 0;
      switch (fType) {
         case kObjectHStack:
         case kObjectHisto:
            izrange = fZaxis->GetNbins();
            izmin = fZaxis->GetFirst();
            izmax = fZaxis->GetLast();
            break;

         case kObjectGraph:
         case kObjectGraph2D:
         case kObjectMultiGraph:
         case kObjectTree:
         default:
            //not implemented
            break;
      }

      if (izmin > 1 || izmax < izrange) {
         fSliderZ->SetRange(izmin,izmax);
         fSliderZ->SetPosition(izmin, izmax);
      } else {
         fSliderZ->SetRange(1,izrange);
         fSliderZ->SetPosition(izmin,izmax);
      }

      fSliderZ->SetScale(5);
      fSliderZ->Connect("PositionChanged()","TFitEditor",this, "DoSliderZMoved()");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot called when the user clicks on an object inside a canvas.
/// Updates pointers to the parent pad and the selected object
/// for fitting (if suitable).

void TFitEditor::SetFitObject(TVirtualPad *pad, TObject *obj, Int_t event)
{
   if (event != kButton1Down) return;

   if ( !obj ) {
      DoNoSelection();
      return;
   }

   // is obj suitable for fitting?
   if (!SetObjectType(obj)) return;

   fParentPad = pad;
   fFitObject = obj;
   ShowObjectName(obj);
   UpdateGUI();

   ConnectSlots();

   TF1* fitFunc = HasFitFunction();

   if (fitFunc)
   {
      //fFuncPars = FuncParams_t( fitFunc->GetNpar() );
      GetParameters(fFuncPars, fitFunc);

      TString tmpStr = fitFunc->GetExpFormula();
      TGLBEntry *en = 0;
      // If the function comes from a C raw function.
      if ( tmpStr.Length() == 0 )
      {
         // Show the name of the function
         fEnteredFunc->SetText(fitFunc->GetName());
         en= fFuncList->FindEntry(fitFunc->GetName());
         // Don't allow edition!
         SetEditable(kFALSE);
      }
      // otherwise, it's got a formula
      else
      {
         // Show the formula
         fEnteredFunc->SetText(fitFunc->GetExpFormula().Data());
         en= fFuncList->FindEntry(fitFunc->GetExpFormula().Data());
         SetEditable(kTRUE);
      }
      // Select the proper entry in the function list
      if (en) fFuncList->Select(en->EntryId());
   }
   else
   { // if there is no fit function in the object
      // Use the selected function in fFuncList
      TGTextLBEntry *te = (TGTextLBEntry *)fFuncList->GetSelectedEntry();
      // Add the text to fEnteredFunc
      if     (te && fNone->GetState() == kButtonDown)
         fEnteredFunc->SetText(te->GetTitle());
      else if (te && fAdd->GetState() == kButtonDown)
      {
         TString tmpStr = fEnteredFunc->GetText();
         tmpStr += '+';
         tmpStr += te->GetTitle();
         fEnteredFunc->SetText(tmpStr);
      }
      else if (te && fNormAdd->GetState() == kButtonDown)
      {
         TString tmpStr = fEnteredFunc->GetText();
         tmpStr += '+';
         tmpStr += te -> GetTitle();
         fEnteredFunc -> SetText(tmpStr);
      }
      else if (te && fConv->GetState() == kButtonDown)
      {
         TString tmpStr = fEnteredFunc->GetText();
         tmpStr += '*';
         tmpStr +=te->GetTitle();
         fEnteredFunc->SetText(tmpStr);
      }
      else if ( !te )
         // If there is no space, an error message is shown:
         // Error in <TString::AssertElement>: out of bounds: i = -1, Length = 0
         // If there is no function selected, then put nothing.
         fEnteredFunc->SetText(" ");
   }
   fEnteredFunc->SelectAll();


   // Update the information about the selected object.
   if (fSetParam->GetState() == kButtonDisabled)
      fSetParam->SetEnabled(kTRUE);
   if (fFitButton->GetState() == kButtonDisabled)
      fFitButton->SetEnabled(kTRUE);
   if (fResetButton->GetState() == kButtonDisabled)
      fResetButton->SetEnabled(kTRUE);
   DoLinearFit();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot called when users close a TCanvas or when the user select
/// no object.

void TFitEditor::DoNoSelection()
{
   if (gROOT->GetListOfCanvases()->IsEmpty()) {
      Terminate();
      return;
   }

   // Minimize user interaction until an object is selected
   DisconnectSlots();
   fParentPad = 0;
   fFitObject = 0;
   fStatusBar->SetText("No selection",0);
   fDataSet->Select(kFP_NOSEL, kFALSE);
   Layout();

   fSetParam->SetEnabled(kFALSE);
   fFitButton->SetEnabled(kFALSE);
   fResetButton->SetEnabled(kFALSE);
   fDrawAdvanced->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// When obj is deleted, clear fFitObject if fFitObject = obj.

void TFitEditor::RecursiveRemove(TObject* obj)
{
   if (obj == fFitObject) {
      fFitObject = 0;
      DisconnectSlots();
      fStatusBar->SetText("No selection",0);
      fDataSet->Select(kFP_NOSEL, kFALSE);
      Layout();

      fFitButton->SetEnabled(kFALSE);
      fResetButton->SetEnabled(kFALSE);
      fSetParam->SetEnabled(kFALSE);

      TQObject::Connect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)",
                        "TFitEditor",this,
                        "SetFitObject(TVirtualPad *, TObject *, Int_t)");
      TQObject::Connect("TCanvas", "Closed()", "TFitEditor", this,
                        "DoNoSelection()");

      DoUpdate();
      return;
   }
   if (obj == fParentPad) {
      fFitObject = 0;
      fParentPad = 0;
      DisconnectSlots();
      fStatusBar->SetText("No selection",0);
      fDataSet->Select(kFP_NOSEL, kFALSE);
      Layout();

      fFitButton->SetEnabled(kFALSE);
      fResetButton->SetEnabled(kFALSE);
      fSetParam->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the list of functions depending on the type of fit
/// selected.

void TFitEditor::FillFunctionList(Int_t)
{
   fFuncList->RemoveAll();
   // Case when the user has selected predefined functions in 1D.
   if ( fTypeFit->GetSelected() == kFP_PRED1D && fDim <= 1 ) {
      // Fill function list combo box.
      fFuncList->AddEntry("gaus" ,  kFP_GAUS);
      fFuncList->AddEntry("gausn",  kFP_GAUSN);
      fFuncList->AddEntry("expo",   kFP_EXPO);
      fFuncList->AddEntry("landau", kFP_LAND);
      fFuncList->AddEntry("landaun",kFP_LANDN);
      fFuncList->AddEntry("pol0",   kFP_POL0);
      fFuncList->AddEntry("pol1",   kFP_POL1);
      fFuncList->AddEntry("pol2",   kFP_POL2);
      fFuncList->AddEntry("pol3",   kFP_POL3);
      fFuncList->AddEntry("pol4",   kFP_POL4);
      fFuncList->AddEntry("pol5",   kFP_POL5);
      fFuncList->AddEntry("pol6",   kFP_POL6);
      fFuncList->AddEntry("pol7",   kFP_POL7);
      fFuncList->AddEntry("pol8",   kFP_POL8);
      fFuncList->AddEntry("pol9",   kFP_POL9);
      fFuncList->AddEntry("cheb0",   kFP_CHEB0);
      fFuncList->AddEntry("cheb1",   kFP_CHEB1);
      fFuncList->AddEntry("cheb2",   kFP_CHEB2);
      fFuncList->AddEntry("cheb3",   kFP_CHEB3);
      fFuncList->AddEntry("cheb4",   kFP_CHEB4);
      fFuncList->AddEntry("cheb5",   kFP_CHEB5);
      fFuncList->AddEntry("cheb6",   kFP_CHEB6);
      fFuncList->AddEntry("cheb7",   kFP_CHEB7);
      fFuncList->AddEntry("cheb8",   kFP_CHEB8);
      fFuncList->AddEntry("cheb9",   kFP_CHEB9);
      fFuncList->AddEntry("user",   kFP_USER);

      // Need to be setted this way, otherwise when the functions
      // are removed, the list doesn't show them.
      TGListBox *lb = fFuncList->GetListBox();
      lb->Resize(lb->GetWidth(), 200);

      // Select Gaus1D by default
      fFuncList->Select(kFP_GAUS);

   }
   // Case for predefined 2D functions
   else if ( fTypeFit->GetSelected() == kFP_PRED2D && fDim == 2 ) {
      fFuncList->AddEntry("xygaus", kFP_XYGAUS);
      fFuncList->AddEntry("bigaus", kFP_BIGAUS);
      fFuncList->AddEntry("xyexpo", kFP_XYEXP);
      fFuncList->AddEntry("xylandau", kFP_XYLAN);
      fFuncList->AddEntry("xylandaun", kFP_XYLANN);

      // Need to be setted this way, otherwise when the functions
      // are removed, the list doesn't show them.x
      TGListBox *lb = fFuncList->GetListBox();
      lb->Resize(lb->GetWidth(), 200);

      // Select Gaus2D by default
      fFuncList->Select(kFP_XYGAUS);
   }
   // Case for user defined functions. References to these functions
   // are kept by the fitpanel, so the information is gathered from
   // there.
   else if ( fTypeFit->GetSelected() == kFP_UFUNC ) {
      Int_t newid = kFP_ALTFUNC;

      // Add system functions
      for (auto f : fSystemFuncs) {
         // Don't include system functions that has been previously
         // used to fit, as those are included under the kFP_PREVFIT
         // section.
         if ( strncmp(f->GetName(), "PrevFit", 7) != 0 ) {
            // If the dimension of the object coincides with the
            // dimension of the function, then include the function in
            // the list. It will also include de function if the
            // dimension of the object is 0 (i.e. a multivariable
            // TTree) as it is currently imposible to know how many
            // dimensions a TF1 coming from a C raw function has.
            if ( f->GetNdim() == fDim || fDim == 0) {
               fFuncList->AddEntry(f->GetName(), newid++);
            }
         }
      }

      // If no function was added
      if ( newid != kFP_ALTFUNC )
         fFuncList->Select(newid-1);
      else if( fDim == 1 ) {
         // Select predefined 1D functions for 1D objects
         fTypeFit->Select(kFP_PRED1D, kTRUE);
      } else if( fDim == 2 ) {
         // Select predefined 2D functions for 2D objects
         fTypeFit->Select(kFP_PRED2D, kTRUE);
      }
   }
   // Case for previously used functions.
   else if ( fTypeFit->GetSelected() == kFP_PREVFIT ) {
      Int_t newid = kFP_ALTFUNC;

      // Look only for those functions used in the selected object
      std::pair<fPrevFitIter, fPrevFitIter> look = fPrevFit.equal_range(fFitObject);
      // Then go over all those functions and add them to the list
      for ( fPrevFitIter it = look.first; it != look.second; ++it ) {
         fFuncList->AddEntry(it->second->GetName(), newid++);
      }

      // If no functions were added.
      if ( newid == kFP_ALTFUNC ) {
         // Remove the entry previous fit from fTypeFit
         fTypeFit->RemoveEntry(kFP_PREVFIT);
         if( fDim == 1 )
            // Select predefined 1D functions for 1D objects
            fTypeFit->Select(kFP_PRED1D, kTRUE);
         else if ( fDim == 2 )
            // Select predefined 2D functions for 2D objects
            fTypeFit->Select(kFP_PRED2D, kTRUE);
         else
            // For more than 2 dimensions, select the user functions.
            fTypeFit->Select(kFP_UFUNC, kTRUE);
      }
      else
         // If there is there are previously used functions, select
         // the last one inserted.
         fFuncList->Select(newid-1, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the list of methods depending on the minimization library
/// selected.

void TFitEditor::FillMinMethodList(Int_t)
{
   fMinMethodList->RemoveAll();

   if ( fLibMinuit->GetState() == kButtonDown )
   {
      fMinMethodList->AddEntry("MIGRAD" ,       kFP_MIGRAD);
      fMinMethodList->AddEntry("SIMPLEX" ,      kFP_SIMPLX);
      fMinMethodList->AddEntry("SCAN" ,         kFP_SCAN);
      fMinMethodList->AddEntry("Combination" ,  kFP_COMBINATION);
      fMinMethodList->Select(kFP_MIGRAD, kFALSE);
      fStatusBar->SetText("MIGRAD",2);
   } else if ( fLibFumili->GetState() == kButtonDown )
   {
      fMinMethodList->AddEntry("FUMILI" , kFP_FUMILI);
      fMinMethodList->Select(kFP_FUMILI, kFALSE);
      fStatusBar->SetText("FUMILI",2);
   } else if ( fLibGSL->GetState() == kButtonDown )
   {
      fMinMethodList->AddEntry("Fletcher-Reeves conjugate gradient" , kFP_GSLFR);
      fMinMethodList->AddEntry("Polak-Ribiere conjugate gradient" ,   kFP_GSLPR);
      fMinMethodList->AddEntry("BFGS conjugate gradient" ,            kFP_BFGS);
      fMinMethodList->AddEntry("BFGS conjugate gradient (Version 2)", kFP_BFGS2);
      fMinMethodList->AddEntry("Levenberg-Marquardt" ,                kFP_GSLLM);
      fMinMethodList->AddEntry("Simulated Annealing" ,                kFP_GSLSA);
      fMinMethodList->Select(kFP_GSLFR, kFALSE);
      fStatusBar->SetText("CONJFR",2);
   } else if ( fLibGenetics->GetState() == kButtonDown )
   {
      if ( gPluginMgr->FindHandler("ROOT::Math::Minimizer","GAlibMin") ) {
         fMinMethodList->AddEntry("GA Lib Genetic Algorithm" , kFP_GALIB);
         fMinMethodList->Select(kFP_GALIB, kFALSE);
      } else if (gPluginMgr->FindHandler("ROOT::Math::Minimizer","Genetic")) {
         fMinMethodList->AddEntry("TMVA Genetic Algorithm" ,   kFP_TMVAGA);
         fMinMethodList->Select(kFP_TMVAGA, kFALSE);
      }
   } else // if ( fLibMinuit2->GetState() == kButtonDown )
   {
      fMinMethodList->AddEntry("MIGRAD" ,       kFP_MIGRAD);
      fMinMethodList->AddEntry("SIMPLEX" ,      kFP_SIMPLX);
      fMinMethodList->AddEntry("FUMILI" ,       kFP_FUMILI);
      fMinMethodList->AddEntry("SCAN" ,         kFP_SCAN);
      fMinMethodList->AddEntry("Combination" ,  kFP_COMBINATION);
      fMinMethodList->Select(kFP_MIGRAD, kFALSE);
      fStatusBar->SetText("MIGRAD",2);
   }
}

void SearchCanvases(TSeqCollection* canvases, std::vector<TObject*>& objects)
{
   // Auxiliary function to recursively search for objects inside the
   // current canvases.

   TIter canvasIter(canvases);
   // Iterate over all the canvases in canvases.
   while(TObject* obj = (TObject*) canvasIter()) {
      // If the object is another canvas, call this function
      // recursively.
      if ( TPad* can = dynamic_cast<TPad*>(obj))
         SearchCanvases(can->GetListOfPrimitives(), objects);
      // Otherwhise, if it's a recognised object, add it to the vector
      else if (    dynamic_cast<TH1*>(obj)
                || dynamic_cast<TGraph*>(obj)
                || dynamic_cast<TGraph2D*>(obj)
                || dynamic_cast<TMultiGraph*>(obj)
                || dynamic_cast<THStack*>(obj)
                || dynamic_cast<TTree*>(obj) ) {
         bool insertNew = true;
         // Be careful no to insert the same element twice.
         for ( std::vector<TObject*>::iterator i = objects.begin(); i != objects.end(); ++i )
            if ( (*i) == obj ) {
               insertNew = false;
               break;
            }
         // If the object is not already in the vector, then insert
         // it.
         if ( insertNew ) objects.push_back(obj);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a combo box with all the possible objects to be fitted.

void TFitEditor::FillDataSetList()
{
   // Get the title of the entry selected, so that we can select it
   // again once the fDataSet has been refilled.
   TGTextLBEntry * entry = (TGTextLBEntry*) fDataSet->GetSelectedEntry();
   TString selEntryStr;
   if ( entry ) {
      selEntryStr = entry->GetTitle();
   }

   // Remove all the elements
   fDataSet->RemoveAll();
   std::vector<TObject*> objects;

   // Get all the objects registered in gDirectory
   if (gDirectory) {
      TList * l = gDirectory->GetList();
      if (l) {
         TIter next(l);
         TObject* obj = NULL;
         while ( (obj = (TObject*) next()) ) {
            // But only if they are of a type recognized by the FitPanel
            if ( dynamic_cast<TH1*>(obj) ||
                 dynamic_cast<TGraph2D*>(obj) ||
                 dynamic_cast<TTree*>(obj) ) {
               objects.push_back(obj);
            }
         }
      }
   }

   // Look for all the drawn objects. The method will take care the
   // same objects are not inserted twice.
   SearchCanvases(gROOT->GetListOfCanvases(), objects);

   // Add all the objects stored in the vector
   int selected = kFP_NOSEL;
   // Add the No selection.
   Int_t newid = kFP_NOSEL;
   fDataSet->AddEntry("No Selection", newid++);
   for ( std::vector<TObject*>::iterator i = objects.begin(); i != objects.end(); ++i ) {
      // Insert the name as the class name followed by the name of the
      // object.
      TString name = (*i)->ClassName(); name.Append("::"); name.Append((*i)->GetName());
      // Check whether the names are the same!
      if ( selEntryStr && name == selEntryStr )
         selected = newid;
      fDataSet->AddEntry(name, newid++);
   }

   // If there was an entry selected (which should be always the case
   // except the first time this method is executed), then make it the
   // selected one again.
   if (entry) {
      fDataSet->Select(selected);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create method list in a combo box.

TGComboBox* TFitEditor::BuildMethodList(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);
   c->AddEntry("Chi-square", kFP_MCHIS);
   c->AddEntry("Binned Likelihood", kFP_MBINL);
   c->AddEntry("Unbinned Likelihood", kFP_MUBIN);
   //c->AddEntry("User", kFP_MUSER);                //for later use
   c->Select(kFP_MCHIS);
   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to advanced option button (opens a dialog).

void TFitEditor::DoAdvancedOptions()
{
   new TAdvancedGraphicsDialog( fClient->GetRoot(), GetMainFrame());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to 'include emtry bins and forse all weights to 1' setting.

void TFitEditor::DoEmptyBinsAllWeights1()
{
   if (fEmptyBinsWghts1->GetState() == kButtonDown)
      if (fAllWeights1->GetState() == kButtonDown)
         fAllWeights1->SetState(kButtonUp, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////

void TFitEditor::DoUseFuncRange()
{
   if ( fUseRange->GetState() == kButtonDown ) {
      if (fNone->GetState() == kButtonDown || fNone->GetState() == kButtonDisabled) {
         // Get the function
         TF1* tmpTF1 = FindFunction();
         if ( !tmpTF1 ) {
            if (GetFitObjectListOfFunctions()) {
               TGTextLBEntry *te = (TGTextLBEntry *)fFuncList->GetSelectedEntry();
               tmpTF1 = (TF1*) GetFitObjectListOfFunctions()->FindObject( te->GetTitle() );
            }
         }
         // If the function has been retrieved, i.e. is a registered function.
         if ( tmpTF1 ) {
            Double_t xmin, ymin, zmin, xmax, ymax, zmax;
            // Get the range
            tmpTF1->GetRange(xmin, ymin, zmin, xmax, ymax, zmax);
            // And set the sliders
            if ( fType != kObjectTree ) {
               fSliderXMin->SetNumber( xmin );
               fSliderXMax->SetNumber( xmax );
               DoNumericSliderXChanged();
               if ( fDim > 1 ) {
                  fSliderYMin->SetNumber( ymin );
                  fSliderYMax->SetNumber( ymax );
                  DoNumericSliderYChanged();
               }
            }
         }
      }
      fUseRange->SetState(kButtonDown);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to 'set all weights to 1' setting.

void TFitEditor::DoAllWeights1()
{
   if (fAllWeights1->GetState() == kButtonDown)
      if (fEmptyBinsWghts1->GetState() == kButtonDown)
         fEmptyBinsWghts1->SetState(kButtonUp, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Close the fit panel.

void TFitEditor::DoClose()
{
   Hide();
}

////////////////////////////////////////////////////////////////////////////////
/// Easy here!

void TFitEditor::DoUpdate()
{
   GetFunctionsFromSystem();
   FillDataSetList();
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a fit with current parameters' settings.

void TFitEditor::DoFit()
{
   if (!fFitObject) return;
   //if (!fParentPad) return;

   // If fNone->GetState() == kButtonDisabled means the function is
   // not editable, i.e. it comes from a raw C function. So in this
   // case, it is editable and we have to check wheather the formula
   // is well built.
   if ( fNone->GetState() != kButtonDisabled && CheckFunctionString(fEnteredFunc->GetText()) )
   {
      // If not, then show an error message and leave.
      new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                   "Error...", "2) Verify the entered function string!",
                   kMBIconStop,kMBOk, 0);
      return;
   }

   // Set the button so that the user cannot use it while fitting, set
   // the mouse to watch type and so on.
   fFitButton->SetState(kButtonEngaged);
   if (gPad && gPad->GetVirtCanvas()) gPad->GetVirtCanvas()->SetCursor(kWatch);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kWatch));

   TVirtualPad *save = nullptr;
   if ( fParentPad ) {
      fParentPad->Disconnect("RangeAxisChanged()");
      save = gPad;
      gPad = fParentPad;
      fParentPad->cd();

      if (fParentPad->GetCanvas())
         fParentPad->GetCanvas()->SetCursor(kWatch);
   }

   // Get the ranges from the sliders
   ROOT::Fit::DataRange drange;
   GetRanges(drange);

   // Create a static pointer to fitFunc. Every second call to the
   // DoFit method, the old fitFunc is deleted. We need not to delete
   // the function after the fitting in case we want to do Advaced
   // graphics. The VirtualFitter need the function to be alived. One
   // problem, after the last fit the function is never deleted, but
   // ROOT's garbage collector will do the job for us.
   static TF1 *fitFunc = nullptr;
   if ( fitFunc ) {
      //std::cout << "TFitEditor::DoFit - deleting fit function " << fitFunc->GetName() << "  " << fitFunc << std::endl;
      delete fitFunc;
   }
   fitFunc = GetFitFunction();

   std::cout << "TFitEditor::DoFit - using function " << fitFunc->GetName() << "  " << fitFunc << std::endl;
   // This assert
   if (!fitFunc) {
      Error("DoFit","This should have never happend, the fitfunc pointer is NULL! - Please Report" );
      return;
   }

   // set parameters from panel in function
   SetParameters(fFuncPars, fitFunc);
   // Get the options stored in the GUI elements.
   ROOT::Math::MinimizerOptions mopts;
   Foption_t fitOpts;
   TString strDrawOpts;
   RetrieveOptions(fitOpts, strDrawOpts, mopts, fitFunc->GetNpar());

   // Call the fit method, depending on the object to fit.
   switch (fType) {
      case kObjectHisto: {

         TH1 *hist = dynamic_cast<TH1*>(fFitObject);
         if (hist)
            ROOT::Fit::FitObject(hist, fitFunc, fitOpts, mopts, strDrawOpts, drange);

         break;
      }
      case kObjectGraph: {

         TGraph *gr = dynamic_cast<TGraph*>(fFitObject);
         if (gr)
            FitObject(gr, fitFunc, fitOpts, mopts, strDrawOpts, drange);
         break;
      }
      case kObjectMultiGraph: {

         TMultiGraph *mg = dynamic_cast<TMultiGraph*>(fFitObject);
         if (mg)
            FitObject(mg, fitFunc, fitOpts, mopts, strDrawOpts, drange);

         break;
      }
      case kObjectGraph2D: {

         TGraph2D *g2d = dynamic_cast<TGraph2D*>(fFitObject);
         if (g2d)
            FitObject(g2d, fitFunc, fitOpts, mopts, strDrawOpts, drange);

         break;
      }
      case kObjectHStack: {
         // N/A
         break;
      }
      case kObjectTree:  {
         // The three is a much more special case. The steps for
         // fitting have to be done manually here until they are
         // properly implemented within a FitObject method in
         // THFitImpl.cxx

         // Retrieve the variables and cuts selected from the current
         // tree.
         TString variables;
         TString cuts;
         GetTreeVarsAndCuts(fDataSet, variables, cuts);

         // This should be straight forward and the return should
         // never be called.
         TTree *tree = dynamic_cast<TTree*>(fFitObject);
         if ( !tree ) return;

         // These method calls are just to set up everything for the
         // fitting. It's taken from another script.
         gROOT->ls();
         tree->Draw(variables,cuts,"goff");

         auto player = tree->GetPlayer();
         if ( !player ) {
            Error("DoFit","Player reference is NULL");
            return;
         }

         auto selector = dynamic_cast<TSelectorDraw *>(player->GetSelector());
         if ( !selector ) {
            Error("DoFit","Selector reference is NULL");
            return;
         }

         // use pointer stored in the tree (not copy the data in)
         unsigned int ndim = player->GetDimension();
         if ( ndim == 0 ) {
            Error("DoFit","NDIM == 0");
            return;
         }

         std::vector<double *> vlist;
         for (unsigned int i = 0; i < ndim; ++i) {
            double * v =  selector->GetVal(i);
            if (v != 0) vlist.push_back(v);
            else
               std::cerr << "pointer for variable " << i << " is zero" << std::endl;
         }
         if (vlist.size() != ndim) {
            Error("DoFit","Vector is not complete");
            return;
         }

         // fill the data
         Long64_t nrows = player->GetSelectedRows();
         if ( !nrows ) {
            Error("DoFit","NROWS == 0");
            return;
         }

         ROOT::Fit::UnBinData * fitdata = new ROOT::Fit::UnBinData(nrows, ndim, vlist.begin());

         for ( int i = 0; i < std::min(int(fitdata->Size()),10); ++i) {
            // print j coordinate
            for (unsigned int j = 0; j < ndim; ++j) {
               printf(" x_%d [%d] = %f  \n", j, i,*(fitdata->Coords(i)+j) );
            }
            printf("\n");
         }


         //TVirtualFitter::SetDefaultFitter("Minuit");
         Foption_t fitOption;
         ROOT::Math::MinimizerOptions minOption;
         fitOption.Verbose=1;

         // After all the set up is performed, then do the Fit!!
         ROOT::Fit::UnBinFit(fitdata, fitFunc, fitOption, minOption);

         break;
      }
   }

   // if SAME is set re-plot the function
   // useful in case histogram was drawn with HIST
   //  and no function will be drawm)
   if (fDrawSame->GetState() == kButtonDown && fitFunc)
      fitFunc->Draw("same");


   // update parameters value shown in dialog
   //if (!fFuncPars) fFuncPars = new Double_t[fitFunc->GetNpar()][3];
   GetParameters(fFuncPars,fitFunc);

   // Save fit data for future use as a PrevFit function.
   TF1* tmpTF1 = copyTF1(fitFunc);
   TString name = TString::Format("PrevFit-%d", (int) fPrevFit.size() + 1);
   if (!strstr(fitFunc->GetName(),"PrevFit"))
      name.Append(TString::Format("-%s", fitFunc->GetName()));
   tmpTF1->SetName(name.Data());
   fPrevFit.emplace(fFitObject, tmpTF1);
   fSystemFuncs.emplace_back( copyTF1(tmpTF1) );

   float xmin = 0.f, xmax = 0.f, ymin = 0.f, ymax = 0.f, zmin = 0.f, zmax = 0.f;
   if ( fParentPad ) {
      fParentPad->Modified();
      // As the range is not changed, save the old values and restore
      // after the GUI has been updated.  It would be more elegant to
      // disconnect the signal from fParentPad, however, this doesn't
      // work for unknown reasons.
      if ( fType != kObjectTree ) fSliderX->GetPosition(xmin, xmax);
      if ( fDim > 1 ) fSliderY->GetPosition(ymin, ymax);
      if ( fDim > 2 ) fSliderZ->GetPosition(zmin, zmax);
      fParentPad->Update();
   }

   // In case the fit method draws something! Set the canvas!
   fParentPad = gPad;
   UpdateGUI();

   // Change the sliders if necessary.
   if ( fParentPad ) {
      if ( fType != kObjectTree ) { fSliderX->SetPosition(xmin, xmax); DoSliderXMoved(); }
      if ( fType != kObjectTree && fDim > 1 ) { fSliderY->SetPosition(ymin, ymax); DoSliderYMoved(); }
      if ( fType != kObjectTree && fDim > 2 ) { fSliderZ->SetPosition(zmin, zmax); DoSliderZMoved(); }
      if (fParentPad->GetCanvas())
         fParentPad->GetCanvas()->SetCursor(kPointer);
      fParentPad->Connect("RangeAxisChanged()", "TFitEditor", this, "UpdateGUI()");

      if (save) gPad = save;
      if (fSetParam->GetState() == kButtonDisabled &&
          fLinearFit->GetState() == kButtonUp)
         fSetParam->SetState(kButtonUp);
   }

   // Restore the Fit button and mouse cursor to their proper state.
   if (gPad && gPad->GetVirtCanvas()) gPad->GetVirtCanvas()->SetCursor(kPointer);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kPointer));
   fFitButton->SetState(kButtonUp);

   if ( !fTypeFit->FindEntry("Prev. Fit") )
      fTypeFit->InsertEntry("Prev. Fit",kFP_PREVFIT, kFP_UFUNC);

   fDrawAdvanced->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Check entered function string.

Int_t TFitEditor::CheckFunctionString(const char *fname)
{
   Int_t rvalue = 0;
   if ( fDim == 1 || fDim == 0 ) {
      TF1 form("tmpCheck", fname);
      // coverity[uninit_use_in_call]
      rvalue = form.IsValid() ? 0 : -1;
   } else if ( fDim == 2 ) {
      TF2 form("tmpCheck", fname);
      // coverity[uninit_use_in_call]
      rvalue = form.IsValid() ? 0 : -1;
   } else if ( fDim == 3 ) {
      TF3 form("tmpCheck", fname);
      // coverity[uninit_use_in_call]
      rvalue = form.IsValid() ? 0 : -1;
   }

   return rvalue;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to addition of predefined functions. It will
/// insert the next selected function with a plus sign so that it
/// doesn't override the current content of the formula.

void TFitEditor::DoAddition(Bool_t on)
{
   static Bool_t first = kFALSE;
   TString s = fEnteredFunc->GetText();
   if (on) {
      if (!first) {
         fSelLabel->SetText(s.Sizeof()>30?s(0,30)+"...":s);
         s += "(0)";
         fEnteredFunc->SetText(s.Data());
         first = kTRUE;
         ((TGCompositeFrame *)fSelLabel->GetParent())->Layout();
      }
   } else {
      first = kFALSE;
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Slot connected to addition of predefined functions. It will
/// insert the next selected function with a plus sign so that it
/// doesn't override the current content of the formula.

void TFitEditor::DoNormAddition(Bool_t on)
{
   /*
   static Bool_t first = kFALSE;
   TString s = fEnteredFunc->GetText();
   if (on) {
      if (!first) {
         fSelLabel->SetText(s.Sizeof()>30?s(0,30)+"...":s);
         fEnteredFunc->SetText(s.Data());
         first = kTRUE;
         ((TGCompositeFrame *)fSelLabel->GetParent())->Layout();
      }
   } else {
      first = kFALSE;
   }*/

   if (on) Info("DoNormAddition","Normalized addition is selected");
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to addition of predefined functions. It will
/// insert the next selected function with a plus sign so that it
/// doesn't override the current content of the formula.

void TFitEditor::DoConvolution(Bool_t on)
{
   /*
   static Bool_t first = kFALSE;
   TString s = fEnteredFunc->GetText();
   if (on) {
      if (!first) {
         fSelLabel->SetText(s.Sizeof()>30?s(0,30)+"...":s);
        // s += "(0)";
         fEnteredFunc->SetText(s.Data());
         first = kTRUE;
         ((TGCompositeFrame *)fSelLabel->GetParent())->Layout();
      }
   } else
      first = kFALSE;*/

   if (on) Info("DoConvolution","Convolution is selected");
}

////////////////////////////////////////////////////////////////////////////////
/// Selects the data set to be fitted

void TFitEditor::DoDataSet(Int_t selected)
{
   if ( selected == kFP_NOSEL ) {
      DoNoSelection();
      return;
   }

   // Get the name and class of the selected object.
   TGTextLBEntry* textEntry = static_cast<TGTextLBEntry*>(fDataSet->GetListBox()->GetEntry(selected));
   if (!textEntry) return;
   TString textEntryStr = textEntry->GetText()->GetString();
   TString name = textEntry->GetText()->GetString()+textEntry->GetText()->First(':')+2;
   TString className = textEntryStr(0,textEntry->GetText()->First(':'));

   // Check the object exists in the ROOT session and it is registered
   TObject* objSelected(0);
   if ( className == "TTree" ) {
      // It's a tree, so the name is before the space (' ')
      TString lookStr;
      if ( name.First(' ') == kNPOS )
         lookStr = name;
      else
         lookStr = name(0, name.First(' '));
      //std::cout << "\t1 SITREE: '" << lookStr << "'" << std::endl;
      objSelected = gROOT->FindObject(lookStr);
   } else {
      // It's not a tree, so the name is the complete string
      //std::cout << "\t1 NOTREE: '" << name << "'" << std::endl;
      objSelected = gROOT->FindObject(name);
   }
   if ( !objSelected )
   {
      //std::cerr << "Object not found! Please report the error! " << std::endl;
      return;
   }

   // If it is a tree, and there are no variables selected, show a dialog
   if ( objSelected->InheritsFrom(TTree::Class()) &&
        name.First(' ') == kNPOS ) {
      char variables[256] = {0}; char cuts[256] = {0};
      strlcpy(variables, "Sin input!", 256);
      new TTreeInput( fClient->GetRoot(), GetMainFrame(), variables, cuts );
      if ( strcmp ( variables, "" ) == 0 ) {
         DoNoSelection();
         return;
      }
      ProcessTreeInput(objSelected, selected, variables, cuts);
   }

   // Search the canvas where the object is drawn, if any
   TPad* currentPad = NULL;
   bool found = false;
   std::queue<TPad*> stPad;
   TIter padIter( gROOT->GetListOfCanvases() );
   while ( TObject* canvas = static_cast<TObject*>(padIter() ) ) {
      if ( dynamic_cast<TPad*>(canvas) )
         stPad.push(dynamic_cast<TPad*>(canvas));
   }

   while ( !stPad.empty() && !found ) {
      currentPad = stPad.front();
      stPad.pop();
      TIter elemIter( currentPad->GetListOfPrimitives() );
      while ( TObject* elem = static_cast<TObject*>(elemIter() ) ) {
         if ( elem == objSelected ) {
            found = true;
            break;
         } else if ( dynamic_cast<TPad*>(elem) )
            stPad.push( dynamic_cast<TPad*>(elem) );
      }
   }

   // Set the proper object and canvas (if found!)
   SetFitObject( found ? currentPad : nullptr, objSelected, kButton1Down);
}

void TFitEditor::ProcessTreeInput(TObject* objSelected, Int_t selected, TString variables, TString cuts)
{
   // If the input is valid, insert the tree with the selections as an entry to fDataSet
   TString entryName = (objSelected)->ClassName(); entryName.Append("::"); entryName.Append((objSelected)->GetName());
   entryName.Append(" (\""); entryName.Append(variables); entryName.Append("\", \"");
   entryName.Append(cuts); entryName.Append("\")");
   Int_t newid = fDataSet->GetNumberOfEntries() + kFP_NOSEL;
   fDataSet->InsertEntry(entryName, newid, selected );
   fDataSet->Select(newid);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to predefined fit function settings.

void TFitEditor::DoFunction(Int_t selected)
{
   TGTextLBEntry *te = (TGTextLBEntry *)fFuncList->GetSelectedEntry();

   // check that selected passesd value is the correct one in the TextEntry
   R__ASSERT( selected == te->EntryId());
   //std::cout << "calling do function " << selected << "  " << te->GetTitle() << " function " << te->EntryId() << std::endl;
   //selected = te->EntryId();

   bool editable = false;
   if (fNone -> GetState() == kButtonDown || fNone->GetState() == kButtonDisabled)
   {
      // Get the function selected and check weather it is a raw C
      // function or not
      TF1* tmpTF1 = FindFunction();
      if ( !tmpTF1 )
      {
         if (GetFitObjectListOfFunctions())
            tmpTF1 = (TF1*) GetFitObjectListOfFunctions()->FindObject( te->GetTitle() );
      }
      if ( tmpTF1 && strcmp(tmpTF1->GetExpFormula(), "") )
      {
         editable = kTRUE;
         fEnteredFunc->SetText(tmpTF1->GetExpFormula());
      }
      else
      {
         if ( selected <= kFP_USER )
            editable = kTRUE;
         else
            editable = kFALSE;
         fEnteredFunc->SetText(te->GetTitle());
      }
      // Once you have the function, set the editable.
      SetEditable(editable);
   }
   else if (fAdd  -> GetState() == kButtonDown)
   {
      // If the add button is down don't replace the fEnteredFunc text
      Int_t np = 0;
      TString s = "";
      if (!strcmp(fEnteredFunc->GetText(), ""))
      {
         fEnteredFunc->SetText(te->GetTitle());
      }
      else
      {
         s = fEnteredFunc->GetTitle();
         TFormula tmp("tmp", fEnteredFunc->GetText());
         np = tmp.GetNpar();
      }
      if (np)
         s += TString::Format("+%s(%d)", te->GetTitle(), np);
      else
         s += TString::Format("%s(%d)", te->GetTitle(), np);
      fEnteredFunc->SetText(s.Data());
      editable = true;
   }
   else if (fNormAdd->GetState() == kButtonDown)
   {
      // If the normadd button is down don't replace the fEnteredFunc text
      Int_t np = 0;
      TString s = "";
      if (!strcmp(fEnteredFunc->GetText(), ""))
      {
         fEnteredFunc->SetText(te->GetTitle());
      }
      else
      {
         s = fEnteredFunc->GetTitle();
         TFormula tmp("tmp", fEnteredFunc->GetText());
         np = tmp.GetNpar();
      }
      if (np)
         s += TString::Format("+%s", te->GetTitle());
      else
         s += TString::Format("%s", te->GetTitle());
      fEnteredFunc->SetText(s.Data());
      //std::cout <<fEnteredFunc->GetText()<<std::endl;
      editable = true;
   }
   else if (fConv->GetState() == kButtonDown)
   {
      // If the normadd button is down don't replace the fEnteredFunc text
      Int_t np = 0;
      TString s = "";
      if (!strcmp(fEnteredFunc->GetText(), ""))
         fEnteredFunc->SetText(te->GetTitle());
      else
      {
         s = fEnteredFunc->GetTitle();
         TFormula tmp("tmp", fEnteredFunc->GetText());
         np = tmp.GetNpar();
      }
      if (np)
         s += TString::Format("*%s", te->GetTitle());
      else
         s += TString::Format("%s", te->GetTitle());
      fEnteredFunc->SetText(s.Data());
      //std::cout <<fEnteredFunc->GetText()<<std::endl;
      editable = true;
   }


   // Get the final name in fEnteredFunc to process the function that
   // it would create
   TString tmpStr = fEnteredFunc->GetText();

   // create TF1 with the passed string. Delete previous one if existing
   if (tmpStr.Contains("pol") || tmpStr.Contains("++")) {
      fLinearFit->SetState(kButtonDown, kTRUE);
   } else {
      fLinearFit->SetState(kButtonUp, kTRUE);
   }

   fEnteredFunc->SelectAll();
   fSelLabel->SetText(tmpStr.Sizeof()>30?tmpStr(0,30)+"...":tmpStr);
   ((TGCompositeFrame *)fSelLabel->GetParent())->Layout();

   // reset function parameters if the number of parameters of the new
   // function is different from the old one!
   TF1* fitFunc = GetFitFunction();
   //std::cout << "TFitEditor::DoFunction - using function " << fitFunc->GetName() << "  " << fitFunc << std::endl;

   if ( fitFunc && (unsigned int) fitFunc->GetNpar() != fFuncPars.size() )
      fFuncPars.clear();
   if ( fitFunc ) {
      //std::cout << "TFitEditor::DoFunction - deleting function " << fitFunc->GetName() << "  " << fitFunc << std::endl;
      delete fitFunc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to entered function in text entry.

void TFitEditor::DoEnteredFunction()
{
   if (!strcmp(fEnteredFunc->GetText(), "")) return;

   // Check if the function is well built
   Int_t ok = CheckFunctionString(fEnteredFunc->GetText());

   if (ok != 0) {
      new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                   "Error...", "3) Verify the entered function string!",
                   kMBIconStop,kMBOk, 0);
      return;
   }

   // And set the label with the entered text if everything is fine.
   TString s = fEnteredFunc->GetText();
   fSelLabel->SetText(s.Sizeof()>30?s(0,30)+"...":s);
   ((TGCompositeFrame *)fSelLabel->GetParent())->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to linear fit settings.

void TFitEditor::DoLinearFit()
{
   if (fLinearFit->GetState() == kButtonDown) {
      //fSetParam->SetState(kButtonDisabled);
      fBestErrors->SetState(kButtonDisabled);
      fImproveResults->SetState(kButtonDisabled);
      fEnableRobust->SetState(kButtonUp);
      //fNoChi2->SetState(kButtonUp);
   } else {
      //fSetParam->SetState(kButtonUp);
      fBestErrors->SetState(kButtonUp);
      fImproveResults->SetState(kButtonUp);
      fEnableRobust->SetState(kButtonDisabled);
      fRobustValue->SetState(kFALSE);
      //fNoChi2->SetState(kButtonDisabled);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to 'no chi2' option settings.

void TFitEditor::DoNoChi2()
{
   //LM: no need to do  operations here
   // if (fLinearFit->GetState() == kButtonUp)
   //    fLinearFit->SetState(kButtonDown, kTRUE);
}
////////////////////////////////////////////////////////////////////////////////
/// Slot connected to 'robust fitting' option settings.

void TFitEditor::DoRobustFit()
{
   if (fEnableRobust->GetState() == kButtonDown)
      fRobustValue->SetState(kTRUE);
   else
      fRobustValue->SetState(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to 'no storing, no drawing' settings.

void TFitEditor::DoNoStoreDrawing()
{
   if (fNoDrawing->GetState() == kButtonUp)
      fNoDrawing->SetState(kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to print option settings.

void TFitEditor::DoPrintOpt(Bool_t on)
{
   // Change the states of the buttons depending of which one is
   // selected.
   TGButton *btn = (TGButton *) gTQSender;
   Int_t id = btn->WidgetId();
   switch (id) {
      case kFP_PDEF:
         if (on) {
            fOptDefault->SetState(kButtonDown);
            fOptVerbose->SetState(kButtonUp);
            fOptQuiet->SetState(kButtonUp);
         }
         fStatusBar->SetText("Prn: DEF",4);
         break;
      case kFP_PVER:
         if (on) {
            fOptVerbose->SetState(kButtonDown);
            fOptDefault->SetState(kButtonUp);
            fOptQuiet->SetState(kButtonUp);
         }
         fStatusBar->SetText("Prn: VER",4);
         break;
      case kFP_PQET:
         if (on) {
            fOptQuiet->SetState(kButtonDown);
            fOptDefault->SetState(kButtonUp);
            fOptVerbose->SetState(kButtonUp);
         }
         fStatusBar->SetText("Prn: QT",4);
      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset all fit parameters.

void TFitEditor::DoReset()
{
   if ( fParentPad ) {
      fParentPad->Modified();
      fParentPad->Update();
   }
   fEnteredFunc->SetText("gaus");

   // To restore temporary points and sliders
   UpdateGUI();

   if (fLinearFit->GetState() == kButtonDown)
      fLinearFit->SetState(kButtonUp, kTRUE);
   if (fBestErrors->GetState() == kButtonDown)
      fBestErrors->SetState(kButtonUp, kFALSE);
   if (fUseRange->GetState() == kButtonDown)
      fUseRange->SetState(kButtonUp, kFALSE);
   if (fAllWeights1->GetState() == kButtonDown)
      fAllWeights1->SetState(kButtonUp, kFALSE);
   if (fEmptyBinsWghts1->GetState() == kButtonDown)
      fEmptyBinsWghts1->SetState(kButtonUp, kFALSE);
   if (fImproveResults->GetState() == kButtonDown)
      fImproveResults->SetState(kButtonUp, kFALSE);
   if (fAdd2FuncList->GetState() == kButtonDown)
      fAdd2FuncList->SetState(kButtonUp, kFALSE);
   if (fUseGradient->GetState() == kButtonDown)
      fUseGradient->SetState(kButtonUp, kFALSE);
   if (fEnableRobust->GetState() == kButtonDown)
      fEnableRobust->SetState(kButtonUp, kFALSE);
   // if (fNoChi2->GetState() == kButtonDown)
   //    fNoChi2->SetState(kButtonUp, kFALSE);
   if (fDrawSame->GetState() == kButtonDown)
      fDrawSame->SetState(kButtonUp, kFALSE);
   if (fNoDrawing->GetState() == kButtonDown)
      fNoDrawing->SetState(kButtonUp, kFALSE);
   if (fNoStoreDrawing->GetState() == kButtonDown)
      fNoStoreDrawing->SetState(kButtonUp, kFALSE);
   fNone->SetState(kButtonDown, kTRUE);
   fFuncList->Select(1, kTRUE);

   // minimization tab
   if (fLibMinuit->GetState() != kButtonDown)
      fLibMinuit->SetState(kButtonDown, kTRUE);
   FillMinMethodList();
   if (fOptDefault->GetState() != kButtonDown)
      fOptDefault->SetState(kButtonDown, kTRUE);
   if (fErrorScale->GetNumber() != ROOT::Math::MinimizerOptions::DefaultErrorDef()) {
      fErrorScale->SetNumber(ROOT::Math::MinimizerOptions::DefaultErrorDef());
      fErrorScale->ReturnPressed();
   }
   if (fTolerance->GetNumber() != ROOT::Math::MinimizerOptions::DefaultTolerance()) {
      fTolerance->SetNumber(ROOT::Math::MinimizerOptions::DefaultTolerance());
      fTolerance->ReturnPressed();
   }
   if (fIterations->GetNumber() != ROOT::Math::MinimizerOptions::DefaultMaxIterations()) {
      fIterations->SetIntNumber(ROOT::Math::MinimizerOptions::DefaultMaxIterations());
      fIterations->ReturnPressed();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Open set parameters dialog.

void TFitEditor::DoSetParameters()
{
   // Get the function.
   TF1* fitFunc = GetFitFunction();
   //std::cout << "TFitEditor::DoSetParameters - using function " << fitFunc->GetName() << "  " << fitFunc << std::endl;

   if (!fitFunc) { Error("DoSetParameters","NUll function"); return; }

   // case of special functions (gaus, expo, etc...) if the function
   // has not defined the parameters yet. For those, don't let the
   // parameters to be all equal to 0, as we can provide some good
   // starting value.
   if (fFuncPars.size() == 0) {
      switch (fType) {
      case kObjectHisto:
         InitParameters( fitFunc, (TH1*)fFitObject) ;
         break;
      case kObjectGraph:
         InitParameters( fitFunc, ((TGraph*)fFitObject));
         break;
      case kObjectMultiGraph:
         InitParameters( fitFunc, ((TMultiGraph*)fFitObject));
         break;
      case kObjectGraph2D:
         InitParameters( fitFunc, ((TGraph2D*)fFitObject));
         break;
      case kObjectHStack:
      case kObjectTree:
      default:
         break;
      }
      // The put these parameters into the fFuncPars structure
      GetParameters(fFuncPars, fitFunc);
   }
   else {
      // Otherwise, put the parameters in the function
      SetParameters(fFuncPars, fitFunc);
   }

   if ( fParentPad ) fParentPad->Disconnect("RangeAxisChanged()");
   Int_t ret = 0;
   /// fit parameter dialog willbe deleted automatically when closed
   new TFitParametersDialog(gClient->GetDefaultRoot(), GetMainFrame(),
                            fitFunc, fParentPad, &ret);

   // Once the parameters are set in the fitfunction, save them.
   GetParameters(fFuncPars, fitFunc);

   // check return code to see if parameters settings have been modified
   // in this case we need to set the B option when fitting
   if (ret) fChangedParams = kTRUE;


   if ( fParentPad ) fParentPad->Connect("RangeAxisChanged()", "TFitEditor", this, "UpdateGUI()");

   if ( fNone->GetState() != kButtonDisabled ) {
      //std::cout << "TFitEditor::DoSetParameters - deleting function " << fitFunc->GetName() << "  " << fitFunc << std::endl;
      delete fitFunc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to range settings on x-axis.

void TFitEditor::DoSliderXMoved()
{
   if ( !fFitObject ) return;

   fSliderXMin->SetNumber( fXaxis->GetBinLowEdge( static_cast<Int_t>( fSliderX->GetMinPosition() ) ) );
   fSliderXMax->SetNumber( fXaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderX->GetMaxPosition() ) ) );

   fUseRange->SetState(kButtonUp);

   DrawSelection();
}

////////////////////////////////////////////////////////////////////////////////
/// Draws the square around the object showing where the limits for
/// fitting are.

void TFitEditor::DrawSelection(bool restore)
{
   static Int_t  px1old, py1old, px2old, py2old; // to remember the square drawn.

   if ( !fParentPad ) return;

   if (restore) {
      px1old = fParentPad->XtoAbsPixel(fParentPad->GetUxmin());
      py1old = fParentPad->YtoAbsPixel(fParentPad->GetUymin());
      px2old = fParentPad->XtoAbsPixel(fParentPad->GetUxmax());
      py2old = fParentPad->YtoAbsPixel(fParentPad->GetUymax());
      return;
   }

   Int_t px1,py1,px2,py2;

   TVirtualPad *save = 0;
   save = gPad;
   gPad = fParentPad;
   gPad->cd();

   Double_t xleft = 0;
   Double_t xright = 0;
   xleft  = fXaxis->GetBinLowEdge((Int_t)((fSliderX->GetMinPosition())+0.5));
   xright = fXaxis->GetBinUpEdge((Int_t)((fSliderX->GetMaxPosition())+0.5));

   Float_t ymin, ymax;
   if ( fDim > 1 )
   {
      ymin = fYaxis->GetBinLowEdge((Int_t)((fSliderY->GetMinPosition())+0.5));//gPad->GetUymin();
      ymax = fYaxis->GetBinUpEdge((Int_t)((fSliderY->GetMaxPosition())+0.5));//gPad->GetUymax();
   }
   else
   {
      ymin = gPad->GetUymin();
      ymax = gPad->GetUymax();
   }

   px1 = gPad->XtoAbsPixel(xleft);
   py1 = gPad->YtoAbsPixel(ymin);
   px2 = gPad->XtoAbsPixel(xright);
   py2 = gPad->YtoAbsPixel(ymax);

   if (gPad->GetCanvas()) gPad->GetCanvas()->FeedbackMode(kTRUE);
   gPad->SetLineWidth(1);
   gPad->SetLineColor(2);

   gVirtualX->DrawBox(px1old, py1old, px2old, py2old, TVirtualX::kHollow);
   gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);

   px1old = px1;
   py1old = py1;
   px2old = px2 ;
   py2old = py2;

   if(save) gPad = save;
}

////////////////////////////////////////////////////////////////////////////////
/// Sincronize the numeric sliders with the graphical one.

void TFitEditor::DoNumericSliderXChanged()
{
   if ( fSliderXMin->GetNumber() > fSliderXMax->GetNumber() ) {
      float xmin, xmax;
      fSliderX->GetPosition(xmin, xmax);
      fSliderXMin->SetNumber( fXaxis->GetBinLowEdge( static_cast<Int_t>( xmin ) ) );
      fSliderXMax->SetNumber( fXaxis->GetBinUpEdge ( static_cast<Int_t>( xmax ) ) );
      return;
   }

   fSliderX->SetPosition(fXaxis->FindBin( fSliderXMin->GetNumber() ),
                         fXaxis->FindBin( fSliderXMax->GetNumber() ));

   fUseRange->SetState(kButtonUp);

   DrawSelection();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to range settings on y-axis.

void TFitEditor::DoSliderYMoved()
{
   if ( !fFitObject ) return;

   fSliderYMin->SetNumber( fYaxis->GetBinLowEdge( static_cast<Int_t>( fSliderY->GetMinPosition() ) ) );
   fSliderYMax->SetNumber( fYaxis->GetBinUpEdge ( static_cast<Int_t>( fSliderY->GetMaxPosition() ) ) );

   fUseRange->SetState(kButtonUp);

   DrawSelection();
}

////////////////////////////////////////////////////////////////////////////////
///syncronize the numeric slider with the graphical one.

void TFitEditor::DoNumericSliderYChanged()
{
   if ( fSliderYMin->GetNumber() > fSliderYMax->GetNumber() ) {
      float ymin, ymax;
      fSliderY->GetPosition(ymin, ymax);
      fSliderYMin->SetNumber( fYaxis->GetBinLowEdge( static_cast<Int_t>( ymin ) ) );
      fSliderYMax->SetNumber( fYaxis->GetBinUpEdge ( static_cast<Int_t>( ymax ) ) );
      return;
   }

   fSliderY->SetPosition( fYaxis->FindBin( fSliderYMin->GetNumber() ),
                          fYaxis->FindBin( fSliderYMax->GetNumber() ));

   fUseRange->SetState(kButtonUp);

   DrawSelection();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to range settings on z-axis.

void TFitEditor::DoSliderZMoved()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Open a dialog for getting a user defined method.

void TFitEditor::DoUserDialog()
{
   new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                "Info", "Dialog of user method is not implemented yet",
                kMBIconAsterisk,kMBOk, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the function to be used in performed fit.

void TFitEditor::SetFunction(const char *function)
{
   fEnteredFunc->SetText(function);
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether the object suitable for fitting and set
/// its type, dimension and method combo box accordingly.

Bool_t TFitEditor::SetObjectType(TObject* obj)
{
   Bool_t set = kFALSE;

   // For each kind of object, set a different status in the fit
   // panel.
   if (obj->InheritsFrom(TGraph::Class())) {
      fType = kObjectGraph;
      set = kTRUE;
      fDim = 1;
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Chi-square", kFP_MCHIS);
      fMethodList->Select(kFP_MCHIS, kFALSE);
      fRobustValue->SetState(kTRUE);
      fRobustValue->GetNumberEntry()->SetToolTipText("Set robust value");
   } else if (obj->InheritsFrom(TGraph2D::Class())) {
      fType = kObjectGraph2D;
      set = kTRUE;
      fDim = 2;
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Chi-square", kFP_MCHIS);
      fMethodList->Select(kFP_MCHIS, kFALSE);
   } else if (obj->InheritsFrom(THStack::Class())) {
      fType = kObjectHStack;
      set = kTRUE;
      TH1 *hist = (TH1 *)((THStack *)obj)->GetHists()->First();
      fDim = hist->GetDimension();
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Chi-square", kFP_MCHIS);
      fMethodList->Select(kFP_MCHIS, kFALSE);
   } else if (obj->InheritsFrom(TTree::Class())) {
      fType = kObjectTree;
      set = kTRUE;
      TString variables, cuts;
      GetTreeVarsAndCuts(fDataSet, variables, cuts);
      fDim = 1;
      for ( int i = 0; i < variables.Length() && fDim <= 2; ++i )
         if ( ':' == variables[i] ) fDim += 1;
      // For any three  of dimension bigger than 2,  set the dimension
      // to 0,  as we cannot infer  the dimension from  the TF1s, it's
      // better to have 0 as reference.
      if ( fDim > 2 ) fDim = 0;
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Unbinned Likelihood", kFP_MUBIN);
      fMethodList->Select(kFP_MUBIN, kFALSE);
   } else if (obj->InheritsFrom(TH1::Class())){
      fType = kObjectHisto;
      set = kTRUE;
      fDim = ((TH1*)obj)->GetDimension();
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Chi-square", kFP_MCHIS);
      fMethodList->AddEntry("Binned Likelihood", kFP_MBINL);
      fMethodList->Select(kFP_MCHIS, kFALSE);
   } else if (obj->InheritsFrom(TMultiGraph::Class())) {
      fType = kObjectMultiGraph;
      set = kTRUE;
      fDim = 1;
      fMethodList->RemoveAll();
      fMethodList->AddEntry("Chi-square", kFP_MCHIS);
      fMethodList->Select(kFP_MCHIS, kFALSE);
      fRobustValue->SetState(kTRUE);
      fRobustValue->GetNumberEntry()->SetToolTipText("Set robust value");
   }

   // Depending on the dimension of the object, allow the
   // visualization of sliders.
   if ( fDim < 2 || fType == kObjectTree )
      fGeneral->HideFrame(fSliderYParent);
   else
      fGeneral->ShowFrame(fSliderYParent);

   if ( fDim < 1 || fType == kObjectTree )
      fGeneral->HideFrame(fSliderXParent);
   else
      fGeneral->ShowFrame(fSliderXParent);

   // And also, depending on the dimension, add predefined functions.
   if ( fDim == 1 ) {
      if ( !fTypeFit->FindEntry("Predef-1D") )
         fTypeFit->InsertEntry("Predef-1D", kFP_PRED1D, kFP_PREVFIT);
   } else {
      if ( fTypeFit->FindEntry("Predef-1D") )
         fTypeFit->RemoveEntry(kFP_PRED1D);
   }

   if ( fDim == 2 ) {
      if ( !fTypeFit->FindEntry("Predef-2D") )
         fTypeFit->InsertEntry("Predef-2D", kFP_PRED2D, kFP_PREVFIT);
   } else {
      if ( fTypeFit->FindEntry("Predef-2D") )
         fTypeFit->RemoveEntry(kFP_PRED2D);
   }

   return set;
}

////////////////////////////////////////////////////////////////////////////////
/// Show object name on the top.

void TFitEditor::ShowObjectName(TObject* obj)
{
   TString name;
   bool isTree = false;

   // Build the string to be compared to look for the object.
   if (obj) {
      name = obj->ClassName();
      name.Append("::");
      name.Append(obj->GetName());
      isTree = strcmp(obj->ClassName(), "TTree") == 0;
   } else {
      name = "No object selected";
   }
   fStatusBar->SetText(name.Data(),0);

   // If the selection was done in the fDataSet combo box, there is no need
   // to search through the list
   TGTextLBEntry* selectedEntry = static_cast<TGTextLBEntry*> ( fDataSet->GetSelectedEntry());
   if ( selectedEntry ) {
      TString selectedName = selectedEntry->GetText()->GetString();
      if ( isTree )
         selectedName = selectedName(0, selectedName.First(' '));
      if ( name.CompareTo(selectedName) == 0 ) {
         Layout();
         return;
      }
   }

   // Search through the list for the object
   Int_t entryId = kFP_NOSEL+1;
   bool found = false;
   while ( TGTextLBEntry* entry = static_cast<TGTextLBEntry*>
           ( fDataSet->GetListBox()->GetEntry(entryId)) ) {
      TString compareName = entry->GetText()->GetString();
      if ( isTree )
         compareName = compareName(0, compareName.First(' '));
      if ( name.CompareTo(compareName) == 0 ) {
         // If the object is found, select it
         fDataSet->Select(entryId, false);
         found = true;
         break;
      }
      entryId += 1;
   }

   // If the object was not found, add it and select it.
   if ( !found ) {
      fDataSet->AddEntry(name.Data(), entryId);
      fDataSet->Select(entryId, kTRUE);
   }

   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Get draw options of the selected object.

Option_t *TFitEditor::GetDrawOption() const
{
   if (!fParentPad) return "";

   TListIter next(fParentPad->GetListOfPrimitives());
   TObject *obj;
   while ((obj = next())) {
      if (obj == fFitObject) return next.GetOption();
   }
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// Set selected minimization library in use.

void TFitEditor::DoLibrary(Bool_t on)
{
   TGButton *bt = (TGButton *)gTQSender;
   Int_t id = bt->WidgetId();

   switch (id) {

      // Depending on the selected library, set the state of the rest
      // of the buttons.
      case kFP_LMIN:
         {
            if (on) {
               fLibMinuit->SetState(kButtonDown);
               fLibMinuit2->SetState(kButtonUp);
               fLibFumili->SetState(kButtonUp);
               if ( fLibGSL->GetState() != kButtonDisabled )
                  fLibGSL->SetState(kButtonUp);
               if ( fLibGenetics->GetState() != kButtonDisabled )
                  fLibGenetics->SetState(kButtonUp);
               fStatusBar->SetText("LIB Minuit", 1);
            }

         }
         break;

      case kFP_LMIN2:
         {
            if (on) {
               fLibMinuit->SetState(kButtonUp);
               fLibMinuit2->SetState(kButtonDown);
               fLibFumili->SetState(kButtonUp);
               if ( fLibGSL->GetState() != kButtonDisabled )
                  fLibGSL->SetState(kButtonUp);
               if ( fLibGenetics->GetState() != kButtonDisabled )
                  fLibGenetics->SetState(kButtonUp);
               fStatusBar->SetText("LIB Minuit2", 1);
            }
         }
         break;

      case kFP_LFUM:
         {
            if (on) {
               fLibMinuit->SetState(kButtonUp);
               fLibMinuit2->SetState(kButtonUp);
               fLibFumili->SetState(kButtonDown);
               if ( fLibGSL->GetState() != kButtonDisabled )
                  fLibGSL->SetState(kButtonUp);
               if ( fLibGenetics->GetState() != kButtonDisabled )
                  fLibGenetics->SetState(kButtonUp);
               fStatusBar->SetText("LIB Fumili", 1);
            }
         }
         break;
      case kFP_LGSL:
         {
            if (on) {
               fLibMinuit->SetState(kButtonUp);
               fLibMinuit2->SetState(kButtonUp);
               fLibFumili->SetState(kButtonUp);
               if ( fLibGSL->GetState() != kButtonDisabled )
                  fLibGSL->SetState(kButtonDown);
               if ( fLibGenetics->GetState() != kButtonDisabled )
                  fLibGenetics->SetState(kButtonUp);
               fStatusBar->SetText("LIB GSL", 1);
            }
         }
         break;
      case kFP_LGAS:
      {
         if (on) {
            fLibMinuit->SetState(kButtonUp);
            fLibMinuit2->SetState(kButtonUp);
            fLibFumili->SetState(kButtonUp);
            if ( fLibGSL->GetState() != kButtonDisabled )
               fLibGSL->SetState(kButtonUp);
            if ( fLibGenetics->GetState() != kButtonDisabled )
               fLibGenetics->SetState(kButtonDown);
            fStatusBar->SetText("LIB Genetics", 1);
         }
      }
      default:
         break;
   }
   FillMinMethodList();
}

////////////////////////////////////////////////////////////////////////////////
/// Set selected minimization method in use.

void TFitEditor::DoMinMethod(Int_t )
{
   if ( fMinMethodList->GetSelected() == kFP_MIGRAD )
      fStatusBar->SetText("MIGRAD",2);
   else if ( fMinMethodList->GetSelected() == kFP_FUMILI)
      fStatusBar->SetText("FUMILI",2);
   else if ( fMinMethodList->GetSelected() == kFP_SIMPLX )
      fStatusBar->SetText("SIMPLEX",2);
   else if ( fMinMethodList->GetSelected() == kFP_SCAN )
      fStatusBar->SetText("SCAN",2);
   else if ( fMinMethodList->GetSelected() == kFP_COMBINATION )
      fStatusBar->SetText("Combination",2);
   else if ( fMinMethodList->GetSelected() == kFP_GSLFR )
      fStatusBar->SetText("CONJFR",2);
   else if ( fMinMethodList->GetSelected() == kFP_GSLPR )
      fStatusBar->SetText("CONJPR",2);
   else if ( fMinMethodList->GetSelected() == kFP_BFGS )
      fStatusBar->SetText("BFGS",2);
   else if ( fMinMethodList->GetSelected() == kFP_BFGS2 )
      fStatusBar->SetText("BFGS2",2);
   else if ( fMinMethodList->GetSelected() == kFP_GSLLM )
      fStatusBar->SetText("GSLLM",2);
   else if ( fMinMethodList->GetSelected() == kFP_GSLSA)
      fStatusBar->SetText("SimAn",2);
   else if ( fMinMethodList->GetSelected() == kFP_TMVAGA )
      fStatusBar->SetText("TMVAGA",2);
   else if ( fMinMethodList->GetSelected() == kFP_GALIB )
      fStatusBar->SetText("GALIB",2);


}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum number of iterations.

void TFitEditor::DoMaxIterations()
{
   Long_t itr = fIterations->GetIntNumber();
   fStatusBar->SetText(Form("Itr: %ld",itr),2);
}

////////////////////////////////////////////////////////////////////////////////
/// Create section title in the GUI.

void TFitEditor::MakeTitle(TGCompositeFrame *parent, const char *title)
{
   TGCompositeFrame *ht = new TGCompositeFrame(parent, 350, 10,
                                               kFixedWidth | kHorizontalFrame);
   ht->AddFrame(new TGLabel(ht, title),
                new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   ht->AddFrame(new TGHorizontal3DLine(ht),
                new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 5, 5, 2, 2));
   parent->AddFrame(ht, new TGLayoutHints(kLHintsTop, 5, 0, 5, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Look in the list of function for TF1. If a TF1 is
/// found in the list of functions, it will be returned

TF1* TFitEditor::HasFitFunction()
{
   // Get the list of functions of the fit object
   TList *lf = GetFitObjectListOfFunctions();
   TF1* func = 0;

   // If it exists
   if ( lf ) {
      // Add the posibility to select previous fit function
      if ( !fTypeFit->FindEntry("Prev. Fit") )
         fTypeFit->InsertEntry("Prev. Fit",kFP_PREVFIT, kFP_UFUNC);

      // Then add all these functions to the fPrefFit structure.
      TObject *obj2;
      TIter next(lf, kIterForward);
      // Go over all the elements in lf
      while ((obj2 = next())) {
         if (obj2->InheritsFrom(TF1::Class())) {
            func = (TF1 *)obj2;
            fPrevFitIter it;
            // No go over all elements in fPrevFit
            for ( it = fPrevFit.begin(); it != fPrevFit.end(); ++it) {
               // To see wheather the object corresponds with fFitObject
               if ( it->first != fFitObject ) continue;
               // And if so, whether the function is already included
               if ( strcmp( func->GetName(), it->second->GetName() ) == 0 )
                  break;
               if ( strcmp( func->GetName(), "PrevFitTMP" ) == 0 )
                  break;
            }
            // Only if the function is not already in fPrevFit, the
            // breaks in the loops would make it to be different to
            // fPrevFit.end() if the function is already stored
            if ( it == fPrevFit.end() ) {
               fPrevFit.emplace(fFitObject, copyTF1(func));
            }
         }
      }

      // Select the PrevFit set
      fTypeFit->Select(kFP_PREVFIT);
      // And fill the function list
      FillFunctionList();
      fDrawAdvanced->SetState(kButtonUp);


   } else {
      // If there is no prev fit functions.
      fTypeFit->Select(kFP_UFUNC);
      // Call FillFunctionList as it might happen that the user is
      // changing from a TTree to another one, and thus the fFuncList
      // if not properly filled
      FillFunctionList();
   }

   fDrawAdvanced->SetState(kButtonDisabled);

   return func;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the fitting options from all the widgets.

void TFitEditor::RetrieveOptions(Foption_t& fitOpts, TString& drawOpts, ROOT::Math::MinimizerOptions& minOpts, Int_t /*npar */)
{
   drawOpts = "";

   fitOpts.Range    = (fUseRange->GetState() == kButtonDown);
   fitOpts.Integral = (fIntegral->GetState() == kButtonDown);
   fitOpts.More     = (fImproveResults->GetState() == kButtonDown);
   fitOpts.Errors   = (fBestErrors->GetState() == kButtonDown);
   fitOpts.Like = (fMethodList->GetSelected() != kFP_MCHIS);

   if (fEmptyBinsWghts1->GetState() == kButtonDown)
      fitOpts.W1 = 2;
   else if (fAllWeights1->GetState() == kButtonDown)
      fitOpts.W1 = 1;

   TString tmpStr = fEnteredFunc->GetText();
   if ( !(fLinearFit->GetState() == kButtonDown) &&
        (tmpStr.Contains("pol") || tmpStr.Contains("++")) )
      fitOpts.Minuit = 1;

   // if ( (int) fFuncPars.size() == npar )
   //    for ( Int_t i = 0; i < npar; ++i )
   //       if ( fFuncPars[i][PAR_MIN] != fFuncPars[i][PAR_MAX] )
   //

   //          //fitOpts.Bound = 1;
   //          break;
   //       }

   if (fChangedParams) {
      //std::cout << "Params have changed setting the Bound option " << std::endl;
      fitOpts.Bound = 1;
      fChangedParams = kFALSE;  // reset
   }

   //fitOpts.Nochisq  = (fNoChi2->GetState() == kButtonDown);
   fitOpts.Nostore  = (fNoStoreDrawing->GetState() == kButtonDown);
   fitOpts.Nograph  = (fNoDrawing->GetState() == kButtonDown);
   fitOpts.Plus     = (fAdd2FuncList->GetState() == kButtonDown);
   fitOpts.Gradient = (fUseGradient->GetState() == kButtonDown);
   fitOpts.Quiet    = ( fOptQuiet->GetState() == kButtonDown );
   fitOpts.Verbose  = ( fOptVerbose->GetState() == kButtonDown );

   if ( !(fType != kObjectGraph) && (fEnableRobust->GetState() == kButtonDown) )
   {
      fitOpts.Robust = 1;
      fitOpts.hRobust = fRobustValue->GetNumber();
   }

   drawOpts = GetDrawOption();

   if ( fLibMinuit->GetState() == kButtonDown )
      minOpts.SetMinimizerType ( "Minuit");
   else if ( fLibMinuit2->GetState() == kButtonDown)
      minOpts.SetMinimizerType ( "Minuit2" );
   else if ( fLibFumili->GetState() == kButtonDown )
      minOpts.SetMinimizerType ("Fumili" );
   else if ( fLibGSL->GetState() == kButtonDown )
      minOpts.SetMinimizerType ("GSLMultiMin" );

   if ( fMinMethodList->GetSelected() == kFP_MIGRAD )
      minOpts.SetMinimizerAlgorithm( "Migrad" );
   else if ( fMinMethodList->GetSelected() == kFP_FUMILI)
      if ( fLibMinuit2->GetState() == kButtonDown )
         minOpts.SetMinimizerAlgorithm( "Fumili2" );
      else
         minOpts.SetMinimizerAlgorithm( "Fumili" );
   else if ( fMinMethodList->GetSelected() == kFP_SIMPLX )
      minOpts.SetMinimizerAlgorithm( "Simplex" );
   else if ( fMinMethodList->GetSelected() == kFP_SCAN )
      minOpts.SetMinimizerAlgorithm( "Scan" );
   else if ( fMinMethodList->GetSelected() == kFP_COMBINATION )
      minOpts.SetMinimizerAlgorithm( "Minimize" );
   else if ( fMinMethodList->GetSelected() == kFP_GSLFR )
      minOpts.SetMinimizerAlgorithm( "conjugatefr" );
   else if ( fMinMethodList->GetSelected() == kFP_GSLPR )
      minOpts.SetMinimizerAlgorithm( "conjugatepr" );
   else if ( fMinMethodList->GetSelected() == kFP_BFGS )
      minOpts.SetMinimizerAlgorithm( "bfgs" );
   else if ( fMinMethodList->GetSelected() == kFP_BFGS2 )
      minOpts.SetMinimizerAlgorithm( "bfgs2" );
   else if ( fMinMethodList->GetSelected() == kFP_GSLLM ) {
      minOpts.SetMinimizerType ("GSLMultiFit" );
      minOpts.SetMinimizerAlgorithm( "" );
   } else if ( fMinMethodList->GetSelected() == kFP_GSLSA) {
      minOpts.SetMinimizerType ("GSLSimAn" );
      minOpts.SetMinimizerAlgorithm( "" );
   } else if ( fMinMethodList->GetSelected() == kFP_TMVAGA) {
      minOpts.SetMinimizerType ("Geneti2c" );
      minOpts.SetMinimizerAlgorithm( "" );
   } else if ( fMinMethodList->GetSelected() == kFP_GALIB) {
      minOpts.SetMinimizerType ("GAlibMin" );
      minOpts.SetMinimizerAlgorithm( "" );
   }

   minOpts.SetErrorDef ( fErrorScale->GetNumber() );
   minOpts.SetTolerance( fTolerance->GetNumber() );
   minOpts.SetMaxIterations(fIterations->GetIntNumber());
   minOpts.SetMaxFunctionCalls(fIterations->GetIntNumber());
}

void TFitEditor::SetEditable(Bool_t state)
{
   // Set the state of some input widgets depending on whether the fit
   // function can be defined by text or if it is an existing one.
   if ( state )
   {
      fEnteredFunc-> SetState(kTRUE);
      fAdd        -> SetState(kButtonUp,  kFALSE);
      fNormAdd    -> SetState(kButtonUp,  kFALSE);
      fConv       -> SetState(kButtonUp,  kFALSE);
      fNone       -> SetState(kButtonDown,kFALSE); // fNone::State is the one used as reference
   }
   else
   {
      fEnteredFunc-> SetState(kFALSE);
      fAdd        -> SetState(kButtonDisabled, kFALSE);
      fNormAdd    -> SetState(kButtonDisabled, kFALSE);
      fConv       -> SetState(kButtonDisabled, kFALSE);
      fNone       -> SetState(kButtonDisabled, kFALSE);
   }
}

void TFitEditor::GetRanges(ROOT::Fit::DataRange& drange)
{
   // Return the ranges selected by the sliders.

   // It's not working for trees as they don't have TAxis.
   if ( fType == kObjectTree ) return;

   if ( fType != kObjectTree ) {
      Int_t ixmin = (Int_t)(fSliderX->GetMinPosition());
      Int_t ixmax = (Int_t)(fSliderX->GetMaxPosition());
      Double_t xmin = fXaxis->GetBinLowEdge(ixmin);
      Double_t xmax = fXaxis->GetBinUpEdge(ixmax);
      drange.AddRange(0,xmin, xmax);
   }

   if ( fDim > 1 ) {
      assert(fYaxis);
      Int_t iymin = (Int_t)(fSliderY->GetMinPosition());
      Int_t iymax = (Int_t)(fSliderY->GetMaxPosition());
      Double_t ymin = fYaxis->GetBinLowEdge(iymin);
      Double_t ymax = fYaxis->GetBinUpEdge(iymax);
      drange.AddRange(1,ymin, ymax);
   }
   if ( fDim > 2 ) {
      assert(fZaxis);
      Int_t izmin = (Int_t)(fSliderZ->GetMinPosition());
      Int_t izmax = (Int_t)(fSliderZ->GetMaxPosition());
      Double_t zmin = fZaxis->GetBinLowEdge(izmin);
      Double_t zmax = fZaxis->GetBinUpEdge(izmax);
      drange.AddRange(2,zmin, zmax);
   }
}

TList* TFitEditor::GetFitObjectListOfFunctions()
{
   // Get the list of functions previously used in the fitobject.

   TList *listOfFunctions = 0;
   if ( fFitObject ) {
      switch (fType) {

      case kObjectHisto:
         listOfFunctions = ((TH1 *)fFitObject)->GetListOfFunctions();
         break;

      case kObjectGraph:
         listOfFunctions = ((TGraph *)fFitObject)->GetListOfFunctions();
         break;

      case kObjectMultiGraph:
         listOfFunctions = ((TMultiGraph *)fFitObject)->GetListOfFunctions();
         break;

      case kObjectGraph2D:
         listOfFunctions = ((TGraph2D *)fFitObject)->GetListOfFunctions();
         break;

      case kObjectHStack:
      case kObjectTree:
      default:
         break;
      }
   }
   return listOfFunctions;
}

void TFitEditor::GetFunctionsFromSystem()
{
   // Looks for all the functions registered in the current ROOT
   // session.

   // First, clean the copies stored in fSystemFunc
   for (auto func : fSystemFuncs)
      delete func;

   fSystemFuncs.clear();

   // Be carefull not to store functions that will be in the
   // predefined section
   const unsigned int nfuncs = 16;
   const char* fnames[nfuncs] = { "gaus" ,   "gausn", "expo", "landau",
                                  "landaun", "pol0",  "pol1", "pol2",
                                  "pol3",    "pol4",  "pol5", "pol6",
                                  "pol7",    "pol8",  "pol9", "user"
   };

   // No go through all the objects registered in gROOT
   TIter functionsIter(gROOT->GetListOfFunctions());
   TObject* obj;
   while( ( obj = (TObject*) functionsIter() ) ) {
      // And if they are TF1s
      if ( TF1* func = dynamic_cast<TF1*>(obj) ) {
         bool addFunction = true;
         // And they are not already registered in fSystemFunc
         for ( unsigned int i = 0; i < nfuncs; ++i ) {
            if ( strcmp( func->GetName(), fnames[i] ) == 0 ) {
               addFunction = false;
               break;
            }
         }
         // Add them.
         if ( addFunction )
            fSystemFuncs.emplace_back( copyTF1(func) );
      }
   }
}

TList* TFitEditor::GetListOfFittingFunctions(TObject* obj)
{
   // This function returns a TList with all the functions used in the
   // FitPanel to fit a given object. If the object passed is NULL,
   // then the object used is the currently selected one. It is
   // important to notice that the FitPanel is still the owner of
   // those functions. This means that the user SHOULD NOT delete any
   // of these functions, as the FitPanel will do so in the
   // destructor.

   if (!obj) obj = fFitObject;

   TList *retList = new TList();

   std::pair<fPrevFitIter, fPrevFitIter> look = fPrevFit.equal_range(obj);
   for ( fPrevFitIter it = look.first; it != look.second; ++it ) {
      retList->Add(it->second);
   }

   return retList;
}

TF1* TFitEditor::GetFitFunction()
{
   // Get the fit function selected or declared in the fiteditor

   TF1 *fitFunc = 0;
   // If the function is not editable ==> it means it is registered in
   // gROOT
   if ( fNone->GetState() == kButtonDisabled )
   {
      // So we find it
      TF1* tmpF1 = FindFunction();
      // And if we don't find it, then it means there is something wrong!
      if ( tmpF1 == 0 )
      {
               new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                            "Error...", "1) Verify the entered function string!",
                            kMBIconStop,kMBOk, 0);
               return 0;
      }

      // Now we make a copy that will be used temporary. The caller of
      // the function should delete the returned function.
      fitFunc = (TF1*)tmpF1->IsA()->New();
      tmpF1->Copy(*fitFunc);
      // Copy the parameters of the function, if and only if the
      // parameters stored does not correspond with the ones of these
      // functions. Perhaps the user has already called
      // DoSetParameters. There is no way to know whether the
      // parameters have been modified, so we check the size of
      // fFuncPars against number of parameters.
      if ( int(fFuncPars.size()) != tmpF1->GetNpar() )
      {
         fitFunc->SetParameters(tmpF1->GetParameters());
         GetParameters(fFuncPars, fitFunc);
      } else {
         SetParameters(fFuncPars, fitFunc);
      }
   }

   // If, we have no function at this point, it means that is is
   // described in fEnteredFunc, so we create it from scratch.
   if ( fitFunc == 0 )
   {
      ROOT::Fit::DataRange drange;
      GetRanges(drange);
      double xmin, xmax, ymin, ymax, zmin, zmax;
      drange.GetRange(xmin, xmax, ymin, ymax, zmin, zmax);

      // Depending of course on the number of dimensions the object
      // has. These commands will raise an error message if the user
      // has not defined the function properly
      if ( fDim == 1 || fDim == 0 )
      {

         fitFunc = new TF1("PrevFitTMP",fEnteredFunc->GetText(), xmin, xmax );
         //std::cout << "GetFitFunction - created function PrevFitTMP " << fEnteredFunc->GetText() << "  " << fitFunc << std::endl;
         if (fNormAdd->IsOn())
         {
            if (fSumFunc) delete fSumFunc;
            fSumFunc = new TF1NormSum(fEnteredFunc->GetText(), xmin, xmax);
            fitFunc  = new TF1("PrevFitTMP", *fSumFunc, xmin, xmax, fSumFunc->GetNpar());
            for (int i = 0; i < fitFunc->GetNpar(); ++i) fitFunc->SetParName(i, fSumFunc->GetParName(i) );
            //std::cout << "create fit normalized function " << fSumFunc << " fitfunc " << fitFunc << std::endl;
         }

         if (fConv -> IsOn())
         {
            if (fConvFunc) delete fConvFunc;
            fConvFunc = new TF1Convolution(fEnteredFunc->GetText());
            fitFunc  = new TF1("PrevFitTMP", *fConvFunc, xmin, xmax, fConvFunc->GetNpar());
            for (int i = 0; i < fitFunc->GetNpar(); ++i) fitFunc->SetParName(i, fConvFunc->GetParName(i) );
            //std::cout << "create fit convolution function " << fSumFunc << " fitfunc " << fitFunc << std::endl;
         }
      }
      else if ( fDim == 2 ) {
         fitFunc = new TF2("PrevFitTMP",fEnteredFunc->GetText(), xmin, xmax, ymin, ymax );
      }
      else if ( fDim == 3 ) {
         fitFunc = new TF3("PrevFitTMP",fEnteredFunc->GetText(), xmin, xmax, ymin, ymax, zmin, zmax );
      }

      // if the function is not a C defined
      if ( fNone->GetState() != kButtonDisabled )
      {
         // and the formulas are the same
         TF1* tmpF1 = FindFunction();
//         if (tmpF1)
            //std::cout << "GetFitFunction: found existing function " << tmpF1 << "  " << tmpF1->GetName() << "  " << tmpF1->GetExpFormula() << std::endl;
//         else
            //std::cout << "GetFitFunction: - no existing function  found " << std::endl;
         if ( tmpF1 != 0 && fitFunc != 0 &&
              strcmp(tmpF1->GetExpFormula(), fEnteredFunc->GetText()) == 0 ) {
            // copy everything from the founction available in gROOT
            //std::cout << "GetFitFunction: copying tmp function in PrevFitTMP " <<  tmpF1->GetName()  << "  "
            //          << tmpF1->GetExpFormula() << std::endl;
            tmpF1->Copy(*fitFunc);
            if ( int(fFuncPars.size()) != tmpF1->GetNpar() )
            {
               GetParameters(fFuncPars, fitFunc);
            }
         }
      }
   }

   return fitFunc;
}
