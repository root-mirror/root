// @(#)Root/tmva $Id$
// Author: Kim Albertsson
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVAGUI                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors :                                                                      *
 *      Kim Albertsson  <kim.albertsson@cern.ch> - LTU & CERN                     *
 *                                                                                *
 * Copyright (c) 2005-2017:                                                       *
 *      CERN, Switzerland                                                         *
 *      LTU, Sweden                                                               *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TMVA/efficienciesMulticlass.h"

// TMVA
#include "TMVA/Config.h"
#include "TMVA/tmvaglob.h"

// ROOT
#include "TControlBar.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TIterator.h"
#include "TKey.h"
#include "TROOT.h"

// STL
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
///
/// Note: This file assumes a certain structure on the input file. The structure
/// is as follows:
///
/// - dataset (TDirectory)
/// - ... some variables, plots ...
/// - Method_XXX (TDirectory)
///   + XXX (TDirectory)
///      * ... some plots ...
///      * MVA_Method_XXX_Test_#classname#
///      * MVA_Method_XXX_Train_#classname#
///      * ... some plots ...
/// - Method_YYY (TDirectory)
///   + YYY (TDirectory)
///      * ... some plots ...
///      * MVA_Method_YYY_Test_#classname#
///      * MVA_Method_YYY_Train_#classname#
///      * ... some plots ...
/// - TestTree (TTree)
///   + ... data...
/// - TrainTree (TTree)
///   + ... data...
///
/// Keeping this in mind makes the main loop in getRocCurves easier to follow :)
///

////////////////////////////////////////////////////////////////////////////////
/// Private class that simplify drawing plots combining information from
/// several methods.
///
/// Each wrapper will manage a canvas and a legend and provide convenience
/// functions to add data to these. It also provides a save function for
/// saving an image representation to disk.
///
/// Feel free to extend this class as you see fit. It is intended as a
/// convenience when showing multiclass roccurves, not a fully general tool.
///
/// Usage:
///   auto p = new EfficiencyPlotWrapper(name, title, dataset, i):
///   for (TGraph * g : listOfGraphs) {
///      p->AddGraph(g);
///      p->AddLegendEntry(methodName);
///   }
///   p->save();
///

class EfficiencyPlotWrapper {
public:
   TCanvas *fCanvas;
   TLegend *fLegend;

   TString fDataset;

   Int_t fColor;

   UInt_t fNumMethods;

   EfficiencyPlotWrapper(TString name, TString title, TString dataset, size_t i);

   Int_t addGraph(TGraph *graph);
   void addLegendEntry(TString methodTitle, TGraph *graph);

   void save();

private:
   Float_t fx0L;
   Float_t fdxL;
   Float_t fy0H;
   Float_t fdyH;

   TCanvas *newEfficiencyCanvas(TString name, TString title, size_t i);
   TLegend *newEfficiencyLegend();
};

using classcanvasmap_t = std::map<TString, EfficiencyPlotWrapper *>;
using roccurvelist_t = std::vector<std::tuple<TString, TString, TGraph *>>;

// Constants
const char *BUTTON_TYPE = "button";

// Private functions
namespace TMVA {
std::vector<TString> getclassnames(TString dataset, TString fin);
roccurvelist_t getRocCurves(TDirectory *binDir, TString methodPrefix, TString graphNameRef);
void plotEfficienciesMulticlass(roccurvelist_t rocCurves, classcanvasmap_t classCanvasMap);
}

////////////////////////////////////////////////////////////////////////////////
/// Private (helper) functions - Implementation
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
///

std::vector<TString> TMVA::getclassnames(TString dataset, TString fin)
{
   TFile *file = TMVA::TMVAGlob::OpenFile(fin);
   TDirectory *dir = (TDirectory *)file->GetDirectory(dataset)->GetDirectory("InputVariables_Id");
   if (!dir) {
      std::cout << "Could not locate directory '" << dataset << "/InputVariables_Id' in file: " << fin << std::endl;
      return {};
   }

   auto classnames = TMVA::TMVAGlob::GetClassNames(dir);
   return classnames;
}

////////////////////////////////////////////////////////////////////////////////
///

roccurvelist_t TMVA::getRocCurves(TDirectory *binDir, TString methodPrefix, TString graphNameRef)
{
   roccurvelist_t rocCurves;

   TList methods;
   UInt_t nm = TMVAGlob::GetListOfMethods(methods, binDir);
   if (nm == 0) {
      cout << "ups .. no methods found in to plot ROC curve for ... give up" << endl;
      return rocCurves;
   }
   //   TIter next(file->GetListOfKeys());
   TIter next(&methods);

   // Loop over all method categories
   TKey *key;
   while ((key = (TKey *)next())) {
      TDirectory *mDir = (TDirectory *)key->ReadObj();
      TList titles;
      TMVAGlob::GetListOfTitles(mDir, titles);

      // Loop over each method within a category
      TIter nextTitle(&titles);
      TKey *titkey;
      TDirectory *titDir;
      while ((titkey = TMVAGlob::NextKey(nextTitle, "TDirectory"))) {
         titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle, titDir);

         // Loop through all plots for the method
         TIter nextKey(titDir->GetListOfKeys());
         TKey *hkey2;
         while ((hkey2 = TMVAGlob::NextKey(nextKey, "TGraph"))) {

            TGraph *h = (TGraph *)hkey2->ReadObj();
            TString hname = h->GetName();
            if (hname.Contains(graphNameRef) && hname.BeginsWith(methodPrefix) && !hname.Contains("Train")) {

               // Extract classname from plot name
               // classname is string after nameref
               Int_t index = hname.Index(graphNameRef) + graphNameRef.Length();
               TString classname = hname(index, hname.Length() - index);

               //std::cout << "Found TGraph " << hname << " with classname " << classname << std::endl;

               rocCurves.push_back(std::make_tuple(methodTitle, classname, h));
            }
         }
      }
   }
   return rocCurves;
}

////////////////////////////////////////////////////////////////////////////////
/// Public functions - Implementation
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Private convenience function.
///
/// Adds a given a list of roc curves provided as n-tuple on the form
///    (methodname, classname, graph)
/// to the canvas corresponding to the classname.
///

void TMVA::plotEfficienciesMulticlass(roccurvelist_t rocCurves, classcanvasmap_t classCanvasMap)
{
   for (auto &item : rocCurves) {

      TString methodTitle = std::get<0>(item);
      TString classname = std::get<1>(item);
      TGraph *h = std::get<2>(item);

      try {
         EfficiencyPlotWrapper *plotWrapper = classCanvasMap.at(classname);
         plotWrapper->addGraph(h);
         plotWrapper->addLegendEntry(methodTitle, h);
      } catch (const std::out_of_range &) {
         cout << Form("ERROR: Class %s discovered among plots but was not found by TMVAMulticlassGui. Skipping.",
                      classname.Data())
              << endl;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Entry point. Called from the TMVAMulticlassGui Buttons
///
/// \param dataset        Dataset to operate on. Should be created by the TMVA Multiclass Factory.
/// \param filename_input Name of the input file procuded by a TMVA Multiclass Factory.
/// \param plotType       Specified what kind of ROC curve to draw. Currently only rejB vs. effS is supported.

void TMVA::efficienciesMulticlass1vsRest(TString dataset, TString filename_input, EEfficiencyPlotType plotType,
                                         Bool_t useTMVAStyle)
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize(useTMVAStyle);
   plotEfficienciesMulticlass1vsRest(dataset, plotType, filename_input);
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Work horse function. Will operate on the currently open file (opened by
/// efficienciesMulticlass).
///
/// \param plotType See effcienciesMulticlass.
/// \param binDir   Directory in the file on which to operate.

void TMVA::plotEfficienciesMulticlass1vsRest(TString dataset, EEfficiencyPlotType plotType, TString filename_input)
{
   // The current multiclass version implements only type 2 - rejB vs effS
   if (plotType != EEfficiencyPlotType::kRejBvsEffS) {
      std::cout << "For multiclass, only rejB vs effS is currently implemented.";
      return;
   }

   // checks if filename_input is already open, and if not opens one
   TFile *file = TMVAGlob::OpenFile(filename_input);
   if (file == nullptr) {
      std::cout << "ERROR: filename \"" << filename_input << "\" is not found.";
      return;
   }
   auto binDir = file->GetDirectory(dataset.Data());

   size_t iPlot = 0;
   auto classnames = getclassnames(dataset, filename_input);
   TString methodPrefix = "MVA_";
   TString graphNameRef = "_rejBvsS_";

   classcanvasmap_t classCanvasMap;
   for (auto &classname : classnames) {
      TString name = Form("roc_%s_vs_rest", classname.Data());
      TString title = Form("ROC Curve %s vs rest", classname.Data());
      EfficiencyPlotWrapper *plotWrapper = new EfficiencyPlotWrapper(name, title, dataset, iPlot++);
      classCanvasMap.emplace(classname.Data(), plotWrapper);
   }

   roccurvelist_t rocCurves = getRocCurves(binDir, methodPrefix, graphNameRef);
   plotEfficienciesMulticlass(rocCurves, classCanvasMap);

   for (auto const &item : classCanvasMap) {
      auto plotWrapper = item.second;
      plotWrapper->save();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Entry point. Called from the TMVAMulticlassGui Buttons
///
/// \param dataset
/// \param fin

void TMVA::efficienciesMulticlass1vs1(TString dataset, TString fin)
{
   std::cout << "--- Running Roc1v1Gui for input file: " << fin << std::endl;

   TMVAGlob::Initialize();

   // create the control bar
   TString title = "1v1 ROC curve comparison";
   TControlBar *cbar = new TControlBar("vertical", title, 50, 50);

   gDirectory->pwd();
   auto classnames = getclassnames(dataset, fin);

   // configure buttons
   for (auto &classname : classnames) {
      cbar->AddButton(Form("Class: %s", classname.Data()),
                      Form("TMVA::plotEfficienciesMulticlass1vs1(\"%s\", \"%s\", \"%s\")", dataset.Data(), fin.Data(),
                           classname.Data()),
                      BUTTON_TYPE);
   }

   cbar->SetTextColor("blue");
   cbar->Show();

   gROOT->SaveContext();
}

////////////////////////////////////////////////////////////////////////////////
/// Generates K-1 plots comparing a given base class against all others (except
/// itself). For each plot, the base class is considered signal and the other
/// class is considered background.
///
/// Given 3 classes in the dataset and providing "Class 0" as the base class
/// this would generate 2 plots comparing
///    - Class 0 vs Class 1, and
///    - Class 0 vs Class 2.
/// For the "Class 0 vs Class 1" plot, events from Class 2 are ignored. For the
/// "Class 0 vs Class 2" plot, events from Class 1 are ignored.
///
/// \param dataset
/// \param fin
/// \param baseClassname name of the class which will be considered signal

void TMVA::plotEfficienciesMulticlass1vs1(TString dataset, TString fin, TString baseClassname)
{

   TMVAGlob::Initialize();

   auto classnames = getclassnames(dataset, fin);
   size_t iPlot = 0;

   TString methodPrefix = "MVA_";
   TString graphNameRef = Form("_1v1rejBvsS_%s_vs_", baseClassname.Data());

   TFile *file = TMVAGlob::OpenFile(fin);
   if (file == nullptr) {
      std::cout << "ERROR: filename \"" << fin << "\" is not found.";
      return;
   }
   auto binDir = file->GetDirectory(dataset.Data());

   classcanvasmap_t classCanvasMap;
   for (auto &classname : classnames) {

      if (classname == baseClassname) {
         continue;
      }

      TString name = Form("1v1roc_%s_vs_%s", baseClassname.Data(), classname.Data());
      TString title = Form("ROC Curve %s (Sig) vs %s (Bkg)", baseClassname.Data(), classname.Data());
      EfficiencyPlotWrapper *plotWrapper = new EfficiencyPlotWrapper(name, title, dataset, iPlot++);
      classCanvasMap.emplace(classname.Data(), plotWrapper);
   }

   roccurvelist_t rocCurves = getRocCurves(binDir, methodPrefix, graphNameRef);
   plotEfficienciesMulticlass(rocCurves, classCanvasMap);

   for (auto const &item : classCanvasMap) {
      auto plotWrapper = item.second;
      plotWrapper->save();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Private class EfficiencyPlotWrapper - Implementation
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Constructs a new canvas + auxiliary data for showing an efficiency plot.
///

EfficiencyPlotWrapper::EfficiencyPlotWrapper(TString name, TString title, TString dataset, size_t i)
{
   // Legend extents (init before calling newEfficiencyLegend...)
   fx0L = 0.107;
   fy0H = 0.899;
   fdxL = 0.457 - fx0L;
   fdyH = 0.22;
   fx0L = 0.15;
   fy0H = 1 - fy0H + fdyH + 0.07;

   fColor = 1;
   fNumMethods = 0;

   fDataset = dataset;

   fCanvas = newEfficiencyCanvas(name, title, i);
   fLegend = newEfficiencyLegend();
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a new graph to the plot. The added graph should contain a single ROC
/// curve.
///

Int_t EfficiencyPlotWrapper::addGraph(TGraph *graph)
{
   graph->SetLineWidth(3);
   graph->SetLineColor(fColor);
   fColor++;
   if (fColor == 5 || fColor == 10 || fColor == 11) {
      fColor++;
   }

   fCanvas->cd();
   graph->Draw("");
   fCanvas->Update();

   ++fNumMethods;

   return fColor;
}

////////////////////////////////////////////////////////////////////////////////
/// WARNING: Uses the current color, thus the correct call ordering is:
///   plotWrapper->addGraph(...);
///   plotWrapper->addLegendEntry(...);
///

void EfficiencyPlotWrapper::addLegendEntry(TString methodTitle, TGraph *graph)
{
   fLegend->AddEntry(graph, methodTitle, "l");

   Float_t dyH_local = fdyH * (Float_t(TMath::Min((UInt_t)10, fNumMethods) - 3.0) / 4.0);
   fLegend->SetY2(fy0H + dyH_local);

   fLegend->Paint();
   fCanvas->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to create new Canvas
///
/// \param name  Name...
/// \param title Title to be displayed on canvas
/// \param i     Index to offset a collection of canvases from each other
///

TCanvas *EfficiencyPlotWrapper::newEfficiencyCanvas(TString name, TString title, size_t i)
{
   TCanvas *c = new TCanvas(name, title, 200 + i * 50, 0 + i * 50, 650, 500);
   // global style settings
   c->SetGrid();
   c->SetTicks();

   // Frame
   TString xtit = "Signal Efficiency";
   TString ytit = "Background Rejection (1 - eff)";
   Double_t x1 = 0.0;
   Double_t x2 = 1.0;
   Double_t y1 = 0.0;
   Double_t y2 = 1.0;

   TH1F *frame = new TH1F(Form("%s_%s", title.Data(), "frame"), title, 500, x1, x2);
   frame->SetMinimum(y1);
   frame->SetMaximum(y2);

   frame->GetXaxis()->SetTitle(xtit);
   frame->GetYaxis()->SetTitle(ytit);
   TMVA::TMVAGlob::SetFrameStyle(frame, 1.0);
   frame->Draw();

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to create new legend.

TLegend *EfficiencyPlotWrapper::newEfficiencyLegend()
{
   TLegend *legend = new TLegend(fx0L, fy0H - fdyH, fx0L + fdxL, fy0H);
   // legend->SetTextSize( 0.05 );
   legend->SetHeader("MVA Method:");
   legend->SetMargin(0.4);
   legend->Draw("");

   return legend;
}

////////////////////////////////////////////////////////////////////////////////
/// Saves the current state of the plot to disk.
///

void EfficiencyPlotWrapper::save()
{
   TString fname = fDataset + "/plots/" + fCanvas->GetName();
   TMVA::TMVAGlob::imgconv(fCanvas, fname);
}