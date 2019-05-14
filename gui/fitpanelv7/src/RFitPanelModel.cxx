/// \file RFitPanelModel.cxx
/// \ingroup WebGui ROOT7
/// \author Sergey Linev <S.Linev@gsi.de>
/// \author Iliana Betsou <Iliana.Betsou@cern.ch>
/// \date 2019-04-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RFitPanelModel.hxx>

#include "TH1.h"
#include "TDirectory.h"
#include "TPluginManager.h"

#include "TF1.h"
#include "TF2.h"


enum EFitPanel {
   kFP_NONE = 0,
   kFP_MIGRAD, kFP_SIMPLX, kFP_SCAN, kFP_COMBINATION,
   kFP_FUMILI, kFP_FUMILI2, kFP_GSLFR, kFP_GSLPR,
   kFP_BFGS, kFP_BFGS2, kFP_GSLLM, kFP_GSLSA,
   kFP_GALIB, kFP_TMVAGA
};

using namespace std::string_literals;


void ROOT::Experimental::RFitPanelModel::RFitFuncParsList::Clear()
{
   pars.clear();
   name.clear();
   haspars = false;
}

void ROOT::Experimental::RFitPanelModel::RFitFuncParsList::GetParameters(TF1 *func)
{
   pars.clear();
   haspars = true;

   for (int n = 0; n < func->GetNpar(); ++n) {
      pars.emplace_back(n, func->GetParName(n));
      auto &par = pars.back();

      par.value = std::to_string(func->GetParameter(n));
      par.error = std::to_string(func->GetParError(n));
      double min, max;
      func->GetParLimits(n, min, max);
      par.min = std::to_string(min);
      par.max = std::to_string(max);

      if ((min >= max) && ((min != 0) || (max != 0)))
         par.fixed = true;
   }
}

void ROOT::Experimental::RFitPanelModel::RFitFuncParsList::SetParameters(TF1 *func)
{
   if (func->GetNpar() != (int) pars.size()) {
      ::Error("RFitFuncParsList::SetParameters", "Mismatch in parameters numbers");
      return;
   }

   for (int n = 0; n < func->GetNpar(); ++n) {
      if (pars[n].name.compare(func->GetParName(n)) != 0) {
         ::Error("RFitFuncParsList::SetParameters", "Mismatch in parameter %d name %s %s", n, pars[n].name.c_str(), func->GetParName(n));
         return;
      }

      func->SetParameter(n, std::stod(pars[n].value));
      func->SetParError(n, std::stod(pars[n].error));
      if (pars[n].fixed) {
         func->FixParameter(n, std::stod(pars[n].value));
      } else {
         func->ReleaseParameter(n);
         double min = std::stod(pars[n].min);
         double max = std::stod(pars[n].max);
         if (min < max)
            func->SetParLimits(n, min, max);
      }
   }
}

///////////////////////////////

TH1* ROOT::Experimental::RFitPanelModel::GetSelectedHistogram(TH1 *hist)
{
   if (fSelectedData == "__hist__") return hist;
   if ((fSelectedData.compare(0,6,"gdir::") != 0) || !gDirectory) return nullptr;

   std::string hname = fSelectedData.substr(6);

   return dynamic_cast<TH1*> (gDirectory->GetList()->FindObject(hname.c_str()));
}


// Configure usage of histogram

bool ROOT::Experimental::RFitPanelModel::SelectHistogram(const std::string &hname, TH1 *hist)
{

   std::string histid;

   fDataSet.clear();
   TH1 *selected = nullptr;

   if (gDirectory) {
      TIter iter(gDirectory->GetList());
      TObject *item = nullptr;

       while ((item = iter()) != nullptr)
         if (item->InheritsFrom(TH1::Class())) {
            std::string dataid = "gdir::"s + item->GetName();

            if (hist && (hist == item)) {
               histid = dataid;
               selected = hist;
            } else if (!hname.empty() && hname.compare(item->GetName())) {
               histid = dataid;
               selected = dynamic_cast<TH1 *> (item);
            }
            fDataSet.emplace_back(dataid, Form("%s::%s", item->ClassName(), item->GetName()));
         }
   }

   if (hist && histid.empty()) {
      selected = hist;
      histid = "__hist__";
      fDataSet.emplace_back(histid, Form("%s::%s", hist->ClassName(), hist->GetName()));
   }

   fSelectedData = histid;

   UpdateRange(selected);

   UpdateFuncList();

   UpdateAdvanced(nullptr);

   return selected != nullptr;
}

void ROOT::Experimental::RFitPanelModel::UpdateRange(TH1 *hist)
{
   fDim = hist ? hist->GetDimension() : 0;

   fMinRangeX = 0.;
   fMaxRangeX = 100.;
   fMinRangeY = 0.;
   fMaxRangeY = 100.;

   if (hist && (fDim > 0)) {
      fMinRangeX = hist->GetXaxis()->GetXmin();
      fMaxRangeX = hist->GetXaxis()->GetXmax();
   }
   if (hist && (fDim > 1)) {
      fMinRangeY = hist->GetYaxis()->GetXmin();
      fMaxRangeY = hist->GetYaxis()->GetXmax();
   }

   // defined values
   fStepX = (fMaxRangeX - fMinRangeX) / 100;
   fRangeX[0] = fMinRangeX;
   fRangeX[1] = fMaxRangeX;

   fStepY = (fMaxRangeY - fMinRangeY) / 100;
   fRangeY[0] = fMinRangeY;
   fRangeY[1] = fMaxRangeY;
}

void ROOT::Experimental::RFitPanelModel::SelectedFunc(const std::string &name, TF1 *func)
{
   fSelectedFunc.clear();
   fFuncPars.Clear();
   if (func) {
      fFuncPars.name = fSelectedFunc = name;
      fFuncPars.GetParameters(func);
   } else {
      fFuncPars.name = "<not exists>";
   }
}


void ROOT::Experimental::RFitPanelModel::UpdateFuncList()
{
   fFuncList.clear();

   if (fDim == 1) {
      fFuncList = { {"gaus"}, {"gausn"}, {"expo"}, {"landau"},{"landaun"},
                    {"pol0"},{"pol1"},{"pol2"},{"pol3"},{"pol4"},{"pol5"},{"pol6"},{"pol7"},{"pol8"},{"pol9"},
                    {"cheb0"}, {"cheb1"}, {"cheb2"}, {"cheb3"}, {"cheb4"}, {"cheb5"}, {"cheb6"}, {"cheb7"}, {"cheb8"}, {"use9"} };
   } else if (fDim == 2) {
      fFuncList = { {"xygaus"}, {"bigaus"}, {"xyexpo"}, {"xylandau"}, {"xylandaun"} };
   }
}


void ROOT::Experimental::RFitPanelModel::Initialize()
{
   // build list of available histograms, as id use name from gdir
   fSelectedData = "";

   // build list of available functions

   // Sub ComboBox for Type Function
   fSelectedFunc = "";
   fDim = 1;
   UpdateFuncList();

   // corresponds when Type == User Func (fSelectedTypeID == 1)

   // ComboBox for General Tab --- Method
   fFitMethods = { {"P", "Chi-square"},
                   {"L", "Log Likelihood"},
                   {"WL", "Binned LogLikelihood"} };
   fFitMethod = "P";

   fLinearFit = false;
   fRobust = false;
   fRobustLevel = 0.95;

   fIntegral = false;
   fAllWeights1 = false;
   fAddToList = false;
   fEmptyBins1 = false;
   fUseGradient = false;

   fSame = false;
   fNoDrawing = false;
   fNoStoreDraw = false;

   // Minimization method
   fLibrary = 0;
   // corresponds to library == 0
   fMethodMinAll = {
         {0, kFP_MIGRAD, "MIGRAD"}, {0, kFP_SIMPLX, "SIMPLEX"}, {0, kFP_SCAN, "SCAN"}, {0, kFP_COMBINATION, "Combination"},
         {1, kFP_MIGRAD, "MIGRAD"}, {1, kFP_SIMPLX, "SIMPLEX"}, {1, kFP_FUMILI2, "FUMILI"}, {1, kFP_SCAN, "SCAN"}, {1, kFP_COMBINATION, "Combination"},
         {2, kFP_FUMILI, "FUMILI"},

         {3, kFP_GSLFR, "Fletcher-Reeves conjugate gradient"},
         {3, kFP_GSLPR, "Polak-Ribiere conjugate gradient"},
         {3, kFP_BFGS,  "BFGS conjugate gradient"},
         {3, kFP_BFGS2, "BFGS conjugate gradient (Version 2)"},
         {3, kFP_GSLLM, "Levenberg-Marquardt"},
         {3, kFP_GSLSA, "Simulated Annealing"}
   };

   fHasGenetics = false;
   if ( gPluginMgr->FindHandler("ROOT::Math::Minimizer","GAlibMin") ) {
      fMethodMinAll.emplace_back(4, kFP_GALIB, "GA Lib Genetic Algorithm");
      fHasGenetics = true;
   }
   if (gPluginMgr->FindHandler("ROOT::Math::Minimizer","Genetic")) {
      fMethodMinAll.emplace_back(4, kFP_TMVAGA, "TMVA Genetic Algorithm");
      fHasGenetics = true;
   }

   fSelectMethodMin = kFP_MIGRAD;

   // fOperation = 0;
   fPrint = 0;
}

/// Update advanced parameters associated with fit function for histogram

void ROOT::Experimental::RFitPanelModel::UpdateAdvanced(TF1 *func)
{
   fContour1.clear();
   fContour2.clear();
   fScan.clear();
   fContourPar1Id = "0";
   fContourPar2Id = "0";
   fScanId = "0";

   fHasAdvanced = (func!=nullptr);

   if (func) {
      for (int n = 0; n < func->GetNpar(); ++n) {
         fContour1.emplace_back(std::to_string(n), func->GetParName(n));
         fContour2.emplace_back(std::to_string(n), func->GetParName(n));
         fScan.emplace_back(std::to_string(n), func->GetParName(n));
      }
      fFuncPars.GetParameters(func); // take func parameters
      fFuncPars.name = "hist::"s + func->GetName(); // clearly mark this as function from histogram
   } else {
      // fFuncPars.Clear();
   }
}


std::string ROOT::Experimental::RFitPanelModel::GetFitOption()
{
   std::string opt = fFitMethod;

   if (fIntegral) opt.append("I");
   if (fUseRange) opt.append("R");
   if (fBestErrors) opt.append("E");
   if (fImproveFitResults) opt.append("M");
   if (fAddToList) opt.append("+");
   if (fUseGradient) opt.append("G");

   if (fEmptyBins1)
      opt.append("WW");
   else if (fAllWeights1)
      opt.append("W");

   if (fNoStoreDraw)
      opt.append("N");
   else if (fNoDrawing)
      opt.append("O");

   return opt;
}

void ROOT::Experimental::RFitPanelModel::GetRanges(ROOT::Fit::DataRange &drange)
{
   if (fDim > 0)  {
      drange.AddRange(0, fRangeX[0], fRangeX[1]);
   }

   if ( fDim > 1 ) {
      drange.AddRange(1, fRangeY[0], fRangeY[1]);
   }

}

void ROOT::Experimental::RFitPanelModel::GetFitOptions(Foption_t &fitOpts)
{
   fitOpts.Range    = fUseRange;
   fitOpts.Integral = fIntegral;
   fitOpts.More     = fImproveFitResults;
   fitOpts.Errors   = fBestErrors;
   fitOpts.Like     = false; // (fMethodList->GetSelected() != kFP_MCHIS);

   if (fEmptyBins1)
      fitOpts.W1 = 2;
   else if (fAllWeights1)
      fitOpts.W1 = 1;

   // TODO: fEnteredFunc->GetText();
   TString tmpStr = ""; // fEnteredFunc->GetText();
   if ( !fLinearFit && (tmpStr.Contains("pol") || tmpStr.Contains("++")) )
      fitOpts.Minuit = 1;

   // TODO: fChangedParams
   bool fChangedParams = false;
   if (fChangedParams) {
      fitOpts.Bound = 1;
      fChangedParams = false;  // reset
   }

   //fitOpts.Nochisq  = (fNoChi2->GetState() == kButtonDown);
   fitOpts.Nostore  = fNoStoreDraw;
   fitOpts.Nograph  = fNoDrawing;
   fitOpts.Plus     = false; // TODO: (fAdd2FuncList->GetState() == kButtonDown);
   fitOpts.Gradient = fUseGradient;
   fitOpts.Quiet    = fPrint == 2;
   fitOpts.Verbose  = fPrint == 1;

   // TODO: only TGraph
   if ( /* !(fType != kObjectGraph) &&  */ fRobust ) {
      fitOpts.Robust = 1;
      fitOpts.hRobust = fRobustLevel;
   }
}

void ROOT::Experimental::RFitPanelModel::GetMinimizerOptions(ROOT::Math::MinimizerOptions &minOpts)
{
   if (fLibrary == 0)
      minOpts.SetMinimizerType ( "Minuit");
   else if (fLibrary == 1)
      minOpts.SetMinimizerType ( "Minuit2" );
   else if (fLibrary == 2)
      minOpts.SetMinimizerType ("Fumili" );
   else if (fLibrary == 3)
      minOpts.SetMinimizerType ("GSLMultiMin" );
   else if (fLibrary == 4)
      minOpts.SetMinimizerType ("Geneti2c" ); // should be handled separately

   switch(fSelectMethodMin) {
      case kFP_MIGRAD:  minOpts.SetMinimizerAlgorithm( "Migrad" ); break;
      case kFP_FUMILI:  minOpts.SetMinimizerAlgorithm( "Fumili" ); break;
      case kFP_FUMILI2: minOpts.SetMinimizerAlgorithm( "Fumili2" ); break;
      case kFP_SIMPLX:  minOpts.SetMinimizerAlgorithm( "Simplex" ); break;
      case kFP_SCAN:    minOpts.SetMinimizerAlgorithm( "Scan" ); break;
      case kFP_COMBINATION: minOpts.SetMinimizerAlgorithm( "Minimize" ); break;
      case kFP_GSLFR:  minOpts.SetMinimizerAlgorithm( "conjugatefr" ); break;
      case kFP_GSLPR:  minOpts.SetMinimizerAlgorithm( "conjugatepr" ); break;
      case kFP_BFGS:   minOpts.SetMinimizerAlgorithm( "bfgs" ); break;
      case kFP_BFGS2:  minOpts.SetMinimizerAlgorithm( "bfgs2" ); break;
      case kFP_GSLLM:
         minOpts.SetMinimizerType ("GSLMultiFit" );
         minOpts.SetMinimizerAlgorithm( "" );
         break;
      case kFP_GSLSA:
         minOpts.SetMinimizerType ("GSLSimAn" );
         minOpts.SetMinimizerAlgorithm( "" );
         break;
      case kFP_TMVAGA:
         minOpts.SetMinimizerType ("Geneti2c" );
         minOpts.SetMinimizerAlgorithm( "" );
         break;
      case kFP_GALIB:
         minOpts.SetMinimizerType ("GAlibMin" );
         minOpts.SetMinimizerAlgorithm( "" );
         break;
      default:
         minOpts.SetMinimizerAlgorithm( "" );
         break;
   }

   minOpts.SetErrorDef (fErrorDef);
   minOpts.SetTolerance(fMaxTolerance);
   minOpts.SetMaxIterations(fMaxIterations);
   minOpts.SetMaxFunctionCalls(fMaxIterations);
}

std::string ROOT::Experimental::RFitPanelModel::GetDrawOption()
{
   return "";
}






