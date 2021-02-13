// Authors: Sergey Linev <S.Linev@gsi.de> Iliana Betsou <Iliana.Betsou@cern.ch>
// Date: 2019-04-11
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RFitPanelModel.hxx>

#include <ROOT/RLogger.hxx>

#include "TH1.h"
#include "TPluginManager.h"
#include "TFitResult.h"
#include "TF1.h"

ROOT::Experimental::RLogChannel &ROOT::Experimental::FitPanelLog() {
   static RLogChannel sLog("ROOT.FitPanel");
   return sLog;
}


enum EFitPanel {
   kFP_NONE = 0,
   kFP_MIGRAD, kFP_SIMPLX, kFP_SCAN, kFP_COMBINATION,
   kFP_FUMILI, kFP_FUMILI2, kFP_GSLFR, kFP_GSLPR,
   kFP_BFGS, kFP_BFGS2, kFP_GSLLM, kFP_GSLSA,
   kFP_GALIB, kFP_TMVAGA,

   kFP_MCHIS, kFP_MBINL, kFP_MUBIN // Fit methods

};

using namespace std::string_literals;


/////////////////////////////////////////////////////////
/// Clear list of patameters

void ROOT::Experimental::RFitPanelModel::RFuncParsList::Clear()
{
   haspars = false;
   id.clear();
   name.clear();
   pars.clear();
}

/////////////////////////////////////////////////////////
/// Extract TF1 parameters, used in editor on client sides

void ROOT::Experimental::RFitPanelModel::RFuncParsList::GetParameters(TF1 *func)
{
   pars.clear();
   haspars = true;
   name = func->GetName();

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

/////////////////////////////////////////////////////////
/// Set parameters to TF1

void ROOT::Experimental::RFitPanelModel::RFuncParsList::SetParameters(TF1 *func)
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

/////////////////////////////////////////////////////////
/// Update range values

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

/////////////////////////////////////////////////////////
/// Check if function id is exists

bool ROOT::Experimental::RFitPanelModel::HasFunction(const std::string &id)
{
   if (id.empty())
      return false;

   for (auto &item : fFuncList)
      if (item.id == id)
         return true;

   return false;
}

/////////////////////////////////////////////////////////
/// Select function

void ROOT::Experimental::RFitPanelModel::SelectedFunc(const std::string &id, TF1 *func)
{
   if (HasFunction(id))
      fSelectedFunc = id;
   else
      fSelectedFunc.clear();

   fFuncPars.Clear();
   if (func) {
      fFuncPars.id = id;
      fFuncPars.GetParameters(func);
   }
}

/////////////////////////////////////////////////////////
/// Initialize model - set some meaningful default values

void ROOT::Experimental::RFitPanelModel::Initialize()
{
   // build list of available histograms, as id use name from gdir
   fSelectedData = "";

   // build list of available functions

   // Sub ComboBox for Type Function
   fSelectedFunc = "";
   fDim = 1;

   fSelectedTab = "General";

   // corresponds when Type == User Func (fSelectedTypeID == 1)

   // ComboBox for General Tab --- Method
   fFitMethods.clear();
   fFitMethod = 0;

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

   fPrint = 0;

   fAdvancedTab = "Contour";

   fContourPoints = 40;
   fScanPoints = 40;
   fConfidenceLevel = 0.683;
   fContourColor = "#c0b6ac"; // #11
   fScanColor = "#0000ff"; // kBlue
   fConfidenceColor = "#c0b6ac"; // #11
}


////////////////////////////////////////////////////////////
/// Update setting dependent from object type

void ROOT::Experimental::RFitPanelModel::SetObjectKind(EFitObjectType kind)
{
   fDataType = kind;

   fFitMethods.clear();
   fFitMethod = 0;
   fRobust = false;

   switch (kind) {
      case kObjectHisto:
         fFitMethods = {{kFP_MCHIS, "Chi-square"}, {kFP_MBINL, "Binned Likelihood"}};
         fFitMethod = kFP_MCHIS;
         break;

      case kObjectGraph:
         fFitMethods = {{kFP_MCHIS, "Chi-square"}};
         fFitMethod = kFP_MCHIS;
         fRobust = true;
         break;

      case kObjectMultiGraph:
         fFitMethods = {{kFP_MCHIS, "Chi-square"}};
         fFitMethod = kFP_MCHIS;
         fRobust = true;
         break;

      case kObjectGraph2D:
         fFitMethods = {{kFP_MCHIS, "Chi-square"}};
         fFitMethod = kFP_MCHIS;
         break;

      case kObjectHStack:
         fFitMethods = {{kFP_MCHIS, "Chi-square"}};
         fFitMethod = kFP_MCHIS;
         break;

      //case kObjectTree:
      //   fFitMethods = {{ kFP_MUBIN, "Unbinned Likelihood" }};
      //   fFitMethod = kFP_MUBIN;
      //   break;

      default:
         break;
   }

}


////////////////////////////////////////////////////////////
/// Update advanced parameters associated with fit function

void ROOT::Experimental::RFitPanelModel::UpdateAdvanced(TFitResult *res)
{
   fAdvancedPars.clear();

   fHasAdvanced = (res!=nullptr);

   auto checkid = [&](std::string &id, const std::string &dflt) {
      if (!res) { id.clear(); return; }
      for (auto &item : fAdvancedPars)
         if (item.key == id) return;
      id = dflt;
   };

   if (res)
      for (unsigned n = 0; n < res->NPar(); ++n)
         fAdvancedPars.emplace_back(std::to_string(n), res->ParName(n));

   checkid(fContourPar1Id, "0");
   checkid(fContourPar2Id, "1");
   checkid(fScanId, "0");
}


ROOT::Fit::DataRange ROOT::Experimental::RFitPanelModel::GetRanges()
{
   ROOT::Fit::DataRange drange;

   if (fDim > 0)
      drange.AddRange(0, fRangeX[0], fRangeX[1]);

   if ( fDim > 1 )
      drange.AddRange(1, fRangeY[0], fRangeY[1]);

   return drange;
}

/////////////////////////////////////////////////////////
/// Provide initialized Foption_t instance

Foption_t ROOT::Experimental::RFitPanelModel::GetFitOptions()
{
   Foption_t fitOpts;
   fitOpts.Range    = fUseRange;
   fitOpts.Integral = fIntegral;
   fitOpts.More     = fImproveFitResults;
   fitOpts.Errors   = fBestErrors;
   fitOpts.Like     = fFitMethod != kFP_MCHIS;

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

   if (fRobust) {
      fitOpts.Robust = 1;
      fitOpts.hRobust = fRobustLevel;
   }

   return fitOpts;
}

/////////////////////////////////////////////////////////
/// Provide initialized ROOT::Math::MinimizerOptions instance

ROOT::Math::MinimizerOptions ROOT::Experimental::RFitPanelModel::GetMinimizerOptions()
{
   ROOT::Math::MinimizerOptions minOpts;

   switch (fLibrary) {
      case 0: minOpts.SetMinimizerType("Minuit"); break;
      case 1: minOpts.SetMinimizerType("Minuit2"); break;
      case 2: minOpts.SetMinimizerType("Fumili"); break;
      case 3: minOpts.SetMinimizerType("GSLMultiMin"); break;
      case 4: minOpts.SetMinimizerType("Geneti2c"); break;
   }

   switch(fSelectMethodMin) {
      case kFP_MIGRAD:  minOpts.SetMinimizerAlgorithm( "Migrad" ); break;
      case kFP_FUMILI: minOpts.SetMinimizerAlgorithm("Fumili"); break;
      case kFP_FUMILI2: minOpts.SetMinimizerAlgorithm("Fumili2"); break;
      case kFP_SIMPLX: minOpts.SetMinimizerAlgorithm("Simplex"); break;
      case kFP_SCAN: minOpts.SetMinimizerAlgorithm("Scan"); break;
      case kFP_COMBINATION: minOpts.SetMinimizerAlgorithm("Minimize"); break;
      case kFP_GSLFR: minOpts.SetMinimizerAlgorithm("conjugatefr"); break;
      case kFP_GSLPR: minOpts.SetMinimizerAlgorithm("conjugatepr"); break;
      case kFP_BFGS: minOpts.SetMinimizerAlgorithm("bfgs"); break;
      case kFP_BFGS2: minOpts.SetMinimizerAlgorithm("bfgs2"); break;

      case kFP_GSLLM:
         minOpts.SetMinimizerType("GSLMultiFit");
         minOpts.SetMinimizerAlgorithm("");
         break;
      case kFP_GSLSA:
         minOpts.SetMinimizerType("GSLSimAn");
         minOpts.SetMinimizerAlgorithm("");
         break;
      case kFP_TMVAGA:
         minOpts.SetMinimizerType("Geneti2c");
         minOpts.SetMinimizerAlgorithm("");
         break;
      case kFP_GALIB:
         minOpts.SetMinimizerType("GAlibMin");
         minOpts.SetMinimizerAlgorithm("");
         break;
      default: minOpts.SetMinimizerAlgorithm(""); break;
   }

   minOpts.SetErrorDef(fErrorDef);
   minOpts.SetTolerance(fMaxTolerance);
   minOpts.SetMaxIterations(fMaxIterations);
   minOpts.SetMaxFunctionCalls(fMaxIterations);

   return minOpts;
}

/////////////////////////////////////////////////////////
/// Retrun draw option - dummy now

TString ROOT::Experimental::RFitPanelModel::GetDrawOption()
{
   TString res;
   return res;
}
