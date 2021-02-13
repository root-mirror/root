// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*! \class TMVA::MethodCrossValidation
\ingroup TMVA
*/
#include "TMVA/MethodCrossValidation.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/CvSplit.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TSystem.h"

REGISTER_METHOD(CrossValidation)

ClassImp(TMVA::MethodCrossValidation);

////////////////////////////////////////////////////////////////////////////////
///

TMVA::MethodCrossValidation::MethodCrossValidation(const TString &jobName, const TString &methodTitle,
                                                   DataSetInfo &theData, const TString &theOption)
   : TMVA::MethodBase(jobName, Types::kCrossValidation, methodTitle, theData, theOption), fSplitExpr(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////

TMVA::MethodCrossValidation::MethodCrossValidation(DataSetInfo &theData, const TString &theWeightFile)
   : TMVA::MethodBase(Types::kCrossValidation, theData, theWeightFile), fSplitExpr(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
///

TMVA::MethodCrossValidation::~MethodCrossValidation(void) {}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodCrossValidation::DeclareOptions()
{
   DeclareOptionRef(fEncapsulatedMethodName, "EncapsulatedMethodName", "");
   DeclareOptionRef(fEncapsulatedMethodTypeName, "EncapsulatedMethodTypeName", "");
   DeclareOptionRef(fNumFolds, "NumFolds", "Number of folds to generate");
   DeclareOptionRef(fOutputEnsembling = TString("None"), "OutputEnsembling",
                    "Combines output from contained methods. If None, no combination is performed. (default None)");
   AddPreDefVal(TString("None"));
   AddPreDefVal(TString("Avg"));
   DeclareOptionRef(fSplitExprString, "SplitExpr", "The expression used to assign events to folds");
}

////////////////////////////////////////////////////////////////////////////////
/// Options that are used ONLY for the READER to ensure backward compatibility.

void TMVA::MethodCrossValidation::DeclareCompatibilityOptions()
{
   MethodBase::DeclareCompatibilityOptions();
}

////////////////////////////////////////////////////////////////////////////////
/// The option string is decoded, for available options see "DeclareOptions".

void TMVA::MethodCrossValidation::ProcessOptions()
{
   Log() << kDEBUG << "ProcessOptions -- fNumFolds: " << fNumFolds << Endl;
   Log() << kDEBUG << "ProcessOptions -- fEncapsulatedMethodName: " << fEncapsulatedMethodName << Endl;
   Log() << kDEBUG << "ProcessOptions -- fEncapsulatedMethodTypeName: " << fEncapsulatedMethodTypeName << Endl;

   if (fSplitExprString != TString("")) {
      fSplitExpr = std::unique_ptr<CvSplitKFoldsExpr>(new CvSplitKFoldsExpr(DataInfo(), fSplitExprString));
   }

   for (UInt_t iFold = 0; iFold < fNumFolds; ++iFold) {
      TString weightfile = GetWeightFileNameForFold(iFold);

      Log() << kINFO << "Reading weightfile: " << weightfile << Endl;
      MethodBase *fold_method = InstantiateMethodFromXML(fEncapsulatedMethodTypeName, weightfile);
      fEncapsulatedMethods.push_back(fold_method);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialisation with defaults for the Method.

void TMVA::MethodCrossValidation::Init(void)
{
   fMulticlassValues = std::vector<Float_t>(DataInfo().GetNClasses());
   fRegressionValues = std::vector<Float_t>(DataInfo().GetNTargets());
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the method, as if it had just been instantiated (forget all training etc.).

void TMVA::MethodCrossValidation::Reset(void) {}

////////////////////////////////////////////////////////////////////////////////
/// \brief Returns filename of weight file for a given fold.
/// \param[in] iFold Ordinal of the fold. Range: 0 to NumFolds exclusive.
///
TString TMVA::MethodCrossValidation::GetWeightFileNameForFold(UInt_t iFold) const
{
   if (iFold >= fNumFolds) {
      Log() << kFATAL << iFold << " out of range. "
            << "Should be < " << fNumFolds << "." << Endl;
   }

   TString foldStr = Form("fold%i", iFold + 1);
   TString fileDir = gSystem->GetDirName(GetWeightFileName());
   TString weightfile = fileDir + "/" + fJobName + "_" + fEncapsulatedMethodName + "_" + foldStr + ".weights.xml";

   return weightfile;
}

////////////////////////////////////////////////////////////////////////////////
/// Call the Optimizer with the set of parameters and ranges that
/// are meant to be tuned.

// std::map<TString,Double_t>  TMVA::MethodCrossValidation::OptimizeTuningParameters(TString fomType, TString fitType)
// {
// }

////////////////////////////////////////////////////////////////////////////////
/// Set the tuning parameters according to the argument.

// void TMVA::MethodCrossValidation::SetTuneParameters(std::map<TString,Double_t> tuneParameters)
// {
// }

////////////////////////////////////////////////////////////////////////////////
///  training.

void TMVA::MethodCrossValidation::Train() {}

////////////////////////////////////////////////////////////////////////////////
/// \brief Reads in a weight file an instantiates the corresponding method
/// \param[in] methodTypeName Canonical name of the method type. E.g. `"BDT"`
///                           for Boosted Decision Trees.
/// \param[in] weightfile File to read method parameters from
TMVA::MethodBase *
TMVA::MethodCrossValidation::InstantiateMethodFromXML(TString methodTypeName, TString weightfile) const
{
   TMVA::MethodBase *m = dynamic_cast<MethodBase *>(
      ClassifierFactory::Instance().Create(std::string(methodTypeName.Data()), DataInfo(), weightfile));

   if (m->GetMethodType() == Types::kCategory) {
      Log() << kFATAL << "MethodCategory not supported for the moment." << Endl;
   }

   TString fileDir = TString(DataInfo().GetName()) + "/" + gConfig().GetIONames().fWeightFileDir;
   m->SetWeightFileDir(fileDir);
   // m->SetModelPersistence(fModelPersistence);
   // m->SetSilentFile(IsSilentFile());
   m->SetAnalysisType(fAnalysisType);
   m->SetupMethod();
   m->ReadStateFromFile();
   // m->SetTestvarName(testvarName);

   return m;
}

////////////////////////////////////////////////////////////////////////////////
/// Write weights to XML.

void TMVA::MethodCrossValidation::AddWeightsXMLTo(void *parent) const
{
   void *wght = gTools().AddChild(parent, "Weights");

   gTools().AddAttr(wght, "JobName", fJobName);
   gTools().AddAttr(wght, "SplitExpr", fSplitExprString);
   gTools().AddAttr(wght, "NumFolds", fNumFolds);
   gTools().AddAttr(wght, "EncapsulatedMethodName", fEncapsulatedMethodName);
   gTools().AddAttr(wght, "EncapsulatedMethodTypeName", fEncapsulatedMethodTypeName);
   gTools().AddAttr(wght, "OutputEnsembling", fOutputEnsembling);

   for (UInt_t iFold = 0; iFold < fNumFolds; ++iFold) {
      TString weightfile = GetWeightFileNameForFold(iFold);

      // TODO: Add a swithch in options for using either split files or only one.
      // TODO: This would store the method inside MethodCrossValidation
      //       Another option is to store the folds as separate files.
      // // Retrieve encap. method for fold n
      // MethodBase * method = InstantiateMethodFromXML(fEncapsulatedMethodTypeName, weightfile);
      //
      // // Serialise encapsulated method for fold n
      // void* foldNode = gTools().AddChild(parent, foldStr);
      // method->WriteStateToXML(foldNode);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reads from the xml file.
///

void TMVA::MethodCrossValidation::ReadWeightsFromXML(void *parent)
{
   gTools().ReadAttr(parent, "JobName", fJobName);
   gTools().ReadAttr(parent, "SplitExpr", fSplitExprString);
   gTools().ReadAttr(parent, "NumFolds", fNumFolds);
   gTools().ReadAttr(parent, "EncapsulatedMethodName", fEncapsulatedMethodName);
   gTools().ReadAttr(parent, "EncapsulatedMethodTypeName", fEncapsulatedMethodTypeName);
   gTools().ReadAttr(parent, "OutputEnsembling", fOutputEnsembling);

   // Read in methods for all folds
   for (UInt_t iFold = 0; iFold < fNumFolds; ++iFold) {
      TString weightfile = GetWeightFileNameForFold(iFold);

      Log() << kINFO << "Reading weightfile: " << weightfile << Endl;
      MethodBase *fold_method = InstantiateMethodFromXML(fEncapsulatedMethodTypeName, weightfile);
      fEncapsulatedMethods.push_back(fold_method);
   }

   // SplitExpr
   if (fSplitExprString != TString("")) {
      fSplitExpr = std::unique_ptr<CvSplitKFoldsExpr>(new CvSplitKFoldsExpr(DataInfo(), fSplitExprString));
   } else {
      Log() << kFATAL << "MethodCrossValidation supports XML reading only for deterministic splitting !" << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read the weights
///

void TMVA::MethodCrossValidation::ReadWeightsFromStream(std::istream & /*istr*/)
{
   Log() << kFATAL << "CrossValidation currently supports only reading from XML." << Endl;
}

////////////////////////////////////////////////////////////////////////////////
///

Double_t TMVA::MethodCrossValidation::GetMvaValue(Double_t *err, Double_t *errUpper)
{
   const Event *ev = GetEvent();

   if (fOutputEnsembling == "None") {
      if (fSplitExpr != nullptr) {
         // K-folds with a deterministic split
         UInt_t iFold = fSplitExpr->Eval(fNumFolds, ev);
         return fEncapsulatedMethods.at(iFold)->GetMvaValue(err, errUpper);
      } else {
         // K-folds with a random split was used
         UInt_t iFold = fEventToFoldMapping.at(Data()->GetEvent());
         return fEncapsulatedMethods.at(iFold)->GetMvaValue(err, errUpper);
      }
   } else if (fOutputEnsembling == "Avg") {
      Double_t val = 0.0;
      for (auto &m : fEncapsulatedMethods) {
         val += m->GetMvaValue(err, errUpper);
      }
      return val / fEncapsulatedMethods.size();
   } else {
      Log() << kFATAL << "Ensembling type " << fOutputEnsembling << " unknown" << Endl;
      return 0; // Cannot happen
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the multiclass MVA response.

const std::vector<Float_t> &TMVA::MethodCrossValidation::GetMulticlassValues()
{
   const Event *ev = GetEvent();

   if (fOutputEnsembling == "None") {
      if (fSplitExpr != nullptr) {
         // K-folds with a deterministic split
         UInt_t iFold = fSplitExpr->Eval(fNumFolds, ev);
         return fEncapsulatedMethods.at(iFold)->GetMulticlassValues();
      } else {
         // K-folds with a random split was used
         UInt_t iFold = fEventToFoldMapping.at(Data()->GetEvent());
         return fEncapsulatedMethods.at(iFold)->GetMulticlassValues();
      }
   } else if (fOutputEnsembling == "Avg") {

      for (auto &e : fMulticlassValues) {
         e = 0;
      }

      for (auto &m : fEncapsulatedMethods) {
         auto methodValues = m->GetMulticlassValues();
         for (size_t i = 0; i < methodValues.size(); ++i) {
            fMulticlassValues[i] += methodValues[i];
         }
      }

      for (auto &e : fMulticlassValues) {
         e /= fEncapsulatedMethods.size();
      }

      return fMulticlassValues;

   } else {
      Log() << kFATAL << "Ensembling type " << fOutputEnsembling << " unknown" << Endl;
      return fMulticlassValues; // Cannot happen
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the regression value generated by the containing methods.

const std::vector<Float_t> &TMVA::MethodCrossValidation::GetRegressionValues()
{
   const Event *ev = GetEvent();

   if (fOutputEnsembling == "None") {
      if (fSplitExpr != nullptr) {
         // K-folds with a deterministic split
         UInt_t iFold = fSplitExpr->Eval(fNumFolds, ev);
         return fEncapsulatedMethods.at(iFold)->GetRegressionValues();
      } else {
         // K-folds with a random split was used
         UInt_t iFold = fEventToFoldMapping.at(Data()->GetEvent());
         return fEncapsulatedMethods.at(iFold)->GetRegressionValues();
      }
   } else if (fOutputEnsembling == "Avg") {

      for (auto &e : fRegressionValues) {
         e = 0;
      }

      for (auto &m : fEncapsulatedMethods) {
         auto methodValues = m->GetRegressionValues();
         for (size_t i = 0; i < methodValues.size(); ++i) {
            fRegressionValues[i] += methodValues[i];
         }
      }

      for (auto &e : fRegressionValues) {
         e /= fEncapsulatedMethods.size();
      }

      return fRegressionValues;

   } else {
      Log() << kFATAL << "Ensembling type " << fOutputEnsembling << " unknown" << Endl;
      return fRegressionValues; // Cannot happen
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::MethodCrossValidation::WriteMonitoringHistosToFile(void) const
{
   // // Used for evaluation, which is outside the life time of MethodCrossEval.
   // Log() << kFATAL << "Method CrossValidation should not be created manually,"
   //                    " only as part of using TMVA::Reader." << Endl;
   // return;
}

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::MethodCrossValidation::GetHelpMessage() const
{
   Log() << kWARNING
         << "Method CrossValidation should not be created manually,"
            " only as part of using TMVA::Reader."
         << Endl;
}

////////////////////////////////////////////////////////////////////////////////
///

const TMVA::Ranking *TMVA::MethodCrossValidation::CreateRanking()
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::MethodCrossValidation::HasAnalysisType(Types::EAnalysisType /*type*/, UInt_t /*numberClasses*/,
                                                    UInt_t /*numberTargets*/)
{
   return kTRUE;
   // if (fEncapsulatedMethods.size() == 0) {return kFALSE;}
   // if (fEncapsulatedMethods.at(0) == nullptr) {return kFALSE;}
   // return fEncapsulatedMethods.at(0)->HasAnalysisType(type, numberClasses, numberTargets);
}

////////////////////////////////////////////////////////////////////////////////
/// Make ROOT-independent C++ class for classifier response (classifier-specific implementation).

void TMVA::MethodCrossValidation::MakeClassSpecific(std::ostream & /*fout*/, const TString & /*className*/) const
{
   Log() << kWARNING << "MakeClassSpecific not implemented for CrossValidation" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Specific class header.

void TMVA::MethodCrossValidation::MakeClassSpecificHeader(std::ostream & /*fout*/, const TString & /*className*/) const
{
   Log() << kWARNING << "MakeClassSpecificHeader not implemented for CrossValidation" << Endl;
}
