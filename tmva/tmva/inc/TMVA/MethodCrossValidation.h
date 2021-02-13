// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMVA_MethodCrossValidation
#define ROOT_TMVA_MethodCrossValidation

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodCrossValidation                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/CvSplit.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MethodBase.h"

#include "TString.h"

#include <iostream>
#include <memory>
#include <vector>
#include <map>

namespace TMVA {

class CrossValidation;
class Ranking;

// Looks for serialised methods of the form methodTitle + "_fold" + iFold;
class MethodCrossValidation : public MethodBase {

   friend CrossValidation;

public:
   // constructor for training and reading
   MethodCrossValidation(const TString &jobName, const TString &methodTitle, DataSetInfo &theData,
                         const TString &theOption = "");

   // constructor for calculating BDT-MVA using previously generatad decision trees
   MethodCrossValidation(DataSetInfo &theData, const TString &theWeightFile);

   virtual ~MethodCrossValidation(void);

   // optimize tuning parameters
   // virtual std::map<TString,Double_t> OptimizeTuningParameters(TString fomType="ROCIntegral", TString
   // fitType="FitGA"); virtual void SetTuneParameters(std::map<TString,Double_t> tuneParameters);

   // training method
   void Train(void);

   // revoke training
   void Reset(void);

   using MethodBase::ReadWeightsFromStream;

   // write weights to file
   void AddWeightsXMLTo(void *parent) const;

   // read weights from file
   void ReadWeightsFromStream(std::istream &istr);
   void ReadWeightsFromXML(void *parent);

   // write method specific histos to target file
   void WriteMonitoringHistosToFile(void) const;

   // calculate the MVA value
   Double_t GetMvaValue(Double_t *err = 0, Double_t *errUpper = 0);
   const std::vector<Float_t> &GetMulticlassValues();
   const std::vector<Float_t> &GetRegressionValues();

   // the option handling methods
   void DeclareOptions();
   void ProcessOptions();

   // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
   void MakeClassSpecific(std::ostream &, const TString &) const;
   void MakeClassSpecificHeader(std::ostream &, const TString &) const;

   void GetHelpMessage() const;

   const Ranking *CreateRanking();
   Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

protected:
   void Init(void);
   void DeclareCompatibilityOptions();

private:
   TString GetWeightFileNameForFold(UInt_t iFold) const;
   MethodBase *InstantiateMethodFromXML(TString methodTypeName, TString weightfile) const;

private:
   TString fEncapsulatedMethodName;
   TString fEncapsulatedMethodTypeName;
   UInt_t fNumFolds;
   TString fOutputEnsembling;

   TString fSplitExprString;
   std::unique_ptr<CvSplitKFoldsExpr> fSplitExpr;

   std::vector<Float_t> fMulticlassValues;
   std::vector<Float_t> fRegressionValues;

   std::vector<MethodBase *> fEncapsulatedMethods;

   // Used for CrossValidation with random splits (not using the
   // CVSplitCrossValisationExpr functionality) to communicate Event to fold
   // mapping.
   std::map<const TMVA::Event *, UInt_t> fEventToFoldMapping;

   // for backward compatibility
   ClassDef(MethodCrossValidation, 0);
};

} // namespace TMVA

#endif
