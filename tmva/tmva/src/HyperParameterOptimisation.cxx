// @(#)root/tmva $Id$
// Author: Omar Zapata, Thomas James Stevenson.

#include "TMVA/HyperParameterOptimisation.h"

#include "TMVA/Configurable.h"
#include "TMVA/CvSplit.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/Types.h"

#include "TMultiGraph.h"
#include "TString.h"

#include <memory>
#include <vector>

/*! \class TMVA::HyperParameterOptimisationResult
\ingroup TMVA

*/

/*! \class TMVA::HyperParameterOptimisation
\ingroup TMVA

*/

//_______________________________________________________________________
TMVA::HyperParameterOptimisationResult::HyperParameterOptimisationResult()
   : fROCAVG(0.0), fROCCurves(std::make_shared<TMultiGraph>())
{
}

//_______________________________________________________________________
TMVA::HyperParameterOptimisationResult::~HyperParameterOptimisationResult()
{
}

//_______________________________________________________________________
TMultiGraph *TMVA::HyperParameterOptimisationResult::GetROCCurves(Bool_t /* fLegend */)
{

    return fROCCurves.get();
}

//_______________________________________________________________________
void TMVA::HyperParameterOptimisationResult::Print() const
{
    TMVA::MsgLogger::EnableOutput();
    TMVA::gConfig().SetSilent(kFALSE);

    MsgLogger fLogger("HyperParameterOptimisation");

    for(UInt_t j=0; j<fFoldParameters.size(); ++j) {
        fLogger<<kHEADER<< "===========================================================" << Endl;
        fLogger<<kINFO<< "Optimisation for " << fMethodName << " fold " << j+1 << Endl;

        for(auto &it : fFoldParameters.at(j)) {
            fLogger<<kINFO<< it.first << "     " << it.second << Endl;
        }
    }

    TMVA::gConfig().SetSilent(kTRUE);

}

//_______________________________________________________________________
TMVA::HyperParameterOptimisation::HyperParameterOptimisation(TMVA::DataLoader *dataloader):Envelope("HyperParameterOptimisation",dataloader),
    fFomType("Separation"),
    fFitType("Minuit"),
    fNumFolds(5),
    fResults(),
    fClassifier(new TMVA::Factory("HyperParameterOptimisation","!V:!ROC:Silent:!ModelPersistence:!Color:!DrawProgressBar:AnalysisType=Classification"))
{
    fFoldStatus=kFALSE;
}

//_______________________________________________________________________
TMVA::HyperParameterOptimisation::~HyperParameterOptimisation()
{
    fClassifier=nullptr;
}

//_______________________________________________________________________
void TMVA::HyperParameterOptimisation::SetNumFolds(UInt_t i)
{
   fNumFolds = i;
   // fDataLoader->MakeKFoldDataSet(fNumFolds);
   fFoldStatus = kFALSE;
}

//_______________________________________________________________________
void TMVA::HyperParameterOptimisation::Evaluate()
{
   for (auto &meth : fMethods) {
      TString methodName = meth.GetValue<TString>("MethodName");
      TString methodTitle = meth.GetValue<TString>("MethodTitle");
      TString methodOptions = meth.GetValue<TString>("MethodOptions");

      CvSplitKFolds split{fNumFolds, "", kFALSE, 0};
      if (!fFoldStatus) {
         fDataLoader->MakeKFoldDataSet(split);
         fFoldStatus = kTRUE;
      }
      fResults.fMethodName = methodName;

      for (UInt_t i = 0; i < fNumFolds; ++i) {
         TString foldTitle = methodTitle;
         foldTitle += "_opt";
         foldTitle += i + 1;

         Event::SetIsTraining(kTRUE);
         fDataLoader->PrepareFoldDataSet(split, i, TMVA::Types::kTraining);

         auto smethod = fClassifier->BookMethod(fDataLoader.get(), methodName, methodTitle, methodOptions);

         auto params = smethod->OptimizeTuningParameters(fFomType, fFitType);
         fResults.fFoldParameters.push_back(params);

         smethod->Data()->DeleteResults(smethod->GetMethodName(), Types::kTraining, Types::kClassification);

         fClassifier->DeleteAllMethods();

         fClassifier->fMethodsMap.clear();
      }
   }
}
