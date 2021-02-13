// @(#)root/tmva $Id$
// Author: Omar Zapata and Sergei Gleyzer

/*! \class TMVA::VariableImportanceResult
\ingroup TMVA
*/

/*! \class TMVA::VariableImportance
\ingroup TMVA
*/

#include "TMVA/VariableImportance.h"

#include "TMVA/Config.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Envelope.h"
#include "TMVA/Factory.h"
#include "TMVA/OptionMap.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"
#include "TMVA/VarTransformHandler.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TStyle.h"

#include <bitset>
#include <memory>
#include <utility>


//number of bits for bitset
#define NBITS          32

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableImportanceResult::VariableImportanceResult():fImportanceValues("VariableImportance"),
                                                           fImportanceHist(nullptr)
{

}

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableImportanceResult::VariableImportanceResult(const VariableImportanceResult &obj)
{
    fImportanceValues = obj.fImportanceValues;
    fImportanceHist   = obj.fImportanceHist;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::VariableImportanceResult::Print() const
{
    TMVA::MsgLogger::EnableOutput();
    TMVA::gConfig().SetSilent(kFALSE);

    MsgLogger fLogger("VariableImportance");
    if(fType==VIType::kShort)
    {
        fLogger<<kINFO<<"Variable Importance Results (Short)"<<Endl;
    }else if(fType==VIType::kAll)
    {
        fLogger<<kINFO<<"Variable Importance Results (All)"<<Endl;
    }else{
        fLogger<<kINFO<<"Variable Importance Results (Random)"<<Endl;
    }

    fImportanceValues.Print();
    TMVA::gConfig().SetSilent(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////

TCanvas* TMVA::VariableImportanceResult::Draw(const TString name) const
{
    TCanvas *c=new TCanvas(name.Data());
    fImportanceHist->Draw("");
    fImportanceHist->GetXaxis()->SetTitle(" Variable Names ");
    fImportanceHist->GetYaxis()->SetTitle(" Importance (%) ");
    c->Draw();
    return c;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableImportance::VariableImportance(TMVA::DataLoader *dataloader):TMVA::Envelope("VariableImportance",dataloader,nullptr),fType(VIType::kShort)
{
    fClassifier=std::unique_ptr<Factory>(new TMVA::Factory("VariableImportance","!V:!ROC:!ModelPersistence:Silent:Color:!DrawProgressBar:AnalysisType=Classification"));
}

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableImportance::~VariableImportance()
{
    fClassifier=nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::VariableImportance::Evaluate()
{

    //NOTE: Put the type of VI Algorithm in the results Print
    if(fType==VIType::kShort)
    {
        EvaluateImportanceShort();
    }else if(fType==VIType::kAll)
    {
        EvaluateImportanceAll();
    }else{
        UInt_t nbits=fDataLoader->GetDefaultDataSetInfo().GetNVariables();
        if(nbits<10)
            Log()<<kERROR<<"Running variable importance with less that 10 varibales in Random mode "<<
                           "can to produce inconsisten results"<<Endl;
        EvaluateImportanceRandom(pow(nbits,2));
    }
    fResults.fType = fType;
    TMVA::MsgLogger::EnableOutput();
    TMVA::gConfig().SetSilent(kFALSE);
    Log()<<kINFO<<"Evaluation done."<<Endl;
    TMVA::gConfig().SetSilent(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////

ULong_t TMVA::VariableImportance::Sum(ULong_t i)
{
    ULong_t sum=0;
    for(ULong_t n=0;n<i;n++) sum+=pow(2,n);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

TH1F* TMVA::VariableImportance::GetImportance(const UInt_t nbits,std::vector<Float_t> &importances,std::vector<TString> &varNames)
{
    TH1F *vihist  = new TH1F("vihist", "", nbits, 0, nbits);

    gStyle->SetOptStat(000000);

    Float_t normalization = 0.0;
    for (UInt_t i = 0; i < nbits; i++) normalization += importances[i];

    Float_t roc = 0.0;

    gStyle->SetTitleXOffset(0.4);
    gStyle->SetTitleXOffset(1.2);


    for (UInt_t i = 1; i < nbits + 1; i++) {
        roc = 100.0 * importances[i - 1] / normalization;
        vihist->GetXaxis()->SetBinLabel(i, varNames[i - 1].Data());
        vihist->SetBinContent(i, roc);
    }

    vihist->LabelsOption("v >", "X");
    vihist->SetBarWidth(0.97);
    vihist->SetFillColor(TColor::GetColor("#006600"));

    vihist->GetXaxis()->SetTitle(" Variable Names ");
    vihist->GetXaxis()->SetTitleSize(0.045);
    vihist->GetXaxis()->CenterTitle();
    vihist->GetXaxis()->SetTitleOffset(1.24);

    vihist->GetYaxis()->SetTitle(" Importance (%)");
    vihist->GetYaxis()->SetTitleSize(0.045);
    vihist->GetYaxis()->CenterTitle();
    vihist->GetYaxis()->SetTitleOffset(1.24);

    vihist->GetYaxis()->SetRangeUser(-7, 50);
    vihist->SetDirectory(0);

    return vihist;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::VariableImportance::EvaluateImportanceShort()
{
   for (auto &meth : fMethods) {
      TString methodName = meth.GetValue<TString>("MethodName");
      TString methodTitle = meth.GetValue<TString>("MethodTitle");
      TString methodOptions = meth.GetValue<TString>("MethodOptions");

      uint32_t x = 0;
      uint32_t y = 0;
      // getting number of variables and variable names from loader
      const UInt_t nbits = fDataLoader->GetDefaultDataSetInfo().GetNVariables();
      std::vector<TString> varNames = fDataLoader->GetDefaultDataSetInfo().GetListOfVariables();

      ULong_t range = Sum(nbits);

      // vector to save importances
      std::vector<Float_t> importances(nbits);
      for (UInt_t i = 0; i < nbits; i++)
         importances[i] = 0;

      Float_t SROC, SSROC; // computed ROC value for every Seed and SubSeed

      x = range;

      std::bitset<NBITS> xbitset(x);
      if (x == 0)
         Log() << kFATAL << "Error: need at least one variable."; // dataloader need at least one variable

      // creating loader for seed
      TMVA::DataLoader *seeddl = new TMVA::DataLoader(xbitset.to_string());

      // adding variables from seed
      for (UInt_t index = 0; index < nbits; index++) {
         if (xbitset[index])
            seeddl->AddVariable(varNames[index], 'F');
      }

      // Loading Dataset
      DataLoaderCopy(seeddl, fDataLoader.get());

      // Booking Seed
      fClassifier->BookMethod(seeddl, methodName, methodTitle, methodOptions);

      // Train/Test/Evaluation
      fClassifier->TrainAllMethods();
      fClassifier->TestAllMethods();
      fClassifier->EvaluateAllMethods();

      // getting ROC
      SROC = fClassifier->GetROCIntegral(xbitset.to_string(), methodTitle);

      delete seeddl;
      fClassifier->DeleteAllMethods();
      fClassifier->fMethodsMap.clear();

      for (uint32_t i = 0; i < NBITS; ++i) {
         if (x & (1 << i)) {
            y = x & ~(1 << i);
            std::bitset<NBITS>  ybitset(y);
            //need at least one variable
            //NOTE: if subssed is zero then is the special case
            //that count in xbitset is 1
            Double_t ny = log(x - y) / 0.693147;
            if (y == 0) {
                importances[ny] = SROC - 0.5;
                continue;
            }

            //creating loader for subseed
            TMVA::DataLoader *subseeddl = new TMVA::DataLoader(ybitset.to_string());
            //adding variables from subseed
            for (UInt_t index = 0; index < nbits; index++) {
                if (ybitset[index]) subseeddl->AddVariable(varNames[index], 'F');
            }

            //Loading Dataset
            DataLoaderCopy(subseeddl,fDataLoader.get());

            //Booking SubSeed
            fClassifier->BookMethod(subseeddl, methodName, methodTitle, methodOptions);

            //Train/Test/Evaluation
            fClassifier->TrainAllMethods();
            fClassifier->TestAllMethods();
            fClassifier->EvaluateAllMethods();

            //getting ROC
            SSROC = fClassifier->GetROCIntegral(ybitset.to_string(), methodTitle);
            importances[ny] += SROC - SSROC;

            delete subseeddl;
            fClassifier->DeleteAllMethods();
            fClassifier->fMethodsMap.clear();
        }
    }
    Float_t normalization = 0.0;
    for (UInt_t i = 0; i < nbits; i++) normalization += importances[i];

    for(UInt_t i=0;i<nbits;i++){
        //adding values
        fResults.fImportanceValues[varNames[i]]=(100.0 * importances[i] / normalization);
        //adding sufix
        fResults.fImportanceValues[varNames[i]]=fResults.fImportanceValues.GetValue<TString>(varNames[i])+" % ";
    }
    fResults.fImportanceHist = std::shared_ptr<TH1F>(GetImportance(nbits,importances,varNames));
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::VariableImportance::EvaluateImportanceRandom(UInt_t seeds)
{
   for (auto &meth : fMethods) {

      TString methodName = meth.GetValue<TString>("MethodName");
      TString methodTitle = meth.GetValue<TString>("MethodTitle");
      TString methodOptions = meth.GetValue<TString>("MethodOptions");

      TRandom3 *rangen = new TRandom3(0); // Random Gen.

      uint32_t x = 0;
      uint32_t y = 0;

      // getting number of variables and variable names from loader
      const UInt_t nbits = fDataLoader->GetDefaultDataSetInfo().GetNVariables();
      std::vector<TString> varNames = fDataLoader->GetDefaultDataSetInfo().GetListOfVariables();

      ULong_t range = pow(2, nbits);

      // vector to save importances
      std::vector<Float_t> importances(nbits);
      Float_t importances_norm = 0;

      for (UInt_t i = 0; i < nbits; i++)
         importances[i] = 0;

      Float_t SROC, SSROC; // computed ROC value for every Seed and SubSeed

      x = range;

      for (UInt_t n = 0; n < seeds; n++) {
         x = rangen->Integer(range);

         std::bitset<NBITS> xbitset(x);
         if (x == 0)
            continue; // dataloader need at least one variable

         // creating loader for seed
         TMVA::DataLoader *seeddl = new TMVA::DataLoader(xbitset.to_string());

         // adding variables from seed
         for (UInt_t index = 0; index < nbits; index++) {
            if (xbitset[index]) seeddl->AddVariable(varNames[index], 'F');
        }

        //Loading Dataset
        DataLoaderCopy(seeddl,fDataLoader.get());

        //Booking Seed
        fClassifier->BookMethod(seeddl, methodName, methodTitle, methodOptions);

        //Train/Test/Evaluation
        fClassifier->TrainAllMethods();
        fClassifier->TestAllMethods();
        fClassifier->EvaluateAllMethods();

        //getting ROC
        SROC = fClassifier->GetROCIntegral(xbitset.to_string(), methodTitle);

        delete seeddl;
        fClassifier->DeleteAllMethods();
        fClassifier->fMethodsMap.clear();

        for (uint32_t i = 0; i < 32; ++i) {
            if (x & (1 << i)) {
                y = x & ~(1 << i);
                std::bitset<NBITS>  ybitset(y);
                //need at least one variable
                //NOTE: if subssed is zero then is the special case
                //that count in xbitset is 1
                Double_t ny = log(x - y) / 0.693147;
                if (y == 0) {
                    importances[ny] = SROC - 0.5;
                    importances_norm += importances[ny];
                    continue;
                }

                //creating loader for subseed
                TMVA::DataLoader *subseeddl = new TMVA::DataLoader(ybitset.to_string());
                //adding variables from subseed
                for (UInt_t index = 0; index < nbits; index++) {
                    if (ybitset[index]) subseeddl->AddVariable(varNames[index], 'F');
                }

                //Loading Dataset
                DataLoaderCopy(subseeddl,fDataLoader.get());

                //Booking SubSeed
                fClassifier->BookMethod(subseeddl, methodName, methodTitle, methodOptions);

                //Train/Test/Evaluation
                fClassifier->TrainAllMethods();
                fClassifier->TestAllMethods();
                fClassifier->EvaluateAllMethods();

                //getting ROC
                SSROC = fClassifier->GetROCIntegral(ybitset.to_string(), methodTitle);
                importances[ny] += SROC - SSROC;

                delete subseeddl;
                fClassifier->DeleteAllMethods();
                fClassifier->fMethodsMap.clear();
            }
        }
    }

    Float_t normalization = 0.0;
    for (UInt_t i = 0; i < nbits; i++) normalization += importances[i];

    for(UInt_t i=0;i<nbits;i++){
        //adding values
        fResults.fImportanceValues[varNames[i]]=(100.0 * importances[i] / normalization);
        //adding sufix
        fResults.fImportanceValues[varNames[i]]=fResults.fImportanceValues.GetValue<TString>(varNames[i])+" % ";
    }
    fResults.fImportanceHist = std::shared_ptr<TH1F>(GetImportance(nbits,importances,varNames));
    delete rangen;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::VariableImportance::EvaluateImportanceAll()
{
   for (auto &meth : fMethods) {
      TString methodName = meth.GetValue<TString>("MethodName");
      TString methodTitle = meth.GetValue<TString>("MethodTitle");
      TString methodOptions = meth.GetValue<TString>("MethodOptions");

      uint32_t x = 0;
      uint32_t y = 0;

      // getting number of variables and variable names from loader
      const UInt_t nbits = fDataLoader->GetDefaultDataSetInfo().GetNVariables();
      std::vector<TString> varNames = fDataLoader->GetDefaultDataSetInfo().GetListOfVariables();

      ULong_t range = pow(2, nbits);

      // vector to save importances
      std::vector<Float_t> importances(nbits);

      // vector to save ROC-Integral values
      std::vector<Float_t> ROC(range);
      ROC[0] = 0.5;
      for (UInt_t i = 0; i < nbits; i++)
         importances[i] = 0;

      Float_t SROC, SSROC; // computed ROC value
      for (x = 1; x < range; x++) {

         std::bitset<NBITS> xbitset(x);
         if (x == 0)
            continue; // dataloader need at least one variable

         // creating loader for seed
         TMVA::DataLoader *seeddl = new TMVA::DataLoader(xbitset.to_string());

         // adding variables from seed
         for (UInt_t index = 0; index < nbits; index++) {
            if (xbitset[index]) seeddl->AddVariable(varNames[index], 'F');
        }

        DataLoaderCopy(seeddl,fDataLoader.get());

        seeddl->PrepareTrainingAndTestTree(fDataLoader->GetDefaultDataSetInfo().GetCut("Signal"), fDataLoader->GetDefaultDataSetInfo().GetCut("Background"), fDataLoader->GetDefaultDataSetInfo().GetSplitOptions());

        //Booking Seed
        fClassifier->BookMethod(seeddl, methodName, methodTitle, methodOptions);

        //Train/Test/Evaluation
        fClassifier->TrainAllMethods();
        fClassifier->TestAllMethods();
        fClassifier->EvaluateAllMethods();

        //getting ROC
        ROC[x] = fClassifier->GetROCIntegral(xbitset.to_string(), methodTitle);

        delete seeddl;
        fClassifier->DeleteAllMethods();
        fClassifier->fMethodsMap.clear();
    }


    for ( x = 0; x <range ; x++)
    {
        SROC=ROC[x];
        for (uint32_t i = 0; i < NBITS; ++i) {
            if (x & (1 << i)) {
                y = x & ~(1 << i);
                std::bitset<NBITS>  ybitset(y);

                Float_t ny = log(x - y) / 0.693147;
                if (y == 0) {
                    importances[ny] = SROC - 0.5;
                    continue;
                }

                //getting ROC
                SSROC = ROC[y];
                importances[ny] += SROC - SSROC;
            }

        }
    }
    Float_t normalization = 0.0;
    for (UInt_t i = 0; i < nbits; i++) normalization += importances[i];

    for(UInt_t i=0;i<nbits;i++){
        //adding values
        fResults.fImportanceValues[varNames[i]]=(100.0 * importances[i] / normalization);
        //adding sufix
        fResults.fImportanceValues[varNames[i]]=fResults.fImportanceValues.GetValue<TString>(varNames[i])+" % ";
    }
    fResults.fImportanceHist = std::shared_ptr<TH1F>(GetImportance(nbits,importances,varNames));
   }
}
