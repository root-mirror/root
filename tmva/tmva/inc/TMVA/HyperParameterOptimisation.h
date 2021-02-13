// @(#)root/tmva $Id$
// Author: Omar Zapata, Thomas James Stevenson.


#ifndef ROOT_TMVA_HyperParameterOptimisation
#define ROOT_TMVA_HyperParameterOptimisation


#include "TString.h"
#include <vector>
#include <map>

#include "TMultiGraph.h"

#include "TMVA/IMethod.h"
#include "TMVA/Configurable.h"
#include "TMVA/Types.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include <TMVA/Results.h>

#include <TMVA/Factory.h>

#include <TMVA/DataLoader.h>

#include <TMVA/Envelope.h>

namespace TMVA {

   class HyperParameterOptimisationResult
   {
     friend class HyperParameterOptimisation;
   private:
      Float_t fROCAVG;
      std::vector<Float_t> fROCs;
      std::vector<Double_t> fSigs;
      std::vector<Double_t> fSeps;
      std::vector<Double_t> fEff01s;
      std::vector<Double_t> fEff10s;
      std::vector<Double_t> fEff30s;
      std::vector<Double_t> fEffAreas;
      std::vector<Double_t> fTrainEff01s;
      std::vector<Double_t> fTrainEff10s;
      std::vector<Double_t> fTrainEff30s;
      std::shared_ptr<TMultiGraph> fROCCurves;
      TString fMethodName;

   public:
       HyperParameterOptimisationResult();
       ~HyperParameterOptimisationResult();

       std::vector<std::map<TString,Double_t> > fFoldParameters;
       
       std::vector<Float_t> GetROCValues(){return fROCs;}
       Float_t GetROCAverage(){return fROCAVG;}
       TMultiGraph *GetROCCurves(Bool_t fLegend=kTRUE);

       void Print() const ;       
//        TCanvas* Draw(const TString name="HyperParameterOptimisation") const;
       
       std::vector<Double_t> GetSigValues(){return fSigs;}
       std::vector<Double_t> GetSepValues(){return fSeps;}
       std::vector<Double_t> GetEff01Values(){return fEff01s;}
       std::vector<Double_t> GetEff10Values(){return fEff10s;}
       std::vector<Double_t> GetEff30Values(){return fEff30s;}
       std::vector<Double_t> GetEffAreaValues(){return fEffAreas;}
       std::vector<Double_t> GetTrainEff01Values(){return fTrainEff01s;}
       std::vector<Double_t> GetTrainEff10Values(){return fTrainEff10s;}
       std::vector<Double_t> GetTrainEff30Values(){return fTrainEff30s;}

   };
    
   class HyperParameterOptimisation : public Envelope {
   public:

       HyperParameterOptimisation(DataLoader *dataloader);
       ~HyperParameterOptimisation();
       
       void SetFitter(TString fitType){fFitType=fitType;}
       TString GetFiiter(){return fFitType;}
       
       
       //Figure of Merit (FOM) default Separation
       void SetFOMType(TString ftype){fFomType=ftype;}
       TString GetFOMType(){return fFitType;}
       
       void SetNumFolds(UInt_t folds);
       UInt_t GetNumFolds(){return fNumFolds;}
       
       virtual void Evaluate();
       const HyperParameterOptimisationResult& GetResults() const {return fResults;}
       
       
   private:
       TString                           fFomType;     //!
       TString                           fFitType;     //!
       UInt_t                            fNumFolds;    //!
       Bool_t                            fFoldStatus;  //!
       HyperParameterOptimisationResult  fResults;     //!
       std::unique_ptr<Factory>          fClassifier;  //!

   public:
       ClassDef(HyperParameterOptimisation,0);  
   };
} 


#endif



