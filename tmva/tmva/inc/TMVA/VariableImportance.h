// @(#)root/tmva $Id$
// Author: Omar Zapata and Sergei Gleyzer. 2016


#ifndef ROOT_TMVA_VariableImportance
#define ROOT_TMVA_VariableImportance


#include "TString.h"
#include <vector>

#include "TMVA/Configurable.h"
#include "TMVA/Types.h"

#include <TMVA/Factory.h>

#include <TMVA/DataLoader.h>

#include <TMVA/OptionMap.h>

#include <TMVA/Envelope.h>

namespace TMVA {

   class VariableImportanceResult
   {
     friend class VariableImportance;
   private:
       OptionMap              fImportanceValues;
       std::shared_ptr<TH1F>  fImportanceHist;
       VIType                 fType {kShort};
   public:
       VariableImportanceResult();
       VariableImportanceResult(const VariableImportanceResult &);
       ~VariableImportanceResult(){fImportanceHist=nullptr;}

       OptionMap &GetImportanceValues(){return fImportanceValues;}
       TH1F *GetImportanceHist(){return fImportanceHist.get();}
       void Print() const ;

       TCanvas* Draw(const TString name="VariableImportance") const;
   };


   class VariableImportance : public Envelope {
   private:
       UInt_t                    fNumFolds = 0;
       VariableImportanceResult  fResults;
       VIType                    fType {kShort};
   public:
       explicit VariableImportance(DataLoader *loader);
       ~VariableImportance();

       virtual void Evaluate();

       void SetType(VIType type){fType=type;}
       VIType GetType(){return fType;}

       const VariableImportanceResult& GetResults() const {return fResults;}//I need to think about this, which is the best way to get the results?
   protected:
       //evaluate the simple case that is removing 1 variable at time
       void EvaluateImportanceShort();
       //evaluate all variables combinations NOTE: use with care in huge datasets with a huge number of variables
       void EvaluateImportanceAll();
       //evaluate randomly given a number of seeds
       void EvaluateImportanceRandom(UInt_t nseeds);

       //method to return a nice histogram with the results ;)
       TH1F* GetImportance(const UInt_t nbits,std::vector<Float_t> &importances,std::vector<TString> &varNames);

       //method to compute the range(number total of operations for every bit configuration)
       ULong_t Sum(ULong_t i);

   private:
       std::unique_ptr<Factory>     fClassifier;
       ClassDef(VariableImportance,0);
   };
}


#endif
