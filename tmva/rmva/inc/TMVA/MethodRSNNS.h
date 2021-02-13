// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata,Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodRSNNS                                                          *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package RSNNS  method based on ROOTR                                  *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodRSNNS
#define ROOT_TMVA_RMethodRSNNS

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodRSNNS                                                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/RMethodBase.h"
#include <vector>

namespace TMVA {

   class Factory;  // DSMTEST
   class Reader;   // DSMTEST
   class DataSetManager;  // DSMTEST
   class Types;
   class MethodRSNNS : public RMethodBase {

   public :

      // constructors
      MethodRSNNS(const TString &jobName,
                  const TString &methodTitle,
                  DataSetInfo &theData,
                  const TString &theOption = "");

      MethodRSNNS(DataSetInfo &dsi,
                  const TString &theWeightFile);


      ~MethodRSNNS(void);
      void     Train();
      // options treatment
      void     Init();
      void     DeclareOptions();
      void     ProcessOptions();
      // create ranking
      const Ranking *CreateRanking()
      {
         return NULL;  // = 0;
      }


      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

      // performs classifier testing
      virtual void     TestClassification();


      Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0);

      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /*parent*/) const {}  // = 0;
      virtual void ReadWeightsFromXML(void * /*wghtnode*/) {} // = 0;
      virtual void ReadWeightsFromStream(std::istream &) {} //= 0;       // backward compatibility

      void ReadModelFromFile();

      // signal/background classification response for all current set of data 
      virtual std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);

   private :
      DataSetManager    *fDataSetManager;     // DSMTEST
      friend class Factory;                   // DSMTEST
      friend class Reader;                    // DSMTEST
   protected:
      UInt_t fMvaCounter;
      std::vector<Float_t> fProbResultForTrainSig;
      std::vector<Float_t> fProbResultForTestSig;

      TString fNetType;//default RMPL
      //RSNNS Options for all NN methods
      TString  fSize;//number of units in the hidden layer(s)
      UInt_t fMaxit;//maximum of iterations to learn

      TString fInitFunc;//the initialization function to use
      TString fInitFuncParams;//the parameters for the initialization function (type 6 see getSnnsRFunctionTable() in RSNNS package)

      TString fLearnFunc;//the learning function to use
      TString fLearnFuncParams;//the parameters for the learning function

      TString fUpdateFunc;//the update function to use
      TString fUpdateFuncParams;//the parameters for the update function

      TString fHiddenActFunc;//the activation function of all hidden units
      Bool_t fShufflePatterns;//should the patterns be shuffled?
      Bool_t fLinOut;//sets the activation function of the output units to linear or logistic

      TString fPruneFunc;//the pruning function to use
      TString fPruneFuncParams;//the parameters for the pruning function. Unlike the
      //other functions, these have to be given in a named list. See
      //the pruning demos for further explanation.
      std::vector<UInt_t>  fFactorNumeric;   //factors creations
      //RSNNS mlp require a numeric factor then background=0 signal=1 from fFactorTrain
      static Bool_t IsModuleLoaded;
      ROOT::R::TRFunctionImport predict;
      ROOT::R::TRFunctionImport mlp;
      ROOT::R::TRFunctionImport asfactor;
      ROOT::R::TRObject         *fModel;
      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRSNNS, 0)
   };
} // namespace TMVA
#endif
