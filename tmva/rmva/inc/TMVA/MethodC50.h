// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata,Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodC50                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package C50  method based on ROOTR                                    *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodC50
#define ROOT_TMVA_RMethodC50

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodC50                                                          //
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
   class MethodC50 : public RMethodBase {

   public :

      // constructors
      MethodC50(const TString &jobName,
                const TString &methodTitle,
                DataSetInfo &theData,
                const TString &theOption = "");

      MethodC50(DataSetInfo &dsi,
                const TString &theWeightFile);


      ~MethodC50(void);
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
      virtual void     MakeClass(const TString &classFileName = TString("")) const;  //required for model persistence
      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /*parent*/) const {} // = 0;
      virtual void ReadWeightsFromXML(void * /*weight*/) {} // = 0;
      virtual void ReadWeightsFromStream(std::istream &) {} //= 0;       // backward compatibility

      // signal/background classification response for all current set of data 
      virtual std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);

      void ReadModelFromFile();
   private :
      DataSetManager    *fDataSetManager;     // DSMTEST
      friend class Factory;                   // DSMTEST
      friend class Reader;                    // DSMTEST
   protected:
      //C5.0 function options
      UInt_t fNTrials;//number of trials with boost enabled
      Bool_t fRules;//A logical: should the tree be decomposed into a rule-based model?

      //Control options see C5.0Control
      Bool_t fControlSubset; //A logical: should the model evaluate groups of discrete predictors for splits?
      UInt_t fControlBands;
      Bool_t fControlWinnow;// A logical: should predictor winnowing (i.e feature selection) be used?
      Bool_t fControlNoGlobalPruning; //A logical to toggle whether the final, global pruning step to simplify the tree.
      Double_t fControlCF; //A number in (0, 1) for the confidence factor.
      UInt_t fControlMinCases;//an integer for the smallest number of samples that must be put in at least two of the splits.
      Bool_t fControlFuzzyThreshold;//A logical toggle to evaluate possible advanced splits of the data. See Quinlan (1993) for details and examples.
      Double_t fControlSample;//A value between (0, .999) that specifies the random proportion of the data should be used to train the model.
      Int_t fControlSeed;//An integer for the random number seed within the C code.
      Bool_t fControlEarlyStopping;// logical to toggle whether the internal method for stopping boosting should be used.

      UInt_t fMvaCounter;
      static Bool_t IsModuleLoaded;

      ROOT::R::TRFunctionImport predict;
      ROOT::R::TRFunctionImport C50;
      ROOT::R::TRFunctionImport C50Control;
      ROOT::R::TRFunctionImport asfactor;
      ROOT::R::TRObject *fModel;
      ROOT::R::TRObject fModelControl;
      std::vector <TString > ListOfVariables;


      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodC50, 0)
   };
} // namespace TMVA
#endif
