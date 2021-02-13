// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata,Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodRXGB                                                           *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package xgboost  method based on ROOTR                                *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodXGB
#define ROOT_TMVA_RMethodXGB

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodRXGB                                                          //
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
   class MethodRXGB: public RMethodBase {

   public :

      // constructors
      MethodRXGB(const TString &jobName,
                 const TString &methodTitle,
                 DataSetInfo &theData,
                 const TString &theOption = "");

      MethodRXGB(DataSetInfo &dsi,
                 const TString &theWeightFile);


      ~MethodRXGB(void);
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

      
      //RXGBfunction options
      //https://github.com/dmlc/xgboost/blob/master/doc/parameter.md
      UInt_t fNRounds;
      Double_t fEta;
      UInt_t fMaxDepth;
      static Bool_t IsModuleLoaded;

      std::vector<UInt_t>  fFactorNumeric;   //factors creations
      //xgboost  require a numeric factor then background=0 signal=1 from fFactorTrain

      ROOT::R::TRFunctionImport predict;
      ROOT::R::TRFunctionImport xgbtrain;
      ROOT::R::TRFunctionImport xgbdmatrix;
      ROOT::R::TRFunctionImport xgbsave;
      ROOT::R::TRFunctionImport xgbload;
      ROOT::R::TRFunctionImport asfactor;
      ROOT::R::TRFunctionImport asmatrix;
      ROOT::R::TRObject *fModel;


      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRXGB, 0)
   };
} // namespace TMVA
#endif
