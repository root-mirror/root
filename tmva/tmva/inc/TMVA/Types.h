// @(#)root/tmva $Id$   
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Types                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      GLobal types (singleton class)                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#ifndef ROOT_TMVA_Types
#define ROOT_TMVA_Types

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Types                                                                //
//                                                                      //
// Singleton class for Global types used by TMVA                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <map>
#if __cplusplus > 199711L
#include <atomic>
#endif

#include "RtypesCore.h"

#include "TString.h"

namespace TMVA {

   typedef UInt_t TMVAVersion_t;

   class MsgLogger;

   // message types for MsgLogger
   // define outside of Types class to facilite access
   enum EMsgType {
      kDEBUG   = 1,
      kVERBOSE = 2,
      kINFO    = 3,
      kWARNING = 4,
      kERROR   = 5,
      kFATAL   = 6,
      kSILENT  = 7,
      kHEADER  = 8
   };

   enum HistType { kMVAType = 0, kProbaType = 1, kRarityType = 2, kCompareType = 3 };

   //Variable Importance type
   enum VIType {kShort=0,kAll=1,kRandom=2};
   
   class Types {

   public:

      // available MVA methods
      enum EMVA {
         kVariable    = 0,
         kCuts           ,
         kLikelihood     ,
         kPDERS          ,
         kHMatrix        ,
         kFisher         ,
         kKNN            ,
         kCFMlpANN       ,
         kTMlpANN        ,
         kBDT            ,
         kDT             ,
         kRuleFit        ,
         kSVM            ,
         kMLP            ,
         kBayesClassifier,
         kFDA            ,
         kBoost          ,
         kPDEFoam        ,
         kLD             ,
         kPlugins        ,
         kCategory       ,
         kDNN            ,
         kDL             ,
         kPyRandomForest ,
         kPyAdaBoost     ,
         kPyGTB          ,
         kPyKeras        ,
         kPyTorch        ,
         kC50            ,
         kRSNNS          ,
         kRSVM           ,
         kRXGB           ,
         kCrossValidation,
         kMaxMethod
      };

      // available variable transformations
      enum EVariableTransform {
         kIdentity = 0,
         kDecorrelated,
         kNormalized,
         kPCA,
         kRearranged,
         kGauss,
         kUniform,
         kMaxVariableTransform
      };

      // type of analysis
      enum EAnalysisType {
         kClassification = 0,
         kRegression,
         kMulticlass,
         kNoAnalysisType,
         kMaxAnalysisType
      };

      enum ESBType {
         kSignal = 0,  // Never change this number - it is elsewhere assumed to be zero !
         kBackground,
         kSBBoth,
         kMaxSBType,
         kTrueType
      };

      enum ETreeType {
         kTraining = 0,
         kTesting,
         kMaxTreeType,  // also used as temporary storage for trees not yet assigned for testing;training... 
         kValidation,   // these are placeholders... currently not used, but could be moved "forward" if
         kTrainingOriginal     // ever needed 
      };

      enum EBoostStage {
         kBoostProcBegin=0,
         kBeforeTraining,
         kBeforeBoosting,
         kAfterBoosting,
         kBoostProcEnd
      };

   public:

      static Types& Instance();
      static void   DestroyInstance();
      ~Types();

      Types::EMVA   GetMethodType( const TString& method ) const;
      TString       GetMethodName( Types::EMVA    method ) const;

      Bool_t        AddTypeMapping(Types::EMVA method, const TString& methodname);

   private:

      Types();
#if __cplusplus > 199711L && !defined _MSC_VER
      static std::atomic<Types*> fgTypesPtr;
#else
      static Types* fgTypesPtr;
#endif

   private:

      std::map<TString, TMVA::Types::EMVA> fStr2type; // types-to-text map
      mutable MsgLogger* fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }

   };
}

#endif
