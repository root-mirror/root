// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Eckhard von Toerne, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataSetFactory                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Eckhard von Toerne <evt@physik.uni-bonn.de> - U. of Bonn, Germany         *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_DataSetFactory
#define ROOT_TMVA_DataSetFactory

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DataSetFactory                                                       //
//                                                                      //
// Class that contains all the data information                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>

#include "TString.h"
#include "TTree.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TMatrixDfwd.h"
#include "TPrincipal.h"
#include "TRandom3.h"

#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/Event.h"

namespace TMVA {

   class DataSet;
   class DataSetInfo;
   class DataInputHandler;
   class TreeInfo;
   class MsgLogger;

   // =============== maybe move these elswhere (e.g. into the tools )

   // =============== functors =======================

   // delete-functor (to be used in e.g. for_each algorithm)
   template<class T>
      struct DeleteFunctor_t
      {
         DeleteFunctor_t& operator()(const T* p) {
            delete p;
            return *this;
         }
      };

   template<class T>
      DeleteFunctor_t<const T> DeleteFunctor()
      {
         return DeleteFunctor_t<const T>();
      }


   template< typename T >
      class Increment {
      T value;
   public:
   Increment( T start ) : value( start ){ }
      T operator()() {
         return value++;
      }
   };



   template <typename F>
      class null_t
      {
      private:
         // returns argF
      public:
         typedef F argument_type;
         F operator()(const F& argF) const 
         {
            return argF;
         }
      };

   template <typename F>
      inline null_t<F> null() {
      return null_t<F>();
   }

   // =========================================================

   class DataSetFactory:public TObject {

      typedef std::vector<Event* >                             EventVector;
      typedef std::vector< EventVector >                        EventVectorOfClasses;
      typedef std::map<Types::ETreeType, EventVectorOfClasses > EventVectorOfClassesOfTreeType;
      typedef std::map<Types::ETreeType, EventVector >          EventVectorOfTreeType;

      typedef std::vector< Double_t >                    ValuePerClass;
      typedef std::map<Types::ETreeType, ValuePerClass > ValuePerClassOfTreeType;

      class EventStats {
      public:
         Int_t    nTrainingEventsRequested;
         Int_t    nTestingEventsRequested;
         Float_t  TrainTestSplitRequested;
         Int_t    nInitialEvents;
         Int_t    nEvBeforeCut;
         Int_t    nEvAfterCut;
         Float_t  nWeEvBeforeCut;
         Float_t  nWeEvAfterCut;
         Double_t nNegWeights;
         Float_t* varAvLength;//->
      EventStats():
         nTrainingEventsRequested(0),
            nTestingEventsRequested(0),
            TrainTestSplitRequested(0),
            nInitialEvents(0),
            nEvBeforeCut(0),
            nEvAfterCut(0),
            nWeEvBeforeCut(0),
            nWeEvAfterCut(0),
            nNegWeights(0),
            varAvLength(0)
               {}
         ~EventStats() { delete[] varAvLength; }
         Float_t cutScaling() const { return Float_t(nEvAfterCut)/nEvBeforeCut; }
      };

      typedef std::vector< int >                            NumberPerClass;
      typedef std::vector< EventStats >                     EvtStatsPerClass;

   public:

      ~DataSetFactory();

      DataSetFactory();

      DataSet* CreateDataSet( DataSetInfo &, DataInputHandler& );
   protected:
     

      DataSet*  BuildInitialDataSet( DataSetInfo&, TMVA::DataInputHandler& );
      DataSet*  BuildDynamicDataSet( DataSetInfo& );

      // ---------- new versions
      void      BuildEventVector ( DataSetInfo& dsi,
                                   DataInputHandler& dataInput,
                                   EventVectorOfClassesOfTreeType& eventsmap,
                                   EvtStatsPerClass& eventCounts);

      DataSet*  MixEvents        ( DataSetInfo& dsi,
                                   EventVectorOfClassesOfTreeType& eventsmap,
                                   EvtStatsPerClass& eventCounts,
                                   const TString& splitMode,
                                   const TString& mixMode,
                                   const TString& normMode,
                                   UInt_t splitSeed);

      void      RenormEvents     ( DataSetInfo& dsi,
                                   EventVectorOfClassesOfTreeType& eventsmap,
                                   const EvtStatsPerClass& eventCounts,
                                   const TString& normMode );

      void      InitOptions      ( DataSetInfo& dsi,
                                   EvtStatsPerClass& eventsmap,
                                   TString& normMode, UInt_t& splitSeed,
                                   TString& splitMode, TString& mixMode);


      // ------------------------

      // auxiliary functions to compute correlations
      TMatrixD* CalcCorrelationMatrix( DataSet*, const UInt_t classNumber );
      TMatrixD* CalcCovarianceMatrix ( DataSet*, const UInt_t classNumber );
      void      CalcMinMax           ( DataSet*, DataSetInfo& dsi );

      // resets branch addresses to current event
      void   ResetBranchAndEventAddresses( TTree* );
      void   ResetCurrentTree() { fCurrentTree = 0; }
      void   ChangeToNewTree( TreeInfo&, const DataSetInfo & );
      Bool_t CheckTTreeFormula( TTreeFormula* ttf, const TString& expression, Bool_t& hasDollar );

      // verbosity
      Bool_t Verbose() { return fVerbose; }

      // data members

      // verbosity
      Bool_t                     fVerbose;           // Verbosity
      TString                    fVerboseLevel;      // VerboseLevel

      // Printing
      Bool_t fCorrelations = kFALSE;          // Whether to print correlations or not
      Bool_t fComputeCorrelations = kFALSE;   // Whether to force computation of correlations or not

      Bool_t                     fScaleWithPreselEff; // how to deal with requested #events in connection with preselection cuts 

      // the event
      TTree*                     fCurrentTree;       // the tree, events are currently read from
      UInt_t                     fCurrentEvtIdx;     // the current event (to avoid reading of the same event)

      // the formulas for reading the original tree
      std::vector<TTreeFormula*> fInputFormulas;   // input variables
      std::vector<std::pair<TTreeFormula*, Int_t>> fInputTableFormulas;    //! input variables expression for arrays
      std::vector<TTreeFormula *> fTargetFormulas; // targets
      std::vector<TTreeFormula*> fCutFormulas;     // cuts
      std::vector<TTreeFormula*> fWeightFormula;   // weights
      std::vector<TTreeFormula*> fSpectatorFormulas; // spectators

      MsgLogger*                 fLogger;          //! message logger
      MsgLogger& Log() const { return *fLogger; }
   public:
      ClassDef(DataSetFactory, 2);
   };
}

#endif
