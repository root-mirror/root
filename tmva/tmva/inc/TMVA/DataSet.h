// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataSet                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_DataSet
#define ROOT_TMVA_DataSet

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DataSet                                                              //
//                                                                      //
// Class that contains all the data information                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>

#include "TNamed.h"
#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"

#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

namespace TMVA {

   class Event;
   class DataSetInfo;
   class MsgLogger;
   class Results;

   class DataSet :public TNamed {

   public:
      DataSet();
      DataSet(const DataSetInfo&);
      virtual ~DataSet();

      void      AddEvent( Event *, Types::ETreeType );

      Long64_t  GetNEvents( Types::ETreeType type = Types::kMaxTreeType ) const;
      Long64_t  GetNTrainingEvents()              const { return GetNEvents(Types::kTraining); }
      Long64_t  GetNTestEvents()                  const { return GetNEvents(Types::kTesting); }

      // const getters
      const Event*    GetEvent()                        const; // returns event without transformations
      const Event*    GetEvent        ( Long64_t ievt ) const { fCurrentEventIdx = ievt; return GetEvent(); } // returns event without transformations
      const Event*    GetTrainingEvent( Long64_t ievt ) const { return GetEvent(ievt, Types::kTraining); }
      const Event*    GetTestEvent    ( Long64_t ievt ) const { return GetEvent(ievt, Types::kTesting); }
      const Event*    GetEvent        ( Long64_t ievt, Types::ETreeType type ) const 
      {
         fCurrentTreeIdx = TreeIndex(type); fCurrentEventIdx = ievt; return GetEvent();
      }




      UInt_t    GetNVariables()   const;
      UInt_t    GetNTargets()     const;
      UInt_t    GetNSpectators()  const;

      void      SetCurrentEvent( Long64_t ievt         ) const { fCurrentEventIdx = ievt; }
      void      SetCurrentType ( Types::ETreeType type ) const { fCurrentTreeIdx = TreeIndex(type); }
      Types::ETreeType GetCurrentType() const;

      void                       SetEventCollection( std::vector<Event*>*, Types::ETreeType, Bool_t deleteEvents = true );
      const std::vector<Event*>& GetEventCollection( Types::ETreeType type = Types::kMaxTreeType ) const;
      const TTree*               GetEventCollectionAsTree();

      Long64_t  GetNEvtSigTest();
      Long64_t  GetNEvtBkgdTest();
      Long64_t  GetNEvtSigTrain();
      Long64_t  GetNEvtBkgdTrain();

      Bool_t    HasNegativeEventWeights() const { return fHasNegativeEventWeights; }

      Results*  GetResults   ( const TString &,
                               Types::ETreeType type,
                               Types::EAnalysisType analysistype );
      void      DeleteResults   ( const TString &,
                                  Types::ETreeType type,
                                  Types::EAnalysisType analysistype );
      void      DeleteAllResults(Types::ETreeType type,
                                 Types::EAnalysisType analysistype);

      void      SetVerbose( Bool_t ) {}

      // sets the number of blocks to which the training set is divided,
      // some of which are given to the Validation sample. As default they belong all to Training set.
      void      DivideTrainingSet( UInt_t blockNum );

      // sets a certrain block from the origin training set to belong to either Training or Validation set
      void      MoveTrainingBlock( Int_t blockInd,Types::ETreeType dest, Bool_t applyChanges = kTRUE );

      void      IncrementNClassEvents( Int_t type, UInt_t classNumber );
      Long64_t  GetNClassEvents      ( Int_t type, UInt_t classNumber );
      void      ClearNClassEvents    ( Int_t type );

      TTree*    GetTree( Types::ETreeType type );

      // accessors for random and importance sampling
      void      InitSampling( Float_t fraction, Float_t weight, UInt_t seed = 0 );
      void      EventResult( Bool_t successful, Long64_t evtNumber = -1 );
      void      CreateSampling() const;

      UInt_t    TreeIndex(Types::ETreeType type) const;

   private:

      // data members
      void DestroyCollection( Types::ETreeType type, Bool_t deleteEvents );

      const DataSetInfo         *fdsi;                       //-> datasetinfo that created this dataset

      std::vector< std::vector<Event*>  > fEventCollection;  // list of events for training/testing/...

      std::vector< std::map< TString, Results* > > fResults; //!  [train/test/...][method-identifier]

      mutable UInt_t             fCurrentTreeIdx;
      mutable Long64_t           fCurrentEventIdx;

      // event sampling
      std::vector<Char_t>        fSampling;                   // random or importance sampling (not all events are taken) !! Bool_t are stored ( no std::vector<bool> taken for speed (performance) issues )
      std::vector<Int_t>         fSamplingNEvents;            // number of events which should be sampled
      std::vector<Float_t>       fSamplingWeight;             // weight change factor [weight is indicating if sampling is random (1.0) or importance (<1.0)] 
      mutable std::vector< std::vector< std::pair< Float_t, Long64_t > > > fSamplingEventList;  // weights and indices for sampling
      mutable std::vector< std::vector< std::pair< Float_t, Long64_t > > > fSamplingSelected;   // selected events
      TRandom3                   *fSamplingRandom;             //-> random generator for sampling


      // further things
      std::vector< std::vector<Long64_t> > fClassEvents;       // number of events of class 0,1,2,... in training[0] 
                                                               // and testing[1] (+validation, trainingoriginal)

      Bool_t                     fHasNegativeEventWeights;     // true if at least one signal or bkg event has negative weight

      mutable MsgLogger*         fLogger;                      //! message logger
      MsgLogger& Log() const { return *fLogger; }
      std::vector<Char_t>        fBlockBelongToTraining;       // when dividing the dataset to blocks, sets whether 
                                                               // the certain block is in the Training set or else 
                                                               // in the validation set 
                                                               // boolean are stored, taken std::vector<Char_t> for performance reasons (instead of std::vector<Bool_t>)
      Long64_t                   fTrainingBlockSize;           // block size into which the training dataset is divided

      void  ApplyTrainingBlockDivision();
      void  ApplyTrainingSetDivision();
   public:
       
       ClassDef(DataSet,1);
   };
}


//_______________________________________________________________________
inline UInt_t TMVA::DataSet::TreeIndex(Types::ETreeType type) const
{
   switch (type) {
   case Types::kMaxTreeType : return fCurrentTreeIdx;
   case Types::kTraining : return 0;
   case Types::kTesting : return 1;
   case Types::kValidation : return 2;
   case Types::kTrainingOriginal : return 3;
   default : return fCurrentTreeIdx;
   }
}

//_______________________________________________________________________
inline TMVA::Types::ETreeType TMVA::DataSet::GetCurrentType() const
{
   switch (fCurrentTreeIdx) {
   case 0: return Types::kTraining;
   case 1: return Types::kTesting;
   case 2: return Types::kValidation;
   case 3: return Types::kTrainingOriginal;
   }
   return Types::kMaxTreeType;
}

//_______________________________________________________________________
inline Long64_t TMVA::DataSet::GetNEvents(Types::ETreeType type) const 
{
   Int_t treeIdx = TreeIndex(type);
   if (fSampling.size() > UInt_t(treeIdx) && fSampling.at(treeIdx)) {
      return fSamplingSelected.at(treeIdx).size();
   }
   return GetEventCollection(type).size();
}

//_______________________________________________________________________
inline const std::vector<TMVA::Event*>& TMVA::DataSet::GetEventCollection( TMVA::Types::ETreeType type ) const
{
   return fEventCollection.at(TreeIndex(type));
}


#endif
