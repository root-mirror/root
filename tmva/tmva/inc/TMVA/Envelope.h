// @(#)root/tmva:$Id$
// Author: Omar Zapata   2016

/*************************************************************************
 * Copyright (C) 2016, Omar Andres Zapata Mesa                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TMVA_Envelope
#define ROOT_TMVA_Envelope

#include <memory>
#include <vector>

#include <TString.h>
#include <TROOT.h>
#include <TStopwatch.h>

#ifndef _MSC_VER
#include <TProcPool.h>
#endif

#include <TMVA/OptionMap.h>
#include <TMVA/Config.h>
#include <TMVA/Tools.h>
#include <TMVA/DataLoader.h>

/*! \class TMVA::Envelope
 * Abstract base class for all high level ml algorithms,
 * you can book ml methods like BDT, MLP. SVM etc..
 * and set a TMVA::DataLoader object to run your code
 * in the overloaded method Evaluate.
\ingroup TMVA

Base class for all machine learning algorithms

*/

namespace TMVA {

      class Envelope:public Configurable
      {
      protected:
         std::vector<OptionMap> fMethods;         //! Booked method information
         std::shared_ptr<DataLoader> fDataLoader; //! data
         std::shared_ptr<TFile> fFile;            //! file to save the results
         Bool_t fModelPersistence;                //! flag to save the trained model
         Bool_t fVerbose;                         //! flag for extra information
         TString fTransformations;                //! List of transformations to test
         Bool_t fSilentFile;                      //! if true dont produce file output
#ifndef _MSC_VER
         TProcPool fWorkers;                      //! procpool object
#endif
         UInt_t fJobs;                            //! number of jobs to run some high level algorithm in parallel
         TStopwatch fTimer;                       //! timer to measute the time.

         Envelope(const TString &name, DataLoader *dataloader = nullptr, TFile *file = nullptr,
                  const TString options = "");

      public:
          /**
           Default destructor
           */
          ~Envelope();

          virtual void BookMethod( TString methodname, TString methodtitle, TString options = "");
          virtual void BookMethod( Types::EMVA method,  TString methodtitle, TString options = "");

          // parse the internal option string
          virtual void ParseOptions();

          Bool_t  IsSilentFile();
          TFile* GetFile();
          void   SetFile(TFile *file);
          Bool_t HasMethod(TString methodname, TString methodtitle);

          DataLoader *GetDataLoader();
          void SetDataLoader(DataLoader *dalaloader);
          Bool_t IsModelPersistence();
          void SetModelPersistence(Bool_t status=kTRUE);
          Bool_t IsVerbose();
          void SetVerbose(Bool_t status);

          /**
            Virtual method to be implemented with your algorithm.
          */
          virtual void Evaluate() = 0;

          std::vector<OptionMap> &GetMethods();

       protected:
          /**
            Utility method to get TMVA::DataInputHandler reference from the DataLoader.
            \return TMVA::DataInputHandler reference.
          */
          DataInputHandler &GetDataLoaderDataInput() { return fDataLoader->DataInput(); }

          /**
            Utility method to get TMVA::DataSetInfo reference from the DataLoader.
            \return TMVA::DataSetInfo reference.
          */
          DataSetInfo &GetDataLoaderDataSetInfo() { return fDataLoader->GetDataSetInfo(); }

          /**
            Utility method to get TMVA::DataSetManager pointer from the DataLoader.
            \return TMVA::DataSetManager pointer.
          */
         DataSetManager *GetDataLoaderDataSetManager()
         {
            return fDataLoader->GetDataSetInfo().GetDataSetManager();
         }

          /**
            Utility method to get base dir directory from current file.
            \return TDirectory* pointer.
          */
          TDirectory *RootBaseDir() { return (TDirectory *)fFile.get(); }

          void WriteDataInformation(TMVA::DataSetInfo &fDataSetInfo, TMVA::Types::EAnalysisType fAnalysisType);

          ClassDef(Envelope, 0);
      };
}

#endif
