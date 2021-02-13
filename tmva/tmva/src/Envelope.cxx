// @(#)root/tmva $Id$
// Author: Omar Zapata, Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <TMVA/Envelope.h>

#include <TMVA/Configurable.h>
#include <TMVA/DataLoader.h>
#include <TMVA/MethodBase.h>
#include <TMVA/OptionMap.h>
#include <TMVA/ResultsClassification.h>
#include <TMVA/Types.h>

#include <TMVA/VariableInfo.h>
#include <TMVA/VariableTransform.h>

#include <TAxis.h>
#include <TFile.h>
#include <TH2.h>

using namespace TMVA;

//_______________________________________________________________________
/**
Constructor for the initialization of Envelopes,
differents Envelopes may needs differents constructors then
this is a generic one protected.
\param name the name algorithm.
\param dataloader TMVA::DataLoader object with the data.
\param file optional file to save the results.
\param options extra options for the algorithm.
*/
Envelope::Envelope(const TString &name, DataLoader *dataloader, TFile *file, const TString options)
   : Configurable(options), fDataLoader(dataloader), fFile(file), fModelPersistence(kTRUE), fVerbose(kFALSE),
     fTransformations("I"), fSilentFile(kFALSE), fJobs(1)
{
    SetName(name.Data());
    // render silent
    if (gTools().CheckForSilentOption(GetOptions()))
       Log().InhibitOutput(); // make sure is silent if wanted to

    fModelPersistence = kTRUE;
    DeclareOptionRef(fVerbose, "V", "Verbose flag");

    DeclareOptionRef(fModelPersistence, "ModelPersistence",
                     "Option to save the trained model in xml file or using serialization");
    DeclareOptionRef(fTransformations, "Transformations", "List of transformations to test; formatting example: "
                                                          "\"Transformations=I;D;P;U;G,D\", for identity, "
                                                          "decorrelation, PCA, Uniform and Gaussianisation followed by "
                                                          "decorrelation transformations");
    DeclareOptionRef(fJobs, "Jobs", "Option to run hign level algorithms in parallel with multi-thread");
}

//_______________________________________________________________________
Envelope::~Envelope()
{
}

//_______________________________________________________________________
/**
Method to see if a file is available to save results
\return Boolean with the status.
*/
Bool_t  Envelope::IsSilentFile(){return fFile==nullptr;}

//_______________________________________________________________________
/**
Method to get the pointer to TFile object.
\return pointer to TFile object.
*/
TFile* Envelope::GetFile(){return fFile.get();}

//_______________________________________________________________________
/**
Method to set the pointer to TFile object,
with a writable file.
\param file pointer to TFile object.
*/
void   Envelope::SetFile(TFile *file){fFile=std::shared_ptr<TFile>(file);}

//_______________________________________________________________________
/**
Method to see if the algorithm should print extra information.
\return Boolean with the status.
*/
Bool_t Envelope::IsVerbose(){return fVerbose;}

//_______________________________________________________________________
/**
Method enable print extra information in the algorithms.
\param status Boolean with the status.
*/
void Envelope::SetVerbose(Bool_t status){fVerbose=status;}

//_______________________________________________________________________
/**
Method get the Booked methods in a option map object.
\return vector of TMVA::OptionMap objects with the information of the Booked method
*/
std::vector<OptionMap> &Envelope::GetMethods()
{
   return fMethods;
}

//_______________________________________________________________________
/**
Method to get the pointer to TMVA::DataLoader object.
\return  pointer to TMVA::DataLoader object.
*/

DataLoader *Envelope::GetDataLoader(){    return fDataLoader.get();}

//_______________________________________________________________________
/**
Method to set the pointer to TMVA::DataLoader object.
\param dataloader pointer to TMVA::DataLoader object.
*/

void Envelope::SetDataLoader(DataLoader *dataloader)
{
   fDataLoader = std::shared_ptr<DataLoader>(dataloader);
}

//_______________________________________________________________________
/**
Method to see if the algorithm model is saved in xml or serialized files.
\return Boolean with the status.
*/
Bool_t TMVA::Envelope::IsModelPersistence(){return fModelPersistence; }

//_______________________________________________________________________
/**
Method enable model persistence, then algorithms model is saved in xml or serialized files.
\param status Boolean with the status.
*/
void TMVA::Envelope::SetModelPersistence(Bool_t status){fModelPersistence=status;}

//_______________________________________________________________________
/**
Method to book the machine learning method to perform the algorithm.
\param method enum TMVA::Types::EMVA with the type of the mva method
\param methodTitle String with the method title.
\param options String with the options for the method.
*/
void TMVA::Envelope::BookMethod(Types::EMVA method, TString methodTitle, TString options){
   BookMethod(Types::Instance().GetMethodName(method), methodTitle, options);
}

//_______________________________________________________________________
/**
Method to book the machine learning method to perform the algorithm.
\param methodName String with the name of the mva method
\param methodTitle String with the method title.
\param options String with the options for the method.
*/
void TMVA::Envelope::BookMethod(TString methodName, TString methodTitle, TString options){
   for (auto &meth : fMethods) {
      if (meth.GetValue<TString>("MethodName") == methodName && meth.GetValue<TString>("MethodTitle") == methodTitle) {
         Log() << kFATAL << "Booking failed since method with title <" << methodTitle << "> already exists "
               << "in with DataSet Name <" << fDataLoader->GetName() << ">  " << Endl;
      }
   }
   OptionMap fMethod;
   fMethod["MethodName"] = methodName;
   fMethod["MethodTitle"] = methodTitle;
   fMethod["MethodOptions"] = options;

   fMethods.push_back(fMethod);
}

//_______________________________________________________________________
/**
Method to parse the internal option string.
*/
void TMVA::Envelope::ParseOptions()
{

   Bool_t silent = kFALSE;
#ifdef WIN32
   // under Windows, switch progress bar and color off by default, as the typical windows shell doesn't handle these
   // (would need different sequences..)
   Bool_t color = kFALSE;
   Bool_t drawProgressBar = kFALSE;
#else
   Bool_t color = !gROOT->IsBatch();
   Bool_t drawProgressBar = kTRUE;
#endif
   DeclareOptionRef(color, "Color", "Flag for coloured screen output (default: True, if in batch mode: False)");
   DeclareOptionRef(drawProgressBar, "DrawProgressBar",
                    "Draw progress bar to display training, testing and evaluation schedule (default: True)");
   DeclareOptionRef(silent, "Silent", "Batch mode: boolean silent flag inhibiting any output from TMVA after the "
                                      "creation of the factory class object (default: False)");

   Configurable::ParseOptions();
   CheckForUnusedOptions();

   if (IsVerbose())
      Log().SetMinType(kVERBOSE);

   // global settings
   gConfig().SetUseColor(color);
   gConfig().SetSilent(silent);
   gConfig().SetDrawProgressBar(drawProgressBar);
}

//_______________________________________________________________________
/**
 * function to check methods booked
 * \param methodname  Method's name.
 * \param methodtitle title associated to the method.
 * \return true if the method was booked.
 */
Bool_t TMVA::Envelope::HasMethod(TString methodname, TString methodtitle)
{
   for (auto &meth : fMethods) {
      if (meth.GetValue<TString>("MethodName") == methodname && meth.GetValue<TString>("MethodTitle") == methodtitle)
         return kTRUE;
   }
   return kFALSE;
}

//_______________________________________________________________________
/**
 * method to save Train/Test information into the output file.
 * \param fDataSetInfo TMVA::DataSetInfo object reference
 * \param fAnalysisType Types::kMulticlass and Types::kRegression
 */
void TMVA::Envelope::WriteDataInformation(TMVA::DataSetInfo &fDataSetInfo, TMVA::Types::EAnalysisType fAnalysisType)
{
   RootBaseDir()->cd();

   if (!RootBaseDir()->GetDirectory(fDataSetInfo.GetName()))
      RootBaseDir()->mkdir(fDataSetInfo.GetName());
   else
      return; // loader is now in the output file, we dont need to save again

   RootBaseDir()->cd(fDataSetInfo.GetName());
   fDataSetInfo.GetDataSet(); // builds dataset (including calculation of correlation matrix)

   // correlation matrix of the default DS
   const TMatrixD *m(0);
   const TH2 *h(0);

   if (fAnalysisType == Types::kMulticlass) {
      for (UInt_t cls = 0; cls < fDataSetInfo.GetNClasses(); cls++) {
         m = fDataSetInfo.CorrelationMatrix(fDataSetInfo.GetClassInfo(cls)->GetName());
         h = fDataSetInfo.CreateCorrelationMatrixHist(
            m, TString("CorrelationMatrix") + fDataSetInfo.GetClassInfo(cls)->GetName(),
            TString("Correlation Matrix (") + fDataSetInfo.GetClassInfo(cls)->GetName() + TString(")"));
         if (h != 0) {
            h->Write();
            delete h;
         }
      }
   } else {
      m = fDataSetInfo.CorrelationMatrix("Signal");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrixS", "Correlation Matrix (signal)");
      if (h != 0) {
         h->Write();
         delete h;
      }

      m = fDataSetInfo.CorrelationMatrix("Background");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrixB", "Correlation Matrix (background)");
      if (h != 0) {
         h->Write();
         delete h;
      }

      m = fDataSetInfo.CorrelationMatrix("Regression");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrix", "Correlation Matrix");
      if (h != 0) {
         h->Write();
         delete h;
      }
   }

   // some default transformations to evaluate
   // NOTE: all transformations are destroyed after this test
   TString processTrfs = "I"; //"I;N;D;P;U;G,D;"

   // plus some user defined transformations
   processTrfs = fTransformations;

   // remove any trace of identity transform - if given (avoid to apply it twice)
   std::vector<TMVA::TransformationHandler *> trfs;
   TransformationHandler *identityTrHandler = 0;

   std::vector<TString> trfsDef = gTools().SplitString(processTrfs, ';');
   std::vector<TString>::iterator trfsDefIt = trfsDef.begin();
   for (; trfsDefIt != trfsDef.end(); ++trfsDefIt) {
      trfs.push_back(new TMVA::TransformationHandler(fDataSetInfo, "Envelope"));
      TString trfS = (*trfsDefIt);

      // Log() << kINFO << Endl;
      Log() << kDEBUG << "current transformation string: '" << trfS.Data() << "'" << Endl;
      TMVA::CreateVariableTransforms(trfS, fDataSetInfo, *(trfs.back()), Log());

      if (trfS.BeginsWith('I'))
         identityTrHandler = trfs.back();
   }

   const std::vector<Event *> &inputEvents = fDataSetInfo.GetDataSet()->GetEventCollection();

   // apply all transformations
   std::vector<TMVA::TransformationHandler *>::iterator trfIt = trfs.begin();

   for (; trfIt != trfs.end(); ++trfIt) {
      // setting a Root dir causes the variables distributions to be saved to the root file
      (*trfIt)->SetRootDir(RootBaseDir()->GetDirectory(fDataSetInfo.GetName())); // every dataloader have its own dir
      (*trfIt)->CalcTransformations(inputEvents);
   }
   if (identityTrHandler)
      identityTrHandler->PrintVariableRanking();

   // clean up
   for (trfIt = trfs.begin(); trfIt != trfs.end(); ++trfIt)
      delete *trfIt;
}
