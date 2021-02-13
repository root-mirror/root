// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard von Toerne, Jan Therhaag
// Updated by: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Factory                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      This is the main MVA steering class: it creates (books) all MVA methods,  *
 *      and guides them through the training, testing and evaluation phases.      *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <stelzer@cern.ch>        - DESY, Germany                  *
 *      Peter Speckmayer <peter.speckmayer@cern.ch> - CERN, Switzerland           *
 *      Jan Therhaag          <Jan.Therhaag@cern.ch>   - U of Bonn, Germany       *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Omar Zapata     <Omar.Zapata@cern.ch>    - UdeA/ITM Colombia              *
 *      Lorenzo Moneta  <Lorenzo.Moneta@cern.ch> - CERN, Switzerland              *
 *      Sergei Gleyzer  <Sergei.Gleyzer@cern.ch> - U of Florida & CERN            *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *      UdeA/ITM, Colombia                                                        *
 *      U. of Florida, USA                                                        *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_Factory
#define ROOT_TMVA_Factory

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Factory                                                              //
//                                                                      //
// This is the main MVA steering class: it creates all MVA methods,     //
// and guides them through the training, testing and evaluation         //
// phases                                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include "TCut.h"

#include "TMVA/Configurable.h"
#include "TMVA/Types.h"
#include "TMVA/DataSet.h"

class TCanvas;
class TDirectory;
class TFile;
class TGraph;
class TH1F;
class TMultiGraph;
class TTree;
namespace TMVA {

   class IMethod;
   class MethodBase;
   class DataInputHandler;
   class DataSetInfo;
   class DataSetManager;
   class DataLoader;
   class ROCCurve;
   class VariableTransformBase;


   class Factory : public Configurable {
      friend class CrossValidation;
   public:

      typedef std::vector<IMethod*> MVector;
      std::map<TString,MVector*>  fMethodsMap;//all methods for every dataset with the same name

      // no default  constructor
      Factory( TString theJobName, TFile* theTargetFile, TString theOption = "" );

      // contructor to work without file
      Factory( TString theJobName, TString theOption = "" );

      // default destructor
      virtual ~Factory();

      // use TName::GetName and define correct name in constructor
      //virtual const char*  GetName() const { return "Factory"; }


      MethodBase* BookMethod( DataLoader *loader, TString theMethodName, TString methodTitle, TString theOption = "" );
      MethodBase* BookMethod( DataLoader *loader, Types::EMVA theMethod,  TString methodTitle, TString theOption = "" );
      MethodBase* BookMethod( DataLoader *, TMVA::Types::EMVA /*theMethod*/,
                              TString /*methodTitle*/,
                              TString /*methodOption*/,
                              TMVA::Types::EMVA /*theComposite*/,
                              TString /*compositeOption = ""*/ ) { return 0; }

      // optimize all booked methods (well, if desired by the method)
      std::map<TString,Double_t> OptimizeAllMethods                 (TString fomType="ROCIntegral", TString fitType="FitGA");
      void OptimizeAllMethodsForClassification(TString fomType="ROCIntegral", TString fitType="FitGA") { OptimizeAllMethods(fomType,fitType); }
      void OptimizeAllMethodsForRegression    (TString fomType="ROCIntegral", TString fitType="FitGA") { OptimizeAllMethods(fomType,fitType); }

      // training for all booked methods
      void TrainAllMethods                 ();
      void TrainAllMethodsForClassification( void ) { TrainAllMethods(); }
      void TrainAllMethodsForRegression    ( void ) { TrainAllMethods(); }

      // testing
      void TestAllMethods();

      // performance evaluation
      void EvaluateAllMethods( void );
      void EvaluateAllVariables(DataLoader *loader, TString options = "" );

      TH1F* EvaluateImportance( DataLoader *loader,VIType vitype, Types::EMVA theMethod,  TString methodTitle, const char *theOption = "" );

      // delete all methods and reset the method vector
      void DeleteAllMethods( void );

      // accessors
      IMethod* GetMethod( const TString& datasetname, const TString& title ) const;
      Bool_t   HasMethod( const TString& datasetname, const TString& title ) const;

      Bool_t Verbose( void ) const { return fVerbose; }
      void SetVerbose( Bool_t v=kTRUE );

      // make ROOT-independent C++ class for classifier response
      // (classifier-specific implementation)
      // If no classifier name is given, help messages for all booked
      // classifiers are printed
      virtual void MakeClass(const TString& datasetname , const TString& methodTitle = "" ) const;

      // prints classifier-specific hepl messages, dedicated to
      // help with the optimisation and configuration options tuning.
      // If no classifier name is given, help messages for all booked
      // classifiers are printed
      void PrintHelpMessage(const TString& datasetname , const TString& methodTitle = "" ) const;

      TDirectory* RootBaseDir() { return (TDirectory*)fgTargetFile; }

      Bool_t IsSilentFile() const { return fSilentFile;}
      Bool_t IsModelPersistence() const { return fModelPersistence; }

      Double_t GetROCIntegral(DataLoader *loader, TString theMethodName, UInt_t iClass = 0,
                              Types::ETreeType type = Types::kTesting);
      Double_t GetROCIntegral(TString datasetname, TString theMethodName, UInt_t iClass = 0,
                              Types::ETreeType type = Types::kTesting);

      // Methods to get a TGraph for an indicated method in dataset.
      // Optional title and axis added with fLegend=kTRUE.
      // Argument iClass used in multiclass settings, otherwise ignored.
      TGraph *GetROCCurve(DataLoader *loader, TString theMethodName, Bool_t setTitles = kTRUE, UInt_t iClass = 0,
                          Types::ETreeType type = Types::kTesting);
      TGraph *GetROCCurve(TString datasetname, TString theMethodName, Bool_t setTitles = kTRUE, UInt_t iClass = 0,
                          Types::ETreeType type = Types::kTesting);

      // Methods to get a TMultiGraph for a given class and all methods in dataset.
      TMultiGraph *GetROCCurveAsMultiGraph(DataLoader *loader, UInt_t iClass, Types::ETreeType type = Types::kTesting);
      TMultiGraph *GetROCCurveAsMultiGraph(TString datasetname, UInt_t iClass, Types::ETreeType type = Types::kTesting);

      // Draw all ROC curves of a given class for all methods in the dataset.
      TCanvas *GetROCCurve(DataLoader *loader, UInt_t iClass = 0, Types::ETreeType type = Types::kTesting);
      TCanvas *GetROCCurve(TString datasetname, UInt_t iClass = 0, Types::ETreeType type = Types::kTesting);

   private:

      // the beautiful greeting message
      void Greetings();

      //evaluate the simple case that is removing 1 variable at time
      TH1F* EvaluateImportanceShort( DataLoader *loader,Types::EMVA theMethod,  TString methodTitle, const char *theOption = "" );
      //evaluate all variables combinations
      TH1F* EvaluateImportanceAll( DataLoader *loader,Types::EMVA theMethod,  TString methodTitle, const char *theOption = "" );
      //evaluate randomly given a number of seeds
      TH1F* EvaluateImportanceRandom( DataLoader *loader,UInt_t nseeds, Types::EMVA theMethod,  TString methodTitle, const char *theOption = "" );

      TH1F* GetImportance(const int nbits,std::vector<Double_t> importances,std::vector<TString> varNames);

      // Helpers for public facing ROC methods
      ROCCurve *GetROC(DataLoader *loader, TString theMethodName, UInt_t iClass = 0,
                       Types::ETreeType type = Types::kTesting);
      ROCCurve *GetROC(TString datasetname, TString theMethodName, UInt_t iClass = 0,
                       Types::ETreeType type = Types::kTesting);

      void WriteDataInformation(DataSetInfo&     fDataSetInfo);

      void SetInputTreesFromEventAssignTrees();

      MethodBase* BookMethodWeightfile(DataLoader *dataloader, TMVA::Types::EMVA methodType, const TString &weightfile);

   private:

      // data members

      TFile*                             fgTargetFile;     //! ROOT output file


      std::vector<TMVA::VariableTransformBase*> fDefaultTrfs;     //! list of transformations on default DataSet

      // cd to local directory
      TString                                   fOptions;         //! option string given by construction (presently only "V")
      TString                                   fTransformations; //! list of transformations to test
      Bool_t                                    fVerbose;         //! verbose mode
      TString                                   fVerboseLevel;    //! verbosity level, controls granularity of logging
      Bool_t                                    fCorrelations;    //! enable to calculate corelations
      Bool_t                                    fROC;             //! enable to calculate ROC values
      Bool_t                                    fSilentFile;      //! used in contructor wihtout file

      TString                                   fJobName;         //! jobname, used as extension in weight file names

      Types::EAnalysisType                      fAnalysisType;    //! the training type
      Bool_t                                    fModelPersistence;//! option to save the trained model in xml file or using serialization


   protected:

      ClassDef(Factory,0);  // The factory creates all MVA methods, and performs their training and testing
   };

} // namespace TMVA

#endif
