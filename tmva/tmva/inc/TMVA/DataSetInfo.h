// // @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataSetInfo                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <speckmay@mail.cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - DESY, Germany                  *
 *                                                                                *
 * Copyright (c) 2008-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      DESY Hamburg, Germany                                                     *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_DataSetInfo
#define ROOT_TMVA_DataSetInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DataSetInfo                                                          //
//                                                                      //
// Class that contains all the data information                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>
#include <vector>
#include <map>

#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TCut.h"
#include "TMatrixDfwd.h"

#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/ClassInfo.h"
#include "TMVA/Event.h"

class TH2;

namespace TMVA {

   class DataSet;
   class VariableTransformBase;
   class MsgLogger;
   class DataSetManager;

   class DataSetInfo : public TObject {

   public:

      enum { kIsArrayVariable = BIT(15) };

      DataSetInfo(const TString& name = "Default");
      virtual ~DataSetInfo();

      virtual const char* GetName() const { return fName.Data(); }

      // the data set
      void        ClearDataSet() const;
      DataSet*    GetDataSet() const;

      // ---
      // the variable data
      // ---
      VariableInfo&     AddVariable( const TString& expression, const TString& title = "", const TString& unit = "", 
                                     Double_t min = 0, Double_t max = 0, char varType='F', 
                                     Bool_t normalized = kTRUE, void* external = 0 );
      VariableInfo&     AddVariable( const VariableInfo& varInfo );

      // NEW: add an array of variables (e.g. for image data)
      void AddVariablesArray(const TString &expression, Int_t size, const TString &title = "", const TString &unit = "",
                             Double_t min = 0, Double_t max = 0, char type = 'F', Bool_t normalized = kTRUE,
                             void *external = 0 );

      VariableInfo&     AddTarget  ( const TString& expression, const TString& title, const TString& unit, 
                                     Double_t min, Double_t max, Bool_t normalized = kTRUE, void* external = 0 );
      VariableInfo&     AddTarget  ( const VariableInfo& varInfo );

      VariableInfo&     AddSpectator ( const TString& expression, const TString& title, const TString& unit, 
                                       Double_t min, Double_t max, char type = 'F', Bool_t normalized = kTRUE, void* external = 0 );
      VariableInfo&     AddSpectator ( const VariableInfo& varInfo );

      ClassInfo*        AddClass   ( const TString& className );

      // accessors

      // general
      std::vector<VariableInfo>&       GetVariableInfos()         { return fVariables; }
      const std::vector<VariableInfo>& GetVariableInfos() const   { return fVariables; }
      VariableInfo&                    GetVariableInfo( Int_t i ) { return fVariables.at(i); }
      const VariableInfo&              GetVariableInfo( Int_t i ) const { return fVariables.at(i); }

      Int_t GetVarArraySize(const TString &expression) const { 
         auto element = fVarArrays.find(expression);
         return (element != fVarArrays.end()) ? element->second : -1;
       }
       Bool_t IsVariableFromArray(Int_t i) const { return GetVariableInfo(i).TestBit(DataSetInfo::kIsArrayVariable);  }

       std::vector<VariableInfo> &GetTargetInfos()
       {
          return fTargets;
       }
       const std::vector<VariableInfo> &GetTargetInfos() const { return fTargets; }
       VariableInfo &GetTargetInfo(Int_t i) { return fTargets.at(i); }
       const VariableInfo &GetTargetInfo(Int_t i) const { return fTargets.at(i); }

       std::vector<VariableInfo> &GetSpectatorInfos() { return fSpectators; }
       const std::vector<VariableInfo> &GetSpectatorInfos() const { return fSpectators; }
       VariableInfo &GetSpectatorInfo(Int_t i) { return fSpectators.at(i); }
       const VariableInfo &GetSpectatorInfo(Int_t i) const { return fSpectators.at(i); }

       UInt_t GetNVariables() const { return fVariables.size(); }
       UInt_t GetNTargets() const { return fTargets.size(); }
       UInt_t GetNSpectators(bool all = kTRUE) const;

       const TString &GetNormalization() const { return fNormalization; }
       void SetNormalization(const TString &norm) { fNormalization = norm; }

       void SetTrainingSumSignalWeights(Double_t trainingSumSignalWeights)
       {
          fTrainingSumSignalWeights = trainingSumSignalWeights;}
      void SetTrainingSumBackgrWeights(Double_t trainingSumBackgrWeights){fTrainingSumBackgrWeights = trainingSumBackgrWeights;}
      void SetTestingSumSignalWeights (Double_t testingSumSignalWeights ){fTestingSumSignalWeights  = testingSumSignalWeights ;}
      void SetTestingSumBackgrWeights (Double_t testingSumBackgrWeights ){fTestingSumBackgrWeights  = testingSumBackgrWeights ;}

      Double_t GetTrainingSumSignalWeights();
      Double_t GetTrainingSumBackgrWeights();
      Double_t GetTestingSumSignalWeights ();
      Double_t GetTestingSumBackgrWeights ();



      // classification information
      Int_t              GetClassNameMaxLength() const;
      Int_t              GetVariableNameMaxLength() const;
      Int_t              GetTargetNameMaxLength() const;
      ClassInfo*         GetClassInfo( Int_t clNum ) const;
      ClassInfo*         GetClassInfo( const TString& name ) const;
      void               PrintClasses() const;
      UInt_t             GetNClasses() const { return fClasses.size(); }
      Bool_t             IsSignal( const Event* ev ) const;
      std::vector<Float_t>* GetTargetsForMulticlass( const Event* ev );
      UInt_t             GetSignalClassIndex(){return fSignalClass;}

      // by variable
      Int_t              FindVarIndex( const TString& )      const;

      // weights
      const TString      GetWeightExpression(Int_t i)      const { return GetClassInfo(i)->GetWeight(); }
      void               SetWeightExpression( const TString& exp, const TString& className = "" );

      // cuts
      const TCut&        GetCut (Int_t i)                         const { return GetClassInfo(i)->GetCut(); }
      const TCut&        GetCut ( const TString& className )      const { return GetClassInfo(className)->GetCut(); }
      void               SetCut ( const TCut& cut, const TString& className );
      void               AddCut ( const TCut& cut, const TString& className );
      Bool_t             HasCuts() const;

      std::vector<TString> GetListOfVariables() const;

      // correlation matrix
      const TMatrixD*    CorrelationMatrix     ( const TString& className ) const;
      void               SetCorrelationMatrix  ( const TString& className, TMatrixD* matrix );
      void               PrintCorrelationMatrix( const TString& className );
      TH2*               CreateCorrelationMatrixHist( const TMatrixD* m,
                                                      const TString& hName,
                                                      const TString& hTitle ) const;

      // options
      void               SetSplitOptions(const TString& so) { fSplitOptions = so; fNeedsRebuilding = kTRUE; }
      const TString&     GetSplitOptions() const { return fSplitOptions; }

      // root dir
      void               SetRootDir(TDirectory* d) { fOwnRootDir = d; }
      TDirectory*        GetRootDir() const { return fOwnRootDir; }

      void               SetMsgType( EMsgType t ) const;

      DataSetManager*   GetDataSetManager(){return fDataSetManager;}
   private:

      TMVA::DataSetManager*            fDataSetManager; // DSMTEST
      void                       SetDataSetManager( DataSetManager* dsm ) { fDataSetManager = dsm; } // DSMTEST
      friend class DataSetManager;  // DSMTEST (datasetmanager test)

   DataSetInfo( const DataSetInfo& ) : TObject() {}

      void PrintCorrelationMatrix( TTree* theTree );

      TString                    fName;              // name of the dataset info object

      mutable DataSet*           fDataSet;           // dataset, owned by this datasetinfo object
      mutable Bool_t             fNeedsRebuilding;   // flag if rebuilding of dataset is needed (after change of cuts, vars, etc.)

      // expressions/formulas
      std::vector<VariableInfo>  fVariables;         // list of variable expressions/internal names
      std::vector<VariableInfo>  fTargets;           // list of targets expressions/internal names
      std::vector<VariableInfo>  fSpectators;        // list of spectators expressions/internal names

      // variable arrays
      std::map<TString, int> fVarArrays;

      // the classes
      mutable std::vector<ClassInfo*> fClasses;      // name and other infos of the classes

      TString                    fNormalization;     //
      TString                    fSplitOptions;      //

      Double_t                   fTrainingSumSignalWeights;
      Double_t                   fTrainingSumBackgrWeights;
      Double_t                   fTestingSumSignalWeights ;
      Double_t                   fTestingSumBackgrWeights ;


      
      TDirectory*                fOwnRootDir;        // ROOT output dir
      Bool_t                     fVerbose;           // Verbosity

      UInt_t                     fSignalClass;       // index of the class with the name signal

      std::vector<Float_t>*      fTargetsForMulticlass;//-> all targets 0 except the one with index==classNumber
      
      mutable MsgLogger*         fLogger;            //! message logger
      MsgLogger& Log() const { return *fLogger; }

   public:
       
       ClassDef(DataSetInfo,1);
   };
}

#endif
