// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Jan Therhaag, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TransformationHandler                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Eckhard v. Toerne  <evt@uni-bonn.de>          - U of Bonn, Germany        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_TransformationHandler
#define ROOT_TMVA_TransformationHandler

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TransformationHandler                                                //
//                                                                      //
// Class that contains all the data information                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TList.h"
#include "TString.h"
#include <vector>

#include "TMVA/DataSetInfo.h"

namespace TMVA {

   class Event;
   class DataSet;
   class Ranking;
   class VariableTransformBase;
   class MsgLogger;

   class TransformationHandler {
   public:

      struct VariableStat {
         Double_t fMean;
         Double_t fRMS;
         Double_t fMin;
         Double_t fMax;
      };

      TransformationHandler( DataSetInfo&, const TString& callerName );
      ~TransformationHandler();

      TString GetName() const;
      TString GetVariableAxisTitle( const VariableInfo& info ) const;

      const Event* Transform(const Event*) const;
      const Event* InverseTransform(const Event*, Bool_t suppressIfNoTargets=true  ) const;

      // overrides the reference classes of all added transformations. Handle with care!!!
      void         SetTransformationReferenceClass( Int_t cls );

      VariableTransformBase* AddTransformation(VariableTransformBase*, Int_t cls );
      const TList& GetTransformationList()   const { return fTransformations; }
      Int_t        GetNumOfTransformations() const { return fTransformations.GetSize(); }
      const std::vector<Event*>* CalcTransformations( const std::vector<Event*>&, Bool_t createNewVector = kFALSE );

      void         CalcStats( const std::vector<Event*>& events );
      void         AddStats ( Int_t k, UInt_t ivar, Double_t mean, Double_t rms, Double_t min, Double_t max );
      Double_t     GetMean  ( Int_t ivar, Int_t cls = -1 ) const;
      Double_t     GetRMS   ( Int_t ivar, Int_t cls = -1 ) const;
      Double_t     GetMin   ( Int_t ivar, Int_t cls = -1 ) const;
      Double_t     GetMax   ( Int_t ivar, Int_t cls = -1 ) const;

      void         WriteToStream ( std::ostream& o ) const;
      void         AddXMLTo      ( void* parent=0 ) const;
      void         ReadFromStream( std::istream& istr );
      void         ReadFromXML   ( void* trfsnode );

      // writer of function code
      void         MakeFunction(std::ostream& fout, const TString& fncName, Int_t part) const;

      // variable ranking
      void         PrintVariableRanking() const;

      // provides string vector giving explicit transformation (only last transform at present -> later show full chain)
      std::vector<TString>* GetTransformationStringsOfLastTransform() const;
      const char*           GetNameOfLastTransform()                  const;

      // modify caller name for output
      void           SetCallerName( const TString& name );
      const TString& GetCallerName() const { return fCallerName; }

      // setting file dir for histograms
      TDirectory*    GetRootDir() const { return fRootBaseDir; }
      void           SetRootDir( TDirectory *d ) { fRootBaseDir = d; }

      void           PlotVariables( const std::vector<Event*>& events, TDirectory* theDirectory = 0 );

   private:

      //      std::vector<TMVA::Event*>* TransformCollection( VariableTransformBase* trf,
      //                                                Int_t cls,
      //                                                std::vector<TMVA::Event*>* events,
      //                                                Bool_t replace ) const;

      const TMVA::VariableInfo& Variable(UInt_t ivar) const { return fDataSetInfo.GetVariableInfos().at(ivar); }
      const TMVA::VariableInfo& Target  (UInt_t itgt) const { return fDataSetInfo.GetTargetInfos()[itgt]; }

      DataSet* Data() { return fDataSetInfo.GetDataSet(); }

      DataSetInfo&          fDataSetInfo;                     // pointer to the datasetinfo
      TList                 fTransformations;                 //! list of transformations
      std::vector< Int_t >  fTransformationsReferenceClasses; //! reference classes for the transformations
      std::vector<std::vector<TMVA::TransformationHandler::VariableStat> >  fVariableStats; // first the variables, then the targets

      Int_t                 fNumC;               // number of categories (#classes +1)

      std::vector<Ranking*> fRanking;            //! ranking object
      TDirectory*           fRootBaseDir;        //! if set put input var hists
      TString               fCallerName;         //! name of the caller for output
      mutable MsgLogger*    fLogger;             //! message logger
      MsgLogger& Log() const { return *fLogger; }
   };
}

#endif
