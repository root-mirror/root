// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : ResultsClassification                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Base-class for result-vectors                                             *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
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

#ifndef ROOT_TMVA_ResultsClassification
#define ROOT_TMVA_ResultsClassification

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ResultsClassification                                                //
//                                                                      //
// Class that is the base-class for a vector of result                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include "TMVA/Results.h"

namespace TMVA {

   class MsgLogger;

   class ResultsClassification :  public Results {

   public:

      ResultsClassification( const DataSetInfo* dsi, TString resultsName );
      ~ResultsClassification();

      // setters : set  score value and type for each single event.
      // note type=TRUE for signal and FLASE for background
      void     SetValue( Float_t value, Int_t ievt, Bool_t type );

      void     Resize( Int_t entries )   {
         fMvaValues.resize( entries );
         fMvaValuesTypes.resize(entries);
      }
      using TObject::Clear;
      virtual void     Clear(Option_t *)                   { fMvaValues.clear(); fMvaValuesTypes.clear(); }

      // getters
      Long64_t GetSize()                  const { return fMvaValues.size(); }
      virtual const std::vector< Float_t >&  operator [] ( Int_t ievt ) const { fRet[0] = fMvaValues[ievt]; return  fRet; }
      std::vector<Float_t>* GetValueVector()    { return &fMvaValues; }
      std::vector<Bool_t>*  GetValueVectorTypes()    { return &fMvaValuesTypes; }

      Types::EAnalysisType  GetAnalysisType() { return Types::kClassification; }


   private:

      std::vector< Float_t >  fMvaValues;   // mva values (Results)
      std::vector< Bool_t>    fMvaValuesTypes;// mva values type(sig/bkg) (Results)
      mutable std::vector< Float_t >  fRet;         // return val
      mutable MsgLogger*      fLogger;      //! message logger
      MsgLogger& Log() const { return *fLogger; }
   protected:

       ClassDef(Results,2);

   };
}

#endif
