// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard von Toerne, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Reader                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Reader class to be used in the user application to interpret the trained  *
 *      MVAs in an analysis context                                               *
 *                                                                                *
 * Authors (alphabetical order):                                                  *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <peter.speckmayer@cern.ch> - CERN, Switzerland           *
 *      Joerg Stelzer <Joerg.Stelzer@cern.ch>    - CERN, Switzerland              *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Eckhard v. Toerne  <evt@uni-bonn.de>          - U of Bonn, Germany        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_Reader
#define ROOT_TMVA_Reader

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Reader                                                               //
//                                                                      //
// Reader class to be used in the user application to interpret the     //
// trained MVAs in an analysis context                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/Configurable.h"
#include "TMVA/Types.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/DataSetManager.h"

#include <vector>
#include <map>
#include <stdexcept>
#include <string>

namespace TMVA {

   class IMethod;
   class MethodBase;
   class DataSetInfo;
   class MethodCuts;

   class Reader : public Configurable {

   public:

      // without prior specification of variables
      Reader( const TString& theOption="", Bool_t verbose = 0 );

      // STL types
      Reader( std::vector<std::string>& varNames, const TString& theOption = "", Bool_t verbose = 0 );
      Reader( const std::string& varNames, const TString& theOption, Bool_t verbose = 0 );  // format: "var1:var2:..."

      // Root types
      Reader( std::vector<TString>& varNames, const TString& theOption = "", Bool_t verbose = 0 );
      Reader( const TString& varNames, const TString& theOption, Bool_t verbose = 0 );  // format: "var1:var2:..."

      virtual ~Reader( void );

      // book MVA method via weight file
      IMethod* BookMVA( const TString& methodTag, const TString& weightfile );
      IMethod* BookMVA( TMVA::Types::EMVA methodType, const char* xmlstr );
      IMethod* FindMVA( const TString& methodTag );
      // special function for Cuts to avoid dynamic_casts in ROOT macros,
      // which are not properly handled by CINT
      MethodCuts* FindCutsMVA( const TString& methodTag );


      // returns the MVA response for given event
      Double_t EvaluateMVA( const std::vector<Float_t> &, const TString& methodTag, Double_t aux = 0 );
      Double_t EvaluateMVA( const std::vector<Double_t>&, const TString& methodTag, Double_t aux = 0 );
      Double_t EvaluateMVA( MethodBase* method,           Double_t aux = 0 );
      Double_t EvaluateMVA( const TString& methodTag,     Double_t aux = 0 );

      // returns error on MVA response for given event
      // NOTE: must be called AFTER "EvaluateMVA(...)" call !
      Double_t GetMVAError() const { return fMvaEventError; }
      Double_t GetMVAErrorLower() const { return fMvaEventError; }
      Double_t GetMVAErrorUpper() const { return fMvaEventErrorUpper; }

      // regression response
      const std::vector< Float_t >& EvaluateRegression( const TString& methodTag, Double_t aux = 0 );
      const std::vector< Float_t >& EvaluateRegression( MethodBase* method, Double_t aux = 0 );
      Float_t  EvaluateRegression( UInt_t tgtNumber, const TString& methodTag, Double_t aux = 0 );

      // multiclass response
      const std::vector< Float_t >& EvaluateMulticlass( const TString& methodTag, Double_t aux = 0 );
      const std::vector< Float_t >& EvaluateMulticlass( MethodBase* method, Double_t aux = 0 );
      Float_t  EvaluateMulticlass( UInt_t clsNumber, const TString& methodTag, Double_t aux = 0 );

      // probability and rarity accessors (see Users Guide for definition of Rarity)
      Double_t GetProba ( const TString& methodTag, Double_t ap_sig=0.5, Double_t mvaVal=-9999999 );
      Double_t GetRarity( const TString& methodTag, Double_t mvaVal=-9999999 );

      // accessors
      virtual const char* GetName() const { return "Reader"; }
      Bool_t   Verbose( void ) const  { return fVerbose; }
      void     SetVerbose( Bool_t v ) { fVerbose = v; }

      const DataSetInfo& DataInfo() const { return fDataSetInfo; }
      DataSetInfo&       DataInfo()       { return fDataSetInfo; }

      void     AddVariable( const TString& expression, Float_t* );
      void     AddVariable( const TString& expression, Int_t* );

      void     AddSpectator( const TString& expression, Float_t* );
      void     AddSpectator( const TString& expression, Int_t* );

   private:

      DataSetManager* fDataSetManager; // DSMTEST


      TString GetMethodTypeFromFile( const TString& filename );

      // this booking method is internal
      IMethod* BookMVA( Types::EMVA method,  const TString& weightfile );

      DataSetInfo fDataSetInfo; // the data set

      DataInputHandler fDataInputHandler;

      // Init Reader class
      void Init( void );

      // Decode Constructor string (or TString) and fill variable name std::vector
      void DecodeVarNames( const std::string& varNames );
      void DecodeVarNames( const TString& varNames );

      void DeclareOptions();

      Bool_t    fVerbose;            // verbosity
      Bool_t    fSilent;             // silent mode
      Bool_t    fColor;              // color mode
      Bool_t    fCalculateError;     // error calculation mode

      Double_t  fMvaEventError;      // per-event error returned by MVA
      Double_t  fMvaEventErrorUpper; // per-event error returned by MVA

      std::map<TString, IMethod*> fMethodMap; // map of methods

      std::vector<Float_t> fTmpEvalVec; // temporary evaluation vector (if user input is v<double>)

      mutable MsgLogger* fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }

      ClassDef(Reader,0); // Interpret the trained MVAs in an analysis context
   };

}

#endif
