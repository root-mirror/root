// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Eckhard von Toerne, Jan Therhaag

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

/*! \class TMVA::Reader
\ingroup TMVA

 The Reader class serves to use the MVAs in a specific analysis context.
 Within an event loop, a vector is filled that corresponds to the variables
 that were used to train the MVA(s) during the training stage. This vector
 is transfered to the Reader, who takes care of interpreting the weight
 file of the MVA of choice, and to return the MVA's output. This is then
 used by the user for further analysis.

 Usage:

~~~ {.cpp}
   // ------ before starting the event loop (eg, in the initialisation step)

   //
   // create TMVA::Reader object
   //
   TMVA::Reader *reader = new TMVA::Reader();

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to
   // those given in the weight file(s) that you use
   Float_t var1, var2, var3, var4;
   reader->AddVariable( "var1", &var1 );
   reader->AddVariable( "var2", &var2 );
   reader->AddVariable( "var3", &var3 );
   reader->AddVariable( "var4", &var4 );

   // book the MVA of your choice (prior training of these methods, ie,
   // existence of the weight files is required)
   reader->BookMVA( "Fisher method",  "weights/Fisher.weights.txt"   );
   reader->BookMVA( "MLP method",     "weights/MLP.weights.txt" );
   // ... etc

   // ------- start your event loop

   for (Long64_t ievt=0; ievt<myTree->GetEntries();ievt++) {

     // fill vector with values of variables computed from those in the tree
     var1 = myvar1;
     var2 = myvar2;
     var3 = myvar3;
     var4 = myvar4;

     // retrieve the corresponding MVA output
     double mvaFi = reader->EvaluateMVA( "Fisher method" );
     double mvaNN = reader->EvaluateMVA( "MLP method"    );

     // do something with these ...., e.g., fill them into your ntuple

   } // end of event loop

   delete reader;
~~~
*/

#include "TMVA/Reader.h"

#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TLeaf.h"
#include "TString.h"
#include "TH1D.h"
#include "TVector.h"
#include "TXMLEngine.h"
#include "TMath.h"

#include <cstdlib>

#include <string>
#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Reader::Reader( const TString& theOption, Bool_t verbose )
: Configurable( theOption ),
   fDataSetManager( NULL ), // DSMTEST
   fDataSetInfo(),
   fVerbose( verbose ),
   fSilent ( kFALSE ),
   fColor  ( kFALSE ),
   fCalculateError(kFALSE),
   fMvaEventError( 0 ),
   fMvaEventErrorUpper( 0 ),
   fLogger ( 0 )
{
   fDataSetManager = new DataSetManager( fDataInputHandler );
   fDataSetManager->AddDataSetInfo(fDataSetInfo);
   fLogger = new MsgLogger(this);
   SetConfigName( GetName() );
   DeclareOptions();
   ParseOptions();

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Reader::Reader( std::vector<TString>& inputVars, const TString& theOption, Bool_t verbose )
   : Configurable( theOption ),
     fDataSetManager( NULL ), // DSMTEST
     fDataSetInfo(),
     fVerbose( verbose ),
     fSilent ( kFALSE ),
     fColor  ( kFALSE ),
     fCalculateError(kFALSE),
     fMvaEventError( 0 ),
     fMvaEventErrorUpper( 0 ),   //zjh
     fLogger ( 0 )
{
   fDataSetManager = new DataSetManager( fDataInputHandler );
   fDataSetManager->AddDataSetInfo(fDataSetInfo);
   fLogger = new MsgLogger(this);
   SetConfigName( GetName() );
   DeclareOptions();
   ParseOptions();

   // arguments: names of input variables (vector)
   //            verbose flag
   for (std::vector<TString>::iterator ivar = inputVars.begin(); ivar != inputVars.end(); ++ivar)
      DataInfo().AddVariable( *ivar );

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Reader::Reader( std::vector<std::string>& inputVars, const TString& theOption, Bool_t verbose )
   : Configurable( theOption ),
     fDataSetManager( NULL ), // DSMTEST
     fDataSetInfo(),
     fVerbose( verbose ),
     fSilent ( kFALSE ),
     fColor  ( kFALSE ),
     fCalculateError(kFALSE),
     fMvaEventError( 0 ),
     fMvaEventErrorUpper( 0 ),
     fLogger ( 0 )
{
   fDataSetManager = new DataSetManager( fDataInputHandler );
   fDataSetManager->AddDataSetInfo(fDataSetInfo);
   fLogger = new MsgLogger(this);
   SetConfigName( GetName() );
   DeclareOptions();
   ParseOptions();

   // arguments: names of input variables (vector)
   //            verbose flag
   for (std::vector<std::string>::iterator ivar = inputVars.begin(); ivar != inputVars.end(); ++ivar)
      DataInfo().AddVariable( ivar->c_str() );

   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Reader::Reader( const std::string& varNames, const TString& theOption, Bool_t verbose )
   : Configurable( theOption ),
     fDataSetManager( NULL ), // DSMTEST
     fDataSetInfo(),
     fVerbose( verbose ),
     fSilent ( kFALSE ),
     fColor  ( kFALSE ),
     fCalculateError(kFALSE),
     fMvaEventError( 0 ),
     fMvaEventErrorUpper( 0 ),
     fLogger ( 0 )
{
   fDataSetManager = new DataSetManager( fDataInputHandler );
   fDataSetManager->AddDataSetInfo(fDataSetInfo);
   fLogger = new MsgLogger(this);
   SetConfigName( GetName() );
   DeclareOptions();
   ParseOptions();

   // arguments: names of input variables given in form: "name1:name2:name3"
   //            verbose flag
   DecodeVarNames(varNames);
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Reader::Reader( const TString& varNames, const TString& theOption, Bool_t verbose )
   : Configurable( theOption ),
     fDataSetManager( NULL ), // DSMTEST
     fDataSetInfo(),
     fVerbose( verbose ),
     fSilent ( kFALSE ),
     fColor  ( kFALSE ),
     fCalculateError(kFALSE),
     fMvaEventError( 0 ),
     fMvaEventErrorUpper( 0 ),
     fLogger ( 0 )
{
   fDataSetManager = new DataSetManager( fDataInputHandler );
   fDataSetManager->AddDataSetInfo(fDataSetInfo);
   fLogger = new MsgLogger(this);
   SetConfigName( GetName() );
   DeclareOptions();
   ParseOptions();

   // arguments: names of input variables given in form: "name1:name2:name3"
   //            verbose flag
   DecodeVarNames(varNames);
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// declaration of configuration options

void TMVA::Reader::DeclareOptions()
{
   if (gTools().CheckForSilentOption( GetOptions() )) Log().InhibitOutput(); // make sure is silent if wanted to

   DeclareOptionRef( fVerbose,        "V",      "Verbose flag" );
   DeclareOptionRef( fColor,          "Color",  "Color flag (default True)" );
   DeclareOptionRef( fSilent,         "Silent", "Boolean silent flag (default False)" );
   DeclareOptionRef( fCalculateError, "Error",  "Calculates errors (default False)" );
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Reader::~Reader( void )
{
   delete fDataSetManager; // DSMTEST

   delete fLogger;

   for (auto it=fMethodMap.begin(); it!=fMethodMap.end(); it++){
      MethodBase * kl = dynamic_cast<TMVA::MethodBase*>(it->second);
      delete kl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation (no member variables)

void TMVA::Reader::Init( void )
{
   if (Verbose()) fLogger->SetMinType( kVERBOSE );

   gConfig().SetUseColor( fColor );
   gConfig().SetSilent  ( fSilent );
}

////////////////////////////////////////////////////////////////////////////////
/// Add a float variable or expression to the reader

void TMVA::Reader::AddVariable( const TString& expression, Float_t* datalink )
{
   DataInfo().AddVariable( expression, "", "", 0, 0, 'F', kFALSE ,(void*)datalink ); // <= should this be F or rather T?
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Reader::AddVariable( const TString& expression, Int_t* datalink )
{
   Log() << kFATAL << "Reader::AddVariable( const TString& expression, Int_t* datalink ), this function is deprecated, please provide all variables to the reader as floats" << Endl;
   // Add an integer variable or expression to the reader
   Log() << kFATAL << "Reader::AddVariable( const TString& expression, Int_t* datalink ), this function is deprecated, please provide all variables to the reader as floats" << Endl;
   DataInfo().AddVariable(expression, "", "", 0, 0, 'I', kFALSE, (void*)datalink ); // <= should this be F or rather T?
}

////////////////////////////////////////////////////////////////////////////////
/// Add a float spectator or expression to the reader

void TMVA::Reader::AddSpectator( const TString& expression, Float_t* datalink )
{
   DataInfo().AddSpectator( expression, "", "", 0, 0, 'F', kFALSE ,(void*)datalink );
}

////////////////////////////////////////////////////////////////////////////////
/// Add an integer spectator or expression to the reader

void TMVA::Reader::AddSpectator( const TString& expression, Int_t* datalink )
{
   DataInfo().AddSpectator(expression, "", "", 0, 0, 'I', kFALSE, (void*)datalink );
}

////////////////////////////////////////////////////////////////////////////////
/// read the method type from the file

TString TMVA::Reader::GetMethodTypeFromFile( const TString& filename )
{
   std::ifstream fin( filename );
   if (!fin.good()) { // file not found --> Error
      Log() << kFATAL << "<BookMVA> fatal error: "
            << "unable to open input weight file: " << filename << Endl;
   }

   TString fullMethodName("");
   if (filename.EndsWith(".xml")) {
      fin.close();
      void* doc      = gTools().xmlengine().ParseFile(filename,gTools().xmlenginebuffersize());// the default buffer size in TXMLEngine::ParseFile is 100k. Starting with ROOT 5.29 one can set the buffer size, see: http://savannah.cern.ch/bugs/?78864. This might be necessary for large XML files
      void* rootnode = gTools().xmlengine().DocGetRootElement(doc); // node "MethodSetup"
      gTools().ReadAttr(rootnode, "Method", fullMethodName);
      gTools().xmlengine().FreeDoc(doc);
   }
   else {
      char buf[512];
      fin.getline(buf,512);
      while (!TString(buf).BeginsWith("Method")) fin.getline(buf,512);
      fullMethodName = TString(buf);
      fin.close();
   }
   TString methodType = fullMethodName(0,fullMethodName.Index("::"));
   if (methodType.Contains(" ")) methodType = methodType(methodType.Last(' ')+1,methodType.Length());
   return methodType;
}

////////////////////////////////////////////////////////////////////////////////
/// read method name from weight file

TMVA::IMethod* TMVA::Reader::BookMVA( const TString& methodTag, const TString& weightfile )
{
   // assert non-existence
   if (fMethodMap.find( methodTag ) != fMethodMap.end())
      Log() << kFATAL << "<BookMVA> method tag \"" << methodTag << "\" already exists!" << Endl;

   TString methodType(GetMethodTypeFromFile(weightfile));

   Log() << kINFO << "Booking \"" << methodTag << "\" of type \"" << methodType << "\" from " << weightfile << "." << Endl;

   MethodBase* method = dynamic_cast<MethodBase*>(this->BookMVA( Types::Instance().GetMethodType(methodType),
                                                                 weightfile ) );
   if( method && method->GetMethodType() == Types::kCategory ){
      MethodCategory *methCat = (dynamic_cast<MethodCategory*>(method));
      if( !methCat )
         Log() << kFATAL << "Method with type kCategory cannot be casted to MethodCategory. /Reader" << Endl;
      methCat->fDataSetManager = fDataSetManager;
   }

   return fMethodMap[methodTag] = method;
}

////////////////////////////////////////////////////////////////////////////////
/// books MVA method from weightfile

TMVA::IMethod* TMVA::Reader::BookMVA( TMVA::Types::EMVA methodType, const TString& weightfile )
{
   IMethod *im =
      ClassifierFactory::Instance().Create(Types::Instance().GetMethodName(methodType).Data(), DataInfo(), weightfile);

   MethodBase *method = (dynamic_cast<MethodBase*>(im));

   if (method==0) return im;

   if( method->GetMethodType() == Types::kCategory ){
      MethodCategory *methCat = (dynamic_cast<MethodCategory*>(method));
      if( !methCat )
         Log() << kERROR << "Method with type kCategory cannot be casted to MethodCategory. /Reader" << Endl;
      methCat->fDataSetManager = fDataSetManager;
   }

   method->SetupMethod();

   // when reading older weight files, they could include options
   // that are not supported any longer
   method->DeclareCompatibilityOptions();

   // read weight file
   method->ReadStateFromFile();

   // check for unused options
   method->CheckSetup();

   Log() << kINFO << "Booked classifier \"" << method->GetMethodName()
         << "\" of type: \"" << method->GetMethodTypeName() << "\"" << Endl;

   return method;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::IMethod* TMVA::Reader::BookMVA( TMVA::Types::EMVA methodType, const char* xmlstr )
{
   // books MVA method from weightfile
   IMethod *im =
      ClassifierFactory::Instance().Create(Types::Instance().GetMethodName(methodType).Data(), DataInfo(), "");

   MethodBase *method = (dynamic_cast<MethodBase*>(im));

   if(!method) return 0;

   if( method->GetMethodType() == Types::kCategory ){
      MethodCategory *methCat = (dynamic_cast<MethodCategory*>(method));
      if( !methCat )
         Log() << kFATAL << "Method with type kCategory cannot be casted to MethodCategory. /Reader" << Endl;
      methCat->fDataSetManager = fDataSetManager;
   }

   method->SetupMethod();

   // when reading older weight files, they could include options
   // that are not supported any longer
   method->DeclareCompatibilityOptions();

   // read weight file
   method->ReadStateFromXMLString( xmlstr );

   // check for unused options
   method->CheckSetup();

   Log() << kINFO << "Booked classifier \"" << method->GetMethodName()
         << "\" of type: \"" << method->GetMethodTypeName() << "\"" << Endl;

   return method;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate a std::vector<float> of input data for a given method
/// The parameter aux is obligatory for the cuts method where it represents the efficiency cutoff

Double_t TMVA::Reader::EvaluateMVA( const std::vector<Float_t>& inputVec, const TString& methodTag, Double_t aux )
{
   // create a temporary event from the vector.
   IMethod* imeth = FindMVA( methodTag );
   MethodBase* meth = dynamic_cast<TMVA::MethodBase*>(imeth);
   if(meth==0) return 0;

   //   Event* tmpEvent=new Event(inputVec, 2); // ToDo resolve magic 2 issue
   Event* tmpEvent=new Event(inputVec, DataInfo().GetNVariables()); // is this the solution?
   for (UInt_t i=0; i<inputVec.size(); i++){
      if (TMath::IsNaN(inputVec[i])) {
         Log() << kERROR << i << "-th variable of the event is NaN --> return MVA value -999, \n that's all I can do, please fix or remove this event." << Endl;
         delete tmpEvent;
         return -999;
      }
   }

   if (meth->GetMethodType() == TMVA::Types::kCuts) {
      TMVA::MethodCuts* mc = dynamic_cast<TMVA::MethodCuts*>(meth);
      if(mc)
         mc->SetTestSignalEfficiency( aux );
   }
   Double_t val = meth->GetMvaValue( tmpEvent, (fCalculateError?&fMvaEventError:0));
   delete tmpEvent;
   return val;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate a std::vector<double> of input data for a given method
/// The parameter aux is obligatory for the cuts method where it represents the efficiency cutoff

Double_t TMVA::Reader::EvaluateMVA( const std::vector<Double_t>& inputVec, const TString& methodTag, Double_t aux )
{
   // performs a copy to float values which are internally used by all methods
   if(fTmpEvalVec.size() != inputVec.size())
      fTmpEvalVec.resize(inputVec.size());

   for (UInt_t idx=0; idx!=inputVec.size(); idx++ )
      fTmpEvalVec[idx]=inputVec[idx];

   return EvaluateMVA( fTmpEvalVec, methodTag, aux );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates MVA for given set of input variables

Double_t TMVA::Reader::EvaluateMVA( const TString& methodTag, Double_t aux )
{
   IMethod* method = 0;

   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it == fMethodMap.end()) {
      Log() << kINFO << "<EvaluateMVA> unknown classifier in map; "
            << "you looked for \"" << methodTag << "\" within available methods: " << Endl;
      for (it = fMethodMap.begin(); it!=fMethodMap.end(); ++it) Log() << "--> " << it->first << Endl;
      Log() << "Check calling string" << kFATAL << Endl;
   }

   else method = it->second;

   MethodBase * kl = dynamic_cast<TMVA::MethodBase*>(method);

   if(kl==0)
      Log() << kFATAL << methodTag << " is not a method" << Endl;

   // check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
   // it is not again checked in each of these subsequent calls..
   const Event* ev = kl->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN --> return MVA value -999, \n that's all I can do, please fix or remove this event." << Endl;
         return -999;
      }
   }
   return this->EvaluateMVA( kl, aux );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates the MVA

Double_t TMVA::Reader::EvaluateMVA( MethodBase* method, Double_t aux )
{
   // the aux value is only needed for MethodCuts: it sets the
   // required signal efficiency
   if (method->GetMethodType() == TMVA::Types::kCuts) {
      TMVA::MethodCuts* mc = dynamic_cast<TMVA::MethodCuts*>(method);
      if(mc)
         mc->SetTestSignalEfficiency( aux );
   }

   return method->GetMvaValue( (fCalculateError?&fMvaEventError:0),
                               (fCalculateError?&fMvaEventErrorUpper:0) );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates MVA for given set of input variables

const std::vector< Float_t >& TMVA::Reader::EvaluateRegression( const TString& methodTag, Double_t aux )
{
   IMethod* method = 0;

   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it == fMethodMap.end()) {
      Log() << kINFO << "<EvaluateMVA> unknown method in map; "
            << "you looked for \"" << methodTag << "\" within available methods: " << Endl;
      for (it = fMethodMap.begin(); it!=fMethodMap.end(); ++it) Log() << "--> " << it->first << Endl;
      Log() << "Check calling string" << kFATAL << Endl;
   }
   else method = it->second;

   MethodBase * kl = dynamic_cast<TMVA::MethodBase*>(method);

   if(kl==0)
      Log() << kFATAL << methodTag << " is not a method" << Endl;
   // check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
   // it is not again checked in each of these subsequent calls..
   const Event* ev = kl->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN, \n regression values might evaluate to .. what do I know. \n sorry this warning is all I can do, please fix or remove this event." << Endl;
      }
   }

   return this->EvaluateRegression( kl, aux );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates the regression MVA
/// check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
/// it is not again checked in each of these subsequent calls.

const std::vector< Float_t >& TMVA::Reader::EvaluateRegression( MethodBase* method, Double_t /*aux*/ )
{
   const Event* ev = method->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN, \n regression values might evaluate to .. what do I know. \n sorry this warning is all I can do, please fix or remove this event." << Endl;
      }
   }
   return method->GetRegressionValues();
}


////////////////////////////////////////////////////////////////////////////////
/// evaluates the regression MVA

Float_t TMVA::Reader::EvaluateRegression( UInt_t tgtNumber, const TString& methodTag, Double_t aux )
{
   try {
      return EvaluateRegression(methodTag, aux).at(tgtNumber);
   }
   catch (std::out_of_range &) {
      Log() << kWARNING << "Regression could not be evaluated for target-number " << tgtNumber << Endl;
      return 0;
   }
}



////////////////////////////////////////////////////////////////////////////////
/// evaluates MVA for given set of input variables

const std::vector< Float_t >& TMVA::Reader::EvaluateMulticlass( const TString& methodTag, Double_t aux )
{
   IMethod* method = 0;

   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it == fMethodMap.end()) {
      Log() << kINFO << "<EvaluateMVA> unknown method in map; "
            << "you looked for \"" << methodTag << "\" within available methods: " << Endl;
      for (it = fMethodMap.begin(); it!=fMethodMap.end(); ++it) Log() << "--> " << it->first << Endl;
      Log() << "Check calling string" << kFATAL << Endl;
   }
   else method = it->second;

   MethodBase * kl = dynamic_cast<TMVA::MethodBase*>(method);

   if(kl==0)
      Log() << kFATAL << methodTag << " is not a method" << Endl;
   // check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
   // it is not again checked in each of these subsequent calls..

   const Event* ev = kl->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN, \n regression values might evaluate to .. what do I know. \n sorry this warning is all I can do, please fix or remove this event." << Endl;
      }
   }

   return this->EvaluateMulticlass( kl, aux );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates the multiclass MVA
/// check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
/// it is not again checked in each of these subsequent calls.

const std::vector< Float_t >& TMVA::Reader::EvaluateMulticlass( MethodBase* method, Double_t /*aux*/ )
{
   const Event* ev = method->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN, \n regression values might evaluate to .. what do I know. \n sorry this warning is all I can do, please fix or remove this event." << Endl;
      }
   }
   return method->GetMulticlassValues();
}


////////////////////////////////////////////////////////////////////////////////
/// evaluates the multiclass MVA

Float_t TMVA::Reader::EvaluateMulticlass( UInt_t clsNumber, const TString& methodTag, Double_t aux )
{
   try {
      return EvaluateMulticlass(methodTag, aux).at(clsNumber);
   }
   catch (std::out_of_range &) {
      Log() << kWARNING << "Multiclass could not be evaluated for class-number " << clsNumber << Endl;
      return 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// return pointer to method with tag "methodTag"

TMVA::IMethod* TMVA::Reader::FindMVA( const TString& methodTag )
{
   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it != fMethodMap.end()) return it->second;
   Log() << kERROR << "Method " << methodTag << " not found!" << Endl;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// special function for Cuts to avoid dynamic_casts in ROOT macros,
/// which are not properly handled by CINT

TMVA::MethodCuts* TMVA::Reader::FindCutsMVA( const TString& methodTag )
{
   return dynamic_cast<MethodCuts*>(FindMVA(methodTag));
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates probability of MVA for given set of input variables

Double_t TMVA::Reader::GetProba( const TString& methodTag,  Double_t ap_sig, Double_t mvaVal )
{
   IMethod* method = 0;
   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it == fMethodMap.end()) {
      for (it = fMethodMap.begin(); it!=fMethodMap.end(); ++it) Log() << "M" << it->first << Endl;
      Log() << kFATAL << "<EvaluateMVA> unknown classifier in map: " << method << "; "
            << "you looked for " << methodTag<< " while the available methods are : " << Endl;
   }
   else method = it->second;

   MethodBase* kl = dynamic_cast<MethodBase*>(method);
   if(kl==0) return -1;
   // check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
   // it is not again checked in each of these subsequent calls..
   const Event* ev = kl->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN --> return MVA value -999, \n that's all I can do, please fix or remove this event." << Endl;
         return -999;
      }
   }

   if (mvaVal == -9999999) mvaVal = kl->GetMvaValue();

   return kl->GetProba( mvaVal, ap_sig );
}

////////////////////////////////////////////////////////////////////////////////
/// evaluates the MVA's rarity

Double_t TMVA::Reader::GetRarity( const TString& methodTag, Double_t mvaVal )
{
   IMethod* method = 0;
   std::map<TString, IMethod*>::iterator it = fMethodMap.find( methodTag );
   if (it == fMethodMap.end()) {
      for (it = fMethodMap.begin(); it!=fMethodMap.end(); ++it) Log() << "M" << it->first << Endl;
      Log() << kFATAL << "<EvaluateMVA> unknown classifier in map: \"" << method << "\"; "
            << "you looked for \"" << methodTag<< "\" while the available methods are : " << Endl;
   }
   else method = it->second;

   MethodBase* kl = dynamic_cast<MethodBase*>(method);
   if(kl==0) return -1;
   // check for NaN in event data:  (note: in the factory, this check was done already at the creation of the datasets, hence
   // it is not again checked in each of these subsequent calls..
   const Event* ev = kl->GetEvent();
   for (UInt_t i=0; i<ev->GetNVariables(); i++){
      if (TMath::IsNaN(ev->GetValue(i))) {
         Log() << kERROR << i << "-th variable of the event is NaN --> return MVA value -999, \n that's all I can do, please fix or remove this event." << Endl;
         return -999;
      }
   }

   if (mvaVal == -9999999) mvaVal = kl->GetMvaValue();

   return kl->GetRarity( mvaVal );
}

// ---------------------------------------------------------------------------------------
// ----- methods related to the decoding of the input variable names ---------------------
// ---------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// decodes "name1:name2:..." form

void TMVA::Reader::DecodeVarNames( const std::string& varNames )
{
   size_t ipos = 0, f = 0;
   while (f != varNames.length()) {
      f = varNames.find( ':', ipos );
      if (f > varNames.length()) f = varNames.length();
      std::string subs = varNames.substr( ipos, f-ipos ); ipos = f+1;
      DataInfo().AddVariable( subs.c_str() );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// decodes "name1:name2:..." form

void TMVA::Reader::DecodeVarNames( const TString& varNames )
{
   TString format;
   Int_t   n = varNames.Length();
   TString format_obj;

   for (int i=0; i< n+1 ; i++) {
      format.Append(varNames(i));
      if (varNames(i) == ':' || i == n) {
         format.Chop();
         format_obj = format;
         format_obj.ReplaceAll("@","");
         DataInfo().AddVariable( format_obj );
         format.Resize(0);
      }
   }
}
