// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableTransformBase                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::VariableTransformBase
\ingroup TMVA
Linear interpolation class.
*/

#include "TMVA/VariableTransformBase.h"

#include "TMVA/Config.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Ranking.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/Version.h"

#include "THashTable.h"
#include "TList.h"
#include "TObjString.h"
#include "TMath.h"
#include "TVectorD.h"

#include <algorithm>
#include <cassert>
#include <exception>
#include <iomanip>
#include <stdexcept>
#include <set>

bool advanced;

ClassImp(TMVA::VariableTransformBase);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

TMVA::VariableTransformBase::VariableTransformBase( DataSetInfo& dsi,
                                                    Types::EVariableTransform tf,
                                                    const TString& trfName )
: TObject(),
   fDsi(dsi),
   fDsiOutput(NULL),
   fTransformedEvent(0),
   fBackTransformedEvent(0),
   fVariableTransform(tf),
   fEnabled( kTRUE ),
   fCreated( kFALSE ),
   fNormalise( kFALSE ),
   fTransformName(trfName),
   fVariableTypesAreCounted(false),
   fNVariables(0),
   fNTargets(0),
   fNSpectators(0),
   fSortGet(kTRUE),
   fTMVAVersion(TMVA_VERSION_CODE),
   fLogger( 0 )
{
   fLogger = new MsgLogger(this, kINFO);
   for (UInt_t ivar = 0; ivar < fDsi.GetNVariables(); ivar++) {
      fVariables.push_back( VariableInfo( fDsi.GetVariableInfo(ivar) ) );
   }
   for (UInt_t itgt = 0; itgt < fDsi.GetNTargets(); itgt++) {
      fTargets.push_back( VariableInfo( fDsi.GetTargetInfo(itgt) ) );
   }
   for (UInt_t ispct = 0; ispct < fDsi.GetNSpectators(); ispct++) {
      fSpectators.push_back( VariableInfo( fDsi.GetSpectatorInfo(ispct) ) );
   }
}

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableTransformBase::~VariableTransformBase()
{
   if (fTransformedEvent!=0)     delete fTransformedEvent;
   if (fBackTransformedEvent!=0) delete fBackTransformedEvent;
   // destructor
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// select the variables/targets/spectators which serve as input to the transformation

void TMVA::VariableTransformBase::SelectInput( const TString& _inputVariables, Bool_t putIntoVariables  )
{
   TString inputVariables = _inputVariables;

   // unselect all variables first
   fGet.clear();

   UInt_t nvars  = GetNVariables();
   UInt_t ntgts  = GetNTargets();
   UInt_t nspcts = GetNSpectators();

   typedef std::set<Int_t> SelectedIndices;

   SelectedIndices varIndices;
   SelectedIndices tgtIndices;
   SelectedIndices spctIndices;

   if (inputVariables == "") // default is all variables and all targets
      {                         //   (the default can be changed by decorating this member function in the implementations)
         inputVariables = "_V_,_T_";
      }

   TList* inList = gTools().ParseFormatLine( inputVariables, "," );
   TListIter inIt(inList);
   while (TObjString* os = (TObjString*)inIt()) {

      TString variables = os->GetString();

      if( variables.BeginsWith("_") && variables.EndsWith("_") ) { // special symbol (keyword)
         variables.Remove( 0,1); // remove first "_"
         variables.Remove( variables.Length()-1,1 ); // remove last "_"

         if( variables.BeginsWith("V") ) {       // variables
            variables.Remove(0,1); // remove "V"
            if( variables.Length() == 0 ){
               for( UInt_t ivar = 0; ivar < nvars; ++ivar ) {
                  fGet.push_back( std::pair<Char_t,UInt_t>('v',ivar) );
                  varIndices.insert( ivar );
               }
            } else {
               UInt_t idx = variables.Atoi();
               if( idx >= nvars )
                  Log() << kFATAL << "You selected variable with index : " << idx << " of only " << nvars << " variables." << Endl;
               fGet.push_back( std::pair<Char_t,UInt_t>('v',idx) );
               varIndices.insert( idx );
            }
         }else if( variables.BeginsWith("T") ) {       // targets
            variables.Remove(0,1); // remove "T"
            if( variables.Length() == 0 ){
               for( UInt_t itgt = 0; itgt < ntgts; ++itgt ) {
                  fGet.push_back( std::pair<Char_t,UInt_t>('t',itgt) );
                  tgtIndices.insert( itgt );
               }
            } else {
               UInt_t idx = variables.Atoi();
               if( idx >= ntgts )
                  Log() << kFATAL << "You selected target with index : " << idx << " of only " << ntgts << " targets." << Endl;
               fGet.push_back( std::pair<Char_t,UInt_t>('t',idx) );
               tgtIndices.insert( idx );
            }
         }else if( variables.BeginsWith("S") ) {       // spectators
            variables.Remove(0,1); // remove "S"
            if( variables.Length() == 0 ){
               for( UInt_t ispct = 0; ispct < nspcts; ++ispct ) {
                  fGet.push_back( std::pair<Char_t,UInt_t>('s',ispct) );
                  spctIndices.insert( ispct );
               }
            } else {
               UInt_t idx = variables.Atoi();
               if( idx >= nspcts )
                  Log() << kFATAL << "You selected spectator with index : " << idx << " of only " << nspcts << " spectators." << Endl;
               fGet.push_back( std::pair<Char_t,UInt_t>('s',idx) );
               spctIndices.insert( idx );
            }
         }else if( TString("REARRANGE").BeginsWith(variables) ) {       // toggle rearrange sorting (take sort order given in the options)
            ToggleInputSortOrder( kFALSE );
            if( !fSortGet )
               Log() << kINFO << "Variable rearrangement set true: Variable order given in transformation option is used for input to transformation!" << Endl;

         }
      }else{ // no keyword, ... user provided variable labels
         Int_t numIndices = varIndices.size()+tgtIndices.size()+spctIndices.size();
         for( UInt_t ivar = 0; ivar < nvars; ++ivar ) { // search all variables
            if( fDsi.GetVariableInfo( ivar ).GetLabel() == variables ) {
               fGet.push_back( std::pair<Char_t,UInt_t>('v',ivar) );
               varIndices.insert( ivar );
               break;
            }
         }
         for( UInt_t itgt = 0; itgt < ntgts; ++itgt ) { // search all targets
            if( fDsi.GetTargetInfo( itgt ).GetLabel() == variables ) {
               fGet.push_back( std::pair<Char_t,UInt_t>('t',itgt) );
               tgtIndices.insert( itgt );
               break;
            }
         }
         for( UInt_t ispct = 0; ispct < nspcts; ++ispct ) { // search all spectators
            if( fDsi.GetSpectatorInfo( ispct ).GetLabel() == variables ) {
               fGet.push_back( std::pair<Char_t,UInt_t>('s',ispct) );
               spctIndices.insert( ispct );
               break;
            }
         }
         Int_t numIndicesEndOfLoop = varIndices.size()+tgtIndices.size()+spctIndices.size();
         if( numIndicesEndOfLoop == numIndices )
            Log() << kWARNING << "Error at parsing the options for the variable transformations: Variable/Target/Spectator '" << variables.Data() << "' not found." << Endl;
         numIndices = numIndicesEndOfLoop;
      }
   }


   if( putIntoVariables ) {
      Int_t idx = 0;
      for( SelectedIndices::iterator it = varIndices.begin(), itEnd = varIndices.end(); it != itEnd; ++it ) {
         fPut.push_back( std::pair<Char_t,UInt_t>('v',idx) );
         ++idx;
      }
      for( SelectedIndices::iterator it = tgtIndices.begin(), itEnd = tgtIndices.end(); it != itEnd; ++it ) {
         fPut.push_back( std::pair<Char_t,UInt_t>('t',idx) );
         ++idx;
      }
      for( SelectedIndices::iterator it = spctIndices.begin(), itEnd = spctIndices.end(); it != itEnd; ++it ) {
         fPut.push_back( std::pair<Char_t,UInt_t>('s',idx) );
         ++idx;
      }
   }else {
      for( SelectedIndices::iterator it = varIndices.begin(), itEnd = varIndices.end(); it != itEnd; ++it ) {
         Int_t idx = (*it);
         fPut.push_back( std::pair<Char_t,UInt_t>('v',idx) );
      }
      for( SelectedIndices::iterator it = tgtIndices.begin(), itEnd = tgtIndices.end(); it != itEnd; ++it ) {
         Int_t idx = (*it);
         fPut.push_back( std::pair<Char_t,UInt_t>('t',idx) );
      }
      for( SelectedIndices::iterator it = spctIndices.begin(), itEnd = spctIndices.end(); it != itEnd; ++it ) {
         Int_t idx = (*it);
         fPut.push_back( std::pair<Char_t,UInt_t>('s',idx) );
      }

      // if sorting is turned on, fGet should have the indices sorted as fPut has them.
      if( fSortGet ) {
         fGet.clear();
         fGet.assign( fPut.begin(), fPut.end() );
      }
   }

   Log() << kHEADER << "Transformation, Variable selection : " << Endl;

   // choose the new dsi for output if present, if not, take the common one
   const DataSetInfo* outputDsiPtr = (fDsiOutput? &(*fDsiOutput) : &fDsi );



   ItVarTypeIdx itGet = fGet.begin(), itGetEnd = fGet.end();
   ItVarTypeIdx itPut = fPut.begin(); // , itPutEnd = fPut.end();
   for( ; itGet != itGetEnd; ++itGet ) {
      TString inputTypeString  = "?";

      Char_t inputType = (*itGet).first;
      Int_t inputIdx  = (*itGet).second;

      TString inputLabel = "NOT FOND";
      if( inputType == 'v' ) {
         inputLabel = fDsi.GetVariableInfo( inputIdx ).GetLabel();
         inputTypeString = "variable";
      }
      else if( inputType == 't' ){
         inputLabel = fDsi.GetTargetInfo( inputIdx ).GetLabel();
         inputTypeString = "target";
      }
      else if( inputType == 's' ){
         inputLabel = fDsi.GetSpectatorInfo( inputIdx ).GetLabel();
         inputTypeString = "spectator";
      }

      TString outputTypeString  = "?";

      Char_t outputType = (*itPut).first;
      Int_t outputIdx  = (*itPut).second;

      TString outputLabel = "NOT FOUND";
      if( outputType == 'v' ) {
         outputLabel = outputDsiPtr->GetVariableInfo( outputIdx ).GetLabel();
         outputTypeString = "variable";
      }
      else if( outputType == 't' ){
         outputLabel = outputDsiPtr->GetTargetInfo( outputIdx ).GetLabel();
         outputTypeString = "target";
      }
      else if( outputType == 's' ){
         outputLabel = outputDsiPtr->GetSpectatorInfo( outputIdx ).GetLabel();
         outputTypeString = "spectator";
      }
      Log() << kINFO << "Input : " << inputTypeString.Data() << " '" << inputLabel.Data() << "'" << " <---> " << "Output : " << outputTypeString.Data() << " '" << outputLabel.Data() << "'" << Endl;
      Log() << kDEBUG << "\t(index=" << inputIdx << ")." << "\t(index=" << outputIdx << ")." << Endl;

      ++itPut;
   }
   //   Log() << kINFO << Endl;
}


////////////////////////////////////////////////////////////////////////////////
/// select the values from the event

Bool_t TMVA::VariableTransformBase::GetInput( const Event* event, std::vector<Float_t>& input, std::vector<Char_t>& mask, Bool_t backTransformation ) const
{
   ItVarTypeIdxConst itEntry;
   ItVarTypeIdxConst itEntryEnd;

   input.clear();
   mask.clear();

   if( backTransformation && !fPut.empty() ){
      itEntry = fPut.begin();
      itEntryEnd = fPut.end();
      input.reserve(fPut.size());
   }
   else {
      itEntry = fGet.begin();
      itEntryEnd = fGet.end();
      input.reserve(fGet.size() );
   }

   Bool_t hasMaskedEntries = kFALSE;
   //   event->Print(std::cout);
   for( ; itEntry != itEntryEnd; ++itEntry ) {
      Char_t type = (*itEntry).first;
      Int_t  idx  = (*itEntry).second;

      try{
         switch( type ) {
         case 'v':
            input.push_back( event->GetValue(idx) );
            break;
         case 't':
            input.push_back( event->GetTarget(idx) );
            break;
         case 's':
            input.push_back( event->GetSpectator(idx) );
            break;
         default:
            Log() << kFATAL << "VariableTransformBase/GetInput : unknown type '" << type << "'." << Endl;
         }
         mask.push_back(kFALSE);
      }
      catch(std::out_of_range& /* excpt */ ){ // happens when an event is transformed which does not yet have the targets calculated (in the application phase)
         input.push_back(0.f);
         mask.push_back(kTRUE);
         hasMaskedEntries = kTRUE;
      }
   }
   return hasMaskedEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// select the values from the event

void TMVA::VariableTransformBase::SetOutput( Event* event, std::vector<Float_t>& output, std::vector<Char_t>& mask, const Event* oldEvent, Bool_t backTransformation ) const
{
   std::vector<Float_t>::iterator itOutput = output.begin();
   std::vector<Char_t>::iterator  itMask   = mask.begin();

   if( oldEvent )
      event->CopyVarValues( *oldEvent );

   try {

      ItVarTypeIdxConst itEntry;
      ItVarTypeIdxConst itEntryEnd;

      if( backTransformation || fPut.empty() ){ // as in GetInput, but the other way round (from fPut for transformation, from fGet for backTransformation)
         itEntry = fGet.begin();
         itEntryEnd = fGet.end();
      }
      else {
         itEntry = fPut.begin();
         itEntryEnd = fPut.end();
      }


      for( ; itEntry != itEntryEnd; ++itEntry ) {

         if( (*itMask) ){ // if the value is masked
            continue;
         }

         Char_t type = (*itEntry).first;
         Int_t  idx  = (*itEntry).second;
         if (itOutput == output.end()) Log() << kFATAL << "Read beyond array boundaries in VariableTransformBase::SetOutput"<<Endl;
         Float_t value = (*itOutput);

         switch( type ) {
         case 'v':
            event->SetVal( idx, value );
            break;
         case 't':
            event->SetTarget( idx, value );
            break;
         case 's':
            event->SetSpectator( idx, value );
            break;
         default:
            Log() << kFATAL << "VariableTransformBase/GetInput : unknown type '" << type << "'." << Endl;
         }
         if( !(*itMask) ) ++itOutput;
         ++itMask;

      }
   }catch( std::exception& except ){
      Log() << kFATAL << "VariableTransformBase/SetOutput : exception/" << except.what() << Endl;
      throw;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// count variables, targets and spectators

void TMVA::VariableTransformBase::CountVariableTypes( UInt_t& nvars, UInt_t& ntgts, UInt_t& nspcts ) const
{
   if( fVariableTypesAreCounted ){
      nvars = fNVariables;
      ntgts = fNTargets;
      nspcts = fNSpectators;
      return;
   }

   nvars = ntgts = nspcts = 0;

   for( ItVarTypeIdxConst itEntry = fGet.begin(), itEntryEnd = fGet.end(); itEntry != itEntryEnd; ++itEntry ) {
      Char_t type = (*itEntry).first;

      switch( type ) {
      case 'v':
         nvars++;
         break;
      case 't':
         ntgts++;
         break;
      case 's':
         nspcts++;
         break;
      default:
         Log() << kFATAL << "VariableTransformBase/GetVariableTypeNumbers : unknown type '" << type << "'." << Endl;
      }
   }

   fNVariables = nvars;
   fNTargets = ntgts;
   fNSpectators = nspcts;

   fVariableTypesAreCounted = true;
}

////////////////////////////////////////////////////////////////////////////////
/// TODO --> adapt to variable,target,spectator selection
/// method to calculate minimum, maximum, mean, and RMS for all
/// variables used in the MVA

void TMVA::VariableTransformBase::CalcNorm( const std::vector<const Event*>& events )
{
   if (!IsCreated()) return;

   const UInt_t nvars = GetNVariables();
   const UInt_t ntgts = GetNTargets();

   UInt_t nevts = events.size();

   TVectorD x2( nvars+ntgts ); x2 *= 0;
   TVectorD x0( nvars+ntgts ); x0 *= 0;
   TVectorD v0( nvars+ntgts ); v0 *= 0;

   Double_t sumOfWeights = 0;
   for (UInt_t ievt=0; ievt<nevts; ievt++) {
      const Event* ev = events[ievt];

      Double_t weight = ev->GetWeight();
      sumOfWeights += weight;
      for (UInt_t ivar=0; ivar<nvars; ivar++) {
         Double_t x = ev->GetValue(ivar);
         if (ievt==0) {
            Variables().at(ivar).SetMin(x);
            Variables().at(ivar).SetMax(x);
         }
         else {
            UpdateNorm( ivar,  x );
         }
         x0(ivar) += x*weight;
         x2(ivar) += x*x*weight;
      }
      for (UInt_t itgt=0; itgt<ntgts; itgt++) {
         Double_t x = ev->GetTarget(itgt);
         if (ievt==0) {
            Targets().at(itgt).SetMin(x);
            Targets().at(itgt).SetMax(x);
         }
         else {
            UpdateNorm( nvars+itgt,  x );
         }
         x0(nvars+itgt) += x*weight;
         x2(nvars+itgt) += x*x*weight;
      }
   }

   if (sumOfWeights <= 0) {
      Log() << kFATAL << " the sum of event weights calculated for your input is == 0"
            << " or exactly: " << sumOfWeights << " there is obviously some problem..."<< Endl;
   }

   // set Mean and RMS
   for (UInt_t ivar=0; ivar<nvars; ivar++) {
      Double_t mean = x0(ivar)/sumOfWeights;

      Variables().at(ivar).SetMean( mean );
      if (x2(ivar)/sumOfWeights - mean*mean < 0) {
         Log() << kFATAL << " the RMS of your input variable " << ivar
               << " evaluates to an imaginary number: sqrt("<< x2(ivar)/sumOfWeights - mean*mean
               <<") .. sometimes related to a problem with outliers and negative event weights"
               << Endl;
      }
      Variables().at(ivar).SetRMS( TMath::Sqrt( x2(ivar)/sumOfWeights - mean*mean) );
   }
   for (UInt_t itgt=0; itgt<ntgts; itgt++) {
      Double_t mean = x0(nvars+itgt)/sumOfWeights;
      Targets().at(itgt).SetMean( mean );
      if (x2(nvars+itgt)/sumOfWeights - mean*mean < 0) {
         Log() << kFATAL << " the RMS of your target variable " << itgt
               << " evaluates to an imaginary number: sqrt(" << x2(nvars+itgt)/sumOfWeights - mean*mean
               <<") .. sometimes related to a problem with outliers and negative event weights"
               << Endl;
      }
      Targets().at(itgt).SetRMS( TMath::Sqrt( x2(nvars+itgt)/sumOfWeights - mean*mean) );
   }
   // calculate variance
   for (UInt_t ievt=0; ievt<nevts; ievt++) {
      const Event* ev = events[ievt];
      Double_t weight = ev->GetWeight();
      for (UInt_t ivar=0; ivar<nvars; ivar++) {
         Double_t x = ev->GetValue(ivar);
         Double_t mean = Variables().at(ivar).GetMean();
         v0(ivar) += weight*(x-mean)*(x-mean);
      }
      for (UInt_t itgt=0; itgt<ntgts; itgt++) {
         Double_t x = ev->GetTarget(itgt);
         Double_t mean = Targets().at(itgt).GetMean();
         v0(nvars+itgt) += weight*(x-mean)*(x-mean);
      }

   }

   // set variance
   for (UInt_t ivar=0; ivar<nvars; ivar++) {
      Double_t variance = v0(ivar)/sumOfWeights;
      Variables().at(ivar).SetVariance( variance );
      Log() << kINFO << "Variable " << Variables().at(ivar).GetExpression() <<" variance = " << variance << Endl;
   }
   for (UInt_t itgt=0; itgt<ntgts; itgt++) {
      Double_t variance = v0(nvars+itgt)/sumOfWeights;
      Targets().at(itgt).SetVariance( variance );
      Log() << kINFO << "Target " << Targets().at(itgt).GetExpression() <<" variance = " << variance << Endl;
   }

   Log() << kVERBOSE << "Set minNorm/maxNorm for variables to: " << Endl;
   Log() << std::setprecision(3);
   for (UInt_t ivar=0; ivar<GetNVariables(); ivar++)
      Log() << "    " << Variables().at(ivar).GetInternalName()
            << "\t: [" << Variables().at(ivar).GetMin() << "\t, " << Variables().at(ivar).GetMax() << "\t] " << Endl;
   Log() << kVERBOSE << "Set minNorm/maxNorm for targets to: " << Endl;
   Log() << std::setprecision(3);
   for (UInt_t itgt=0; itgt<GetNTargets(); itgt++)
      Log() << "    " << Targets().at(itgt).GetInternalName()
            << "\t: [" << Targets().at(itgt).GetMin() << "\t, " << Targets().at(itgt).GetMax() << "\t] " << Endl;
   Log() << std::setprecision(5); // reset to better value
}

////////////////////////////////////////////////////////////////////////////////
/// TODO --> adapt to variable,target,spectator selection
/// default transformation output
/// --> only indicate that transformation occurred

std::vector<TString>* TMVA::VariableTransformBase::GetTransformationStrings( Int_t /*cls*/ ) const
{
   std::vector<TString>* strVec = new std::vector<TString>;
   for (UInt_t ivar=0; ivar<GetNVariables(); ivar++) {
      strVec->push_back( Variables()[ivar].GetLabel() + "_[transformed]");
   }

   return strVec;
}

////////////////////////////////////////////////////////////////////////////////
/// TODO --> adapt to variable,target,spectator selection
/// update min and max of a given variable (target) and a given transformation method

void TMVA::VariableTransformBase::UpdateNorm ( Int_t ivar,  Double_t x )
{
   Int_t nvars = fDsi.GetNVariables();
   if( ivar < nvars ){
      if (x < Variables().at(ivar).GetMin()) Variables().at(ivar).SetMin(x);
      if (x > Variables().at(ivar).GetMax()) Variables().at(ivar).SetMax(x);
   }else{
      if (x < Targets().at(ivar-nvars).GetMin()) Targets().at(ivar-nvars).SetMin(x);
      if (x > Targets().at(ivar-nvars).GetMax()) Targets().at(ivar-nvars).SetMax(x);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description the transformation (write out info of selected variables)

void TMVA::VariableTransformBase::AttachXMLTo(void* parent)
{
   void* selxml = gTools().AddChild(parent, "Selection");

   void* inpxml = gTools().AddChild(selxml, "Input");
   gTools().AddAttr(inpxml, "NInputs", fGet.size() );

   // choose the new dsi for output if present, if not, take the common one
   const DataSetInfo* outputDsiPtr = (fDsiOutput? fDsiOutput : &fDsi );

   for( ItVarTypeIdx itGet = fGet.begin(), itGetEnd = fGet.end(); itGet != itGetEnd; ++itGet ) {
      UInt_t idx  = (*itGet).second;
      Char_t type = (*itGet).first;

      TString label      = "";
      TString expression = "";
      TString typeString = "";
      switch( type ){
      case 'v':
         typeString = "Variable";
         label = fDsi.GetVariableInfo( idx ).GetLabel();
         expression = fDsi.GetVariableInfo( idx ).GetExpression();
         break;
      case 't':
         typeString = "Target";
         label = fDsi.GetTargetInfo( idx ).GetLabel();
         expression = fDsi.GetTargetInfo( idx ).GetExpression();
         break;
      case 's':
         typeString = "Spectator";
         label = fDsi.GetSpectatorInfo( idx ).GetLabel();
         expression = fDsi.GetSpectatorInfo( idx ).GetExpression();
         break;
      default:
         Log() << kFATAL << "VariableTransformBase/AttachXMLTo unknown variable type '" << type << "'." << Endl;
      }

      void* idxxml = gTools().AddChild(inpxml, "Input");
      //      gTools().AddAttr(idxxml, "Index", idx);
      gTools().AddAttr(idxxml, "Type",  typeString);
      gTools().AddAttr(idxxml, "Label", label);
      gTools().AddAttr(idxxml, "Expression", expression);
   }


   void* outxml = gTools().AddChild(selxml, "Output");
   gTools().AddAttr(outxml, "NOutputs", fPut.size() );

   for( ItVarTypeIdx itPut = fPut.begin(), itPutEnd = fPut.end(); itPut != itPutEnd; ++itPut ) {
      UInt_t idx  = (*itPut).second;
      Char_t type = (*itPut).first;

      TString label = "";
      TString expression = "";
      TString typeString = "";
      switch( type ){
      case 'v':
         typeString = "Variable";
         label = outputDsiPtr->GetVariableInfo( idx ).GetLabel();
         expression = outputDsiPtr->GetVariableInfo( idx ).GetExpression();
         break;
      case 't':
         typeString = "Target";
         label = outputDsiPtr->GetTargetInfo( idx ).GetLabel();
         expression = outputDsiPtr->GetTargetInfo( idx ).GetExpression();
         break;
      case 's':
         typeString = "Spectator";
         label = outputDsiPtr->GetSpectatorInfo( idx ).GetLabel();
         expression = outputDsiPtr->GetSpectatorInfo( idx ).GetExpression();
         break;
      default:
         Log() << kFATAL << "VariableTransformBase/AttachXMLTo unknown variable type '" << type << "'." << Endl;
      }

      void* idxxml = gTools().AddChild(outxml, "Output");
      //      gTools().AddAttr(idxxml, "Index", idx);
      gTools().AddAttr(idxxml, "Type",  typeString);
      gTools().AddAttr(idxxml, "Label", label);
      gTools().AddAttr(idxxml, "Expression", expression);
   }


}

////////////////////////////////////////////////////////////////////////////////
/// Read the input variables from the XML node

void TMVA::VariableTransformBase::ReadFromXML( void* selnode )
{
   void* inpnode = gTools().GetChild( selnode );
   void* outnode = gTools().GetNextChild( inpnode );

   UInt_t nvars  = GetNVariables();
   UInt_t ntgts  = GetNTargets();
   UInt_t nspcts = GetNSpectators();

   // read inputs
   fGet.clear();

   UInt_t nInputs = 0;
   gTools().ReadAttr(inpnode, "NInputs", nInputs);

   void* ch = gTools().GetChild( inpnode );
   while(ch) {
      TString typeString = "";
      TString label      = "";
      TString expression = "";

      gTools().ReadAttr(ch, "Type",  typeString);
      gTools().ReadAttr(ch, "Label", label);
      gTools().ReadAttr(ch, "Expression", expression);

      if( typeString == "Variable"  ){
         for( UInt_t ivar = 0; ivar < nvars; ++ivar ) { // search all variables
            if( fDsi.GetVariableInfo( ivar ).GetLabel() == label ||
                fDsi.GetVariableInfo( ivar ).GetExpression() == expression) {
               fGet.push_back( std::pair<Char_t,UInt_t>('v',ivar) );
               break;
            }
         }
      }else if( typeString == "Target"    ){
         for( UInt_t itgt = 0; itgt < ntgts; ++itgt ) { // search all targets
            if( fDsi.GetTargetInfo( itgt ).GetLabel() == label ||
                fDsi.GetTargetInfo( itgt ).GetExpression() == expression ) {
               fGet.push_back( std::pair<Char_t,UInt_t>('t',itgt) );
               break;
            }
         }
      }else if( typeString == "Spectator" ){
         for( UInt_t ispct = 0; ispct < nspcts; ++ispct ) { // search all spectators
            if( fDsi.GetSpectatorInfo( ispct ).GetLabel() == label ||
                fDsi.GetSpectatorInfo( ispct ).GetExpression() == expression ) {
               fGet.push_back( std::pair<Char_t,UInt_t>('s',ispct) );
               break;
            }
         }
      }else{
         Log() << kFATAL << "VariableTransformationBase/ReadFromXML : unknown type '" << typeString << "'." << Endl;
      }
      ch = gTools().GetNextChild( ch );
   }

   assert( nInputs == fGet.size() );

   // read outputs
   fPut.clear();

   UInt_t nOutputs = 0;
   gTools().ReadAttr(outnode, "NOutputs", nOutputs);

   void* chOut = gTools().GetChild( outnode );
   while(chOut) {
      TString typeString = "";
      TString label      = "";
      TString expression = "";

      gTools().ReadAttr(chOut, "Type",  typeString);
      gTools().ReadAttr(chOut, "Label", label);
      gTools().ReadAttr(chOut, "Expression", expression);

      if( typeString == "Variable"  ){
         for( UInt_t ivar = 0; ivar < nvars; ++ivar ) { // search all variables
            if( fDsi.GetVariableInfo( ivar ).GetLabel() == label ||
                fDsi.GetVariableInfo( ivar ).GetExpression() == expression ) {
               fPut.push_back( std::pair<Char_t,UInt_t>('v',ivar) );
               break;
            }
         }
      }else if( typeString == "Target"    ){
         for( UInt_t itgt = 0; itgt < ntgts; ++itgt ) { // search all targets
            if( fDsi.GetTargetInfo( itgt ).GetLabel() == label ||
                fDsi.GetTargetInfo( itgt ).GetExpression() == expression ) {
               fPut.push_back( std::pair<Char_t,UInt_t>('t',itgt) );
               break;
            }
         }
      }else if( typeString == "Spectator" ){
         for( UInt_t ispct = 0; ispct < nspcts; ++ispct ) { // search all spectators
            if( fDsi.GetSpectatorInfo( ispct ).GetLabel() == label ||
                fDsi.GetSpectatorInfo( ispct ).GetExpression() == expression ) {
               fPut.push_back( std::pair<Char_t,UInt_t>('s',ispct) );
               break;
            }
         }
      }else{
         Log() << kFATAL << "VariableTransformationBase/ReadFromXML : unknown type '" << typeString << "'." << Endl;
      }
      chOut = gTools().GetNextChild( chOut );
   }

   assert( nOutputs == fPut.size() );
}

////////////////////////////////////////////////////////////////////////////////
/// getinput and setoutput equivalent

void TMVA::VariableTransformBase::MakeFunction( std::ostream& fout, const TString& /*fncName*/, Int_t part,
                                                UInt_t /*trCounter*/, Int_t /*cls*/ )
{
   if( part == 0 ){ // definitions
      fout << std::endl;
      fout << "   // define the indices of the variables which are transformed by this transformation" << std::endl;
      fout << "   static std::vector<int> indicesGet;" << std::endl;
      fout << "   static std::vector<int> indicesPut;" << std::endl << std::endl;
      fout << "   if ( indicesGet.empty() ) {" << std::endl;
      fout << "      indicesGet.reserve(fNvars);" << std::endl;

      for( ItVarTypeIdxConst itEntry = fGet.begin(), itEntryEnd = fGet.end(); itEntry != itEntryEnd; ++itEntry ) {
         Char_t type = (*itEntry).first;
         Int_t  idx  = (*itEntry).second;

         switch( type ) {
         case 'v':
            fout << "      indicesGet.push_back( " << idx << ");" << std::endl;
            break;
         case 't':
            Log() << kWARNING << "MakeClass doesn't work with transformation of targets. The results will be wrong!" << Endl;
            break;
         case 's':
            Log() << kWARNING << "MakeClass doesn't work with transformation of spectators. The results will be wrong!" << Endl;
            break;
         default:
            Log() << kFATAL << "VariableTransformBase/GetInput : unknown type '" << type << "'." << Endl;
         }
      }
      fout << "   }" <<  std::endl;
      fout << "   if ( indicesPut.empty() ) {" << std::endl;
      fout << "      indicesPut.reserve(fNvars);" << std::endl;

      for( ItVarTypeIdxConst itEntry = fPut.begin(), itEntryEnd = fPut.end(); itEntry != itEntryEnd; ++itEntry ) {
         Char_t type = (*itEntry).first;
         Int_t  idx  = (*itEntry).second;

         switch( type ) {
         case 'v':
            fout << "      indicesPut.push_back( " << idx << ");" << std::endl;
            break;
         case 't':
            Log() << kWARNING << "MakeClass doesn't work with transformation of targets. The results will be wrong!" << Endl;
            break;
         case 's':
            Log() << kWARNING << "MakeClass doesn't work with transformation of spectators. The results will be wrong!" << Endl;
            break;
         default:
            Log() << kFATAL << "VariableTransformBase/PutInput : unknown type '" << type << "'." << Endl;
         }
      }

      fout << "   }" <<  std::endl;
      fout << std::endl;

   }else if( part == 1){
   }
}
