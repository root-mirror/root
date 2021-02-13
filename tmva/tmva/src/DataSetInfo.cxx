// @(#)root/tmva $Id$
// Author: Joerg Stelzer, Peter Speckmeier

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DataSetInfo                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <speckmay@mail.cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - DESY, Germany                  *
 *                                                                                *
 * Copyright (c) 2008:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      DESY Hamburg, Germany                                                     *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::DataSetInfo
\ingroup TMVA

Class that contains all the data information.

*/

#include <vector>

#include "TEventList.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TROOT.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/Event.h"

#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::DataSetInfo::DataSetInfo(const TString& name)
   : TObject(),
     fDataSetManager(NULL),
     fName(name),
     fDataSet( 0 ),
     fNeedsRebuilding( kTRUE ),
     fVariables(),
     fTargets(),
     fSpectators(),
     fClasses( 0 ),
     fNormalization( "NONE" ),
     fSplitOptions(""),
     fTrainingSumSignalWeights(-1),
     fTrainingSumBackgrWeights(-1),
     fTestingSumSignalWeights (-1),
     fTestingSumBackgrWeights (-1),
     fOwnRootDir(0),
     fVerbose( kFALSE ),
     fSignalClass(0),
     fTargetsForMulticlass(0),
     fLogger( new MsgLogger("DataSetInfo", kINFO) )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::DataSetInfo::~DataSetInfo()
{
   ClearDataSet();

   for(UInt_t i=0, iEnd = fClasses.size(); i<iEnd; ++i) {
      delete fClasses[i];
   }

   delete fTargetsForMulticlass;

   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::DataSetInfo::ClearDataSet() const
{
   if(fDataSet!=0) { delete fDataSet; fDataSet=0; }
}

////////////////////////////////////////////////////////////////////////////////

void
TMVA::DataSetInfo::SetMsgType( EMsgType t ) const
{
   fLogger->SetMinType(t);
}

////////////////////////////////////////////////////////////////////////////////

TMVA::ClassInfo* TMVA::DataSetInfo::AddClass( const TString& className )
{
   ClassInfo* theClass = GetClassInfo(className);
   if (theClass) return theClass;


   fClasses.push_back( new ClassInfo(className) );
   fClasses.back()->SetNumber(fClasses.size()-1);

   //Log() << kHEADER << Endl;

   Log() << kHEADER << Form("[%s] : ",fName.Data()) << "Added class \"" << className << "\""<< Endl;

   Log() << kDEBUG <<"\t with internal class number " << fClasses.back()->GetNumber() << Endl;


   if (className == "Signal") fSignalClass = fClasses.size()-1;  // store the signal class index ( for comparison reasons )

   return fClasses.back();
}

////////////////////////////////////////////////////////////////////////////////

TMVA::ClassInfo* TMVA::DataSetInfo::GetClassInfo( const TString& name ) const
{
   for (std::vector<ClassInfo*>::iterator it = fClasses.begin(); it < fClasses.end(); ++it) {
      if ((*it)->GetName() == name) return (*it);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::ClassInfo* TMVA::DataSetInfo::GetClassInfo( Int_t cls ) const
{
   try {
      return fClasses.at(cls);
   }
   catch(...) {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::DataSetInfo::PrintClasses() const
{
   for (UInt_t cls = 0; cls < GetNClasses() ; cls++) {
      Log() << kINFO << Form("Dataset[%s] : ",fName.Data()) << "Class index : " << cls << "  name : " << GetClassInfo(cls)->GetName() << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::DataSetInfo::IsSignal( const TMVA::Event* ev ) const
{
   return (ev->GetClass()  == fSignalClass);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Float_t>*  TMVA::DataSetInfo::GetTargetsForMulticlass( const TMVA::Event* ev )
{
   if( !fTargetsForMulticlass ) fTargetsForMulticlass = new std::vector<Float_t>( GetNClasses() );
   //   fTargetsForMulticlass->resize( GetNClasses() );
   fTargetsForMulticlass->assign( GetNClasses(), 0.0 );
   fTargetsForMulticlass->at( ev->GetClass() ) = 1.0;
   return fTargetsForMulticlass;
}


////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::DataSetInfo::HasCuts() const
{
   Bool_t hasCuts = kFALSE;
   for (std::vector<ClassInfo*>::iterator it = fClasses.begin(); it < fClasses.end(); ++it) {
      if( TString((*it)->GetCut()) != TString("") ) hasCuts = kTRUE;
   }
   return hasCuts;
}

////////////////////////////////////////////////////////////////////////////////

const TMatrixD* TMVA::DataSetInfo::CorrelationMatrix( const TString& className ) const
{
   ClassInfo* ptr = GetClassInfo(className);
   return ptr?ptr->GetCorrelationMatrix():0;
}

////////////////////////////////////////////////////////////////////////////////
/// add a variable (can be a complex expression) to the set of
/// variables used in the MV analysis

TMVA::VariableInfo& TMVA::DataSetInfo::AddVariable( const TString& expression,
                                                    const TString& title,
                                                    const TString& unit,
                                                    Double_t min, Double_t max,
                                                    char varType,
                                                    Bool_t normalized,
                                                    void* external )
{
   TString regexpr = expression; // remove possible blanks
   regexpr.ReplaceAll(" ", "" );
   fVariables.push_back(VariableInfo( regexpr, title, unit,
                                      fVariables.size()+1, varType, external, min, max, normalized ));
   fNeedsRebuilding = kTRUE;
   return fVariables.back();
}

////////////////////////////////////////////////////////////////////////////////
/// add variable with given VariableInfo

TMVA::VariableInfo& TMVA::DataSetInfo::AddVariable( const VariableInfo& varInfo){
   fVariables.push_back(VariableInfo( varInfo ));
   fNeedsRebuilding = kTRUE;
   return fVariables.back();
}

////////////////////////////////////////////////////////////////////////////////
/// add an  array of variables identified by an expression corresponding to an array entry in the tree

void TMVA::DataSetInfo::AddVariablesArray(const TString &expression, Int_t size, const TString &title, const TString &unit,
                                                   Double_t min, Double_t max, char varType, Bool_t normalized,
                                                   void *external)
{
   TString regexpr = expression; // remove possible blanks
   regexpr.ReplaceAll(" ", "");
   fVariables.reserve(fVariables.size() + size);
   for (int i = 0; i < size; ++i) {
      TString newTitle = title + TString::Format("[%d]", i);

      fVariables.emplace_back(regexpr, newTitle, unit, fVariables.size() + 1, varType, external, min, max, normalized);
      // set corresponding bit indicating is a variable from an array
      fVariables.back().SetBit(kIsArrayVariable);
      TString newVarName = fVariables.back().GetInternalName() + TString::Format("[%d]", i);
      fVariables.back().SetInternalName(newVarName);
   }
   fVarArrays[regexpr] = size;
   fNeedsRebuilding = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// add a variable (can be a complex expression) to the set of
/// variables used in the MV analysis

TMVA::VariableInfo& TMVA::DataSetInfo::AddTarget( const TString& expression,
                                                  const TString& title,
                                                  const TString& unit,
                                                  Double_t min, Double_t max,
                                                  Bool_t normalized,
                                                  void* external )
{
   TString regexpr = expression; // remove possible blanks
   regexpr.ReplaceAll(" ", "" );
   char type='F';
   fTargets.push_back(VariableInfo( regexpr, title, unit,
                                    fTargets.size()+1, type, external, min,
                                    max, normalized ));
   fNeedsRebuilding = kTRUE;
   return fTargets.back();
}

////////////////////////////////////////////////////////////////////////////////
/// add target with given VariableInfo

TMVA::VariableInfo& TMVA::DataSetInfo::AddTarget( const VariableInfo& varInfo){
   fTargets.push_back(VariableInfo( varInfo ));
   fNeedsRebuilding = kTRUE;
   return fTargets.back();
}

////////////////////////////////////////////////////////////////////////////////
/// add a spectator (can be a complex expression) to the set of spectator variables used in
/// the MV analysis

TMVA::VariableInfo& TMVA::DataSetInfo::AddSpectator( const TString& expression,
                                                     const TString& title,
                                                     const TString& unit,
                                                     Double_t min, Double_t max, char type,
                                                     Bool_t normalized, void* external )
{
   TString regexpr = expression; // remove possible blanks
   regexpr.ReplaceAll(" ", "" );
   fSpectators.push_back(VariableInfo( regexpr, title, unit,
                                       fSpectators.size()+1, type, external, min, max, normalized ));
   fNeedsRebuilding = kTRUE;
   return fSpectators.back();
}

////////////////////////////////////////////////////////////////////////////////
/// add spectator with given VariableInfo

TMVA::VariableInfo& TMVA::DataSetInfo::AddSpectator( const VariableInfo& varInfo){
   fSpectators.push_back(VariableInfo( varInfo ));
   fNeedsRebuilding = kTRUE;
   return fSpectators.back();
}

////////////////////////////////////////////////////////////////////////////////
/// find variable by name

Int_t TMVA::DataSetInfo::FindVarIndex(const TString& var) const
{
   for (UInt_t ivar=0; ivar<GetNVariables(); ivar++)
      if (var == GetVariableInfo(ivar).GetInternalName()) return ivar;

   for (UInt_t ivar=0; ivar<GetNVariables(); ivar++)
      Log() << kINFO  << Form("Dataset[%s] : ",fName.Data()) <<  GetVariableInfo(ivar).GetInternalName() << Endl;

   Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << "<FindVarIndex> Variable \'" << var << "\' not found." << Endl;

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// set the weight expressions for the classes
/// if class name is specified, set only for this class
/// if class name is unknown, register new class with this name

void TMVA::DataSetInfo::SetWeightExpression( const TString& expr, const TString& className )
{
   if (className != "") {
      TMVA::ClassInfo* ci = AddClass(className);
      ci->SetWeight( expr );
   }
   else {
      // no class name specified, set weight for all classes
      if (fClasses.empty()) {
         Log() << kWARNING << Form("Dataset[%s] : ",fName.Data()) << "No classes registered yet, cannot specify weight expression!" << Endl;
      }
      for (std::vector<ClassInfo*>::iterator it = fClasses.begin(); it < fClasses.end(); ++it) {
         (*it)->SetWeight( expr );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::DataSetInfo::SetCorrelationMatrix( const TString& className, TMatrixD* matrix )
{
   GetClassInfo(className)->SetCorrelationMatrix(matrix);
}

////////////////////////////////////////////////////////////////////////////////
/// set the cut for the classes

void TMVA::DataSetInfo::SetCut( const TCut& cut, const TString& className )
{
   if (className == "") {  // if no className has been given set the cut for all the classes
      for (std::vector<ClassInfo*>::iterator it = fClasses.begin(); it < fClasses.end(); ++it) {
         (*it)->SetCut( cut );
      }
   }
   else {
      TMVA::ClassInfo* ci = AddClass(className);
      ci->SetCut( cut );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// set the cut for the classes

void TMVA::DataSetInfo::AddCut( const TCut& cut, const TString& className )
{
   if (className == "") {  // if no className has been given set the cut for all the classes
      for (std::vector<ClassInfo*>::iterator it = fClasses.begin(); it < fClasses.end(); ++it) {
         const TCut& oldCut = (*it)->GetCut();
         (*it)->SetCut( oldCut+cut );
      }
   }
   else {
      TMVA::ClassInfo* ci = AddClass(className);
      ci->SetCut( ci->GetCut()+cut );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// returns list of variables

std::vector<TString> TMVA::DataSetInfo::GetListOfVariables() const
{
   std::vector<TString> vNames;
   std::vector<TMVA::VariableInfo>::const_iterator viIt = GetVariableInfos().begin();
   for(;viIt != GetVariableInfos().end(); ++viIt) vNames.push_back( (*viIt).GetInternalName() );

   return vNames;
}

////////////////////////////////////////////////////////////////////////////////
/// calculates the correlation matrices for signal and background,
/// prints them to standard output, and fills 2D histograms

void TMVA::DataSetInfo::PrintCorrelationMatrix( const TString& className )
{

   Log() << kHEADER //<< Form("Dataset[%s] : ",fName.Data())
    << "Correlation matrix (" << className << "):" << Endl;
   gTools().FormattedOutput( *CorrelationMatrix( className ), GetListOfVariables(), Log() );
}

////////////////////////////////////////////////////////////////////////////////

TH2* TMVA::DataSetInfo::CreateCorrelationMatrixHist( const TMatrixD* m,
                                                     const TString&  hName,
                                                     const TString&  hTitle ) const
{
   if (m==0) return 0;

   const UInt_t nvar = GetNVariables();

   // workaround till the TMatrix templates are commonly used
   // this keeps backward compatibility
   TMatrixF* tm = new TMatrixF( nvar, nvar );
   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      for (UInt_t jvar=0; jvar<nvar; jvar++) {
         (*tm)(ivar, jvar) = (*m)(ivar,jvar);
      }
   }

   TH2F* h2 = new TH2F( *tm );
   h2->SetNameTitle( hName, hTitle );

   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      h2->GetXaxis()->SetBinLabel( ivar+1, GetVariableInfo(ivar).GetTitle() );
      h2->GetYaxis()->SetBinLabel( ivar+1, GetVariableInfo(ivar).GetTitle() );
   }

   // present in percent, and round off digits
   // also, use absolute value of correlation coefficient (ignore sign)
   h2->Scale( 100.0  );
   for (UInt_t ibin=1; ibin<=nvar; ibin++) {
      for (UInt_t jbin=1; jbin<=nvar; jbin++) {
         h2->SetBinContent( ibin, jbin, Int_t(h2->GetBinContent( ibin, jbin )) );
      }
   }

   // style settings
   const Float_t labelSize = 0.055;
   h2->SetStats( 0 );
   h2->GetXaxis()->SetLabelSize( labelSize );
   h2->GetYaxis()->SetLabelSize( labelSize );
   h2->SetMarkerSize( 1.5 );
   h2->SetMarkerColor( 0 );
   h2->LabelsOption( "d" ); // diagonal labels on x axis
   h2->SetLabelOffset( 0.011 );// label offset on x axis
   h2->SetMinimum( -100.0 );
   h2->SetMaximum( +100.0 );

   // -------------------------------------------------------------------------------------
   // just in case one wants to change the position of the color palette axis
   // -------------------------------------------------------------------------------------
   //     gROOT->SetStyle("Plain");
   //     TStyle* gStyle = gROOT->GetStyle( "Plain" );
   //     gStyle->SetPalette( 1, 0 );
   //     TPaletteAxis* paletteAxis
   //                   = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject( "palette" );
   // -------------------------------------------------------------------------------------

   Log() << kDEBUG << Form("Dataset[%s] : ",fName.Data()) << "Created correlation matrix as 2D histogram: " << h2->GetName() << Endl;

   return h2;
}

////////////////////////////////////////////////////////////////////////////////
/// returns data set

TMVA::DataSet* TMVA::DataSetInfo::GetDataSet() const
{
   if (fDataSet==0 || fNeedsRebuilding) {
      if(fDataSet!=0) ClearDataSet();
      //      fDataSet = DataSetManager::Instance().CreateDataSet(GetName()); //DSMTEST replaced by following lines
      if( !fDataSetManager )
         Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << "DataSetManager has not been set in DataSetInfo (GetDataSet() )." << Endl;
      fDataSet = fDataSetManager->CreateDataSet(GetName());

      fNeedsRebuilding = kFALSE;
   }
   return fDataSet;
}

////////////////////////////////////////////////////////////////////////////////

UInt_t TMVA::DataSetInfo::GetNSpectators(bool all) const
{
   if(all)
      return fSpectators.size();
   UInt_t nsp(0);
   for(std::vector<VariableInfo>::const_iterator spit=fSpectators.begin(); spit!=fSpectators.end(); ++spit) {
      if(spit->GetVarType()!='C') nsp++;
   }
   return nsp;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TMVA::DataSetInfo::GetClassNameMaxLength() const
{
   Int_t maxL = 0;
   for (UInt_t cl = 0; cl < GetNClasses(); cl++) {
      if (TString(GetClassInfo(cl)->GetName()).Length() > maxL) maxL = TString(GetClassInfo(cl)->GetName()).Length();
   }

   return maxL;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TMVA::DataSetInfo::GetVariableNameMaxLength() const
{
   Int_t maxL = 0;
   for (UInt_t i = 0; i < GetNVariables(); i++) {
      if (TString(GetVariableInfo(i).GetExpression()).Length() > maxL) maxL = TString(GetVariableInfo(i).GetExpression()).Length();
   }

   return maxL;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TMVA::DataSetInfo::GetTargetNameMaxLength() const
{
   Int_t maxL = 0;
   for (UInt_t i = 0; i < GetNTargets(); i++) {
      if (TString(GetTargetInfo(i).GetExpression()).Length() > maxL) maxL = TString(GetTargetInfo(i).GetExpression()).Length();
   }

   return maxL;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::DataSetInfo::GetTrainingSumSignalWeights(){
   if (fTrainingSumSignalWeights<0) Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << " asking for the sum of training signal event weights which is not initialized yet" << Endl;
   return fTrainingSumSignalWeights;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::DataSetInfo::GetTrainingSumBackgrWeights(){
   if (fTrainingSumBackgrWeights<0) Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << " asking for the sum of training backgr event weights which is not initialized yet" << Endl;
   return fTrainingSumBackgrWeights;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::DataSetInfo::GetTestingSumSignalWeights (){
   if (fTestingSumSignalWeights<0) Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << " asking for the sum of testing signal event weights which is not initialized yet" << Endl;
   return fTestingSumSignalWeights ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::DataSetInfo::GetTestingSumBackgrWeights (){
   if (fTestingSumBackgrWeights<0) Log() << kFATAL << Form("Dataset[%s] : ",fName.Data()) << " asking for the sum of testing backgr event weights which is not initialized yet" << Endl;
   return fTestingSumBackgrWeights ;
}
