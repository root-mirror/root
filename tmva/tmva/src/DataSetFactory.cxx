// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Eckhard von Toerne, Helge Voss

/*****************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis  *
 * Package: TMVA                                                             *
 * Class  : DataSetFactory                                                   *
 * Web    : http://tmva.sourceforge.net                                      *
 *                                                                           *
 * Description:                                                              *
 *      Implementation (see header for description)                          *
 *                                                                           *
 * Authors (alphabetical):                                                   *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland         *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland      *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - MSU, USA                  *
 *      Eckhard von Toerne <evt@physik.uni-bonn.de>  - U. of Bonn, Germany   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany *
 *                                                                           *
 * Copyright (c) 2009:                                                       *
 *      CERN, Switzerland                                                    *
 *      MPI-K Heidelberg, Germany                                            *
 *      U. of Bonn, Germany                                                  *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted according to the terms listed in LICENSE      *
 * (http://tmva.sourceforge.net/LICENSE)                                     *
 *****************************************************************************/

/*! \class TMVA::DataSetFactory
\ingroup TMVA

Class that contains all the data information

*/

#include <assert.h>

#include <map>
#include <vector>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>

#include "TMVA/DataSetFactory.h"

#include "TEventList.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/Configurable.h"
#include "TMVA/VariableIdentityTransform.h"
#include "TMVA/VariableDecorrTransform.h"
#include "TMVA/VariablePCATransform.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/Event.h"

#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

using namespace std;

//TMVA::DataSetFactory* TMVA::DataSetFactory::fgInstance = 0;

namespace TMVA {
   // calculate the largest common divider
   // this function is not happy if numbers are negative!
   Int_t LargestCommonDivider(Int_t a, Int_t b)
   {
      if (a<b) {Int_t tmp = a; a=b; b=tmp; } // achieve a>=b
      if (b==0) return a;
      Int_t fullFits = a/b;
      return LargestCommonDivider(b,a-b*fullFits);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::DataSetFactory::DataSetFactory() :
   fVerbose(kFALSE),
   fVerboseLevel(TString("Info")),
   fScaleWithPreselEff(0),
   fCurrentTree(0),
   fCurrentEvtIdx(0),
   fInputFormulas(0),
   fLogger( new MsgLogger("DataSetFactory", kINFO) )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::DataSetFactory::~DataSetFactory()
{
   std::vector<TTreeFormula*>::const_iterator formIt;

   for (formIt = fInputFormulas.begin()    ; formIt!=fInputFormulas.end()    ; ++formIt) if (*formIt) delete *formIt;
   for (formIt = fTargetFormulas.begin()   ; formIt!=fTargetFormulas.end()   ; ++formIt) if (*formIt) delete *formIt;
   for (formIt = fCutFormulas.begin()      ; formIt!=fCutFormulas.end()      ; ++formIt) if (*formIt) delete *formIt;
   for (formIt = fWeightFormula.begin()    ; formIt!=fWeightFormula.end()    ; ++formIt) if (*formIt) delete *formIt;
   for (formIt = fSpectatorFormulas.begin(); formIt!=fSpectatorFormulas.end(); ++formIt) if (*formIt) delete *formIt;

   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// steering the creation of a new dataset

TMVA::DataSet* TMVA::DataSetFactory::CreateDataSet( TMVA::DataSetInfo& dsi,
                                                    TMVA::DataInputHandler& dataInput )
{
   // build the first dataset from the data input
   DataSet * ds = BuildInitialDataSet( dsi, dataInput );

   if (ds->GetNEvents() > 1 && fComputeCorrelations ) {
      CalcMinMax(ds,dsi);

      // from the the final dataset build the correlation matrix
      for (UInt_t cl = 0; cl< dsi.GetNClasses(); cl++) {
         const TString className = dsi.GetClassInfo(cl)->GetName();
         dsi.SetCorrelationMatrix( className, CalcCorrelationMatrix( ds, cl ) );
         if (fCorrelations) {
            dsi.PrintCorrelationMatrix(className);
         }
      }
      //Log() << kHEADER <<  Endl;
      Log() << kHEADER << Form("[%s] : ",dsi.GetName()) << " " << Endl << Endl;
   }

   return ds;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::DataSet* TMVA::DataSetFactory::BuildDynamicDataSet( TMVA::DataSetInfo& dsi )
{
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "Build DataSet consisting of one Event with dynamically changing variables" << Endl;
   DataSet* ds = new DataSet(dsi);

   // create a DataSet with one Event which uses dynamic variables
   // (pointers to variables)
   if(dsi.GetNClasses()==0){
      dsi.AddClass( "data" );
      dsi.GetClassInfo( "data" )->SetNumber(0);
   }

   std::vector<Float_t*>* evdyn = new std::vector<Float_t*>(0);

   std::vector<VariableInfo>& varinfos = dsi.GetVariableInfos();

   if (varinfos.empty())
      Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName()) << "Dynamic data set cannot be built, since no variable informations are present. Apparently no variables have been set. This should not happen, please contact the TMVA authors." << Endl;

   std::vector<VariableInfo>::iterator it = varinfos.begin(), itEnd=varinfos.end();
   for (;it!=itEnd;++it) {
      Float_t* external=(Float_t*)(*it).GetExternalLink();
      if (external==0)
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "The link to the external variable is NULL while I am trying to build a dynamic data set. In this case fTmpEvent from MethodBase HAS TO BE USED in the method to get useful values in variables." << Endl;
      else evdyn->push_back (external);
   }

   std::vector<VariableInfo>& spectatorinfos = dsi.GetSpectatorInfos();
   it = spectatorinfos.begin();
   for (;it!=spectatorinfos.end();++it) evdyn->push_back( (Float_t*)(*it).GetExternalLink() );

   TMVA::Event * ev = new Event((const std::vector<Float_t*>*&)evdyn, varinfos.size());
   std::vector<Event*>* newEventVector = new std::vector<Event*>;
   newEventVector->push_back(ev);

   ds->SetEventCollection(newEventVector, Types::kTraining);
   ds->SetCurrentType( Types::kTraining );
   ds->SetCurrentEvent( 0 );

   delete newEventVector;
   return ds;
}

////////////////////////////////////////////////////////////////////////////////
/// if no entries, than create a DataSet with one Event which uses
/// dynamic variables (pointers to variables)

TMVA::DataSet*
TMVA::DataSetFactory::BuildInitialDataSet( DataSetInfo& dsi,
                                           DataInputHandler& dataInput )
{
   if (dataInput.GetEntries()==0) return BuildDynamicDataSet( dsi );
   // -------------------------------------------------------------------------

   // register the classes in the datasetinfo-object
   // information comes from the trees in the dataInputHandler-object
   std::vector< TString >* classList = dataInput.GetClassList();
   for (std::vector<TString>::iterator it = classList->begin(); it< classList->end(); ++it) {
      dsi.AddClass( (*it) );
   }
   delete classList;

   EvtStatsPerClass eventCounts(dsi.GetNClasses());
   TString normMode;
   TString splitMode;
   TString mixMode;
   UInt_t  splitSeed;

   InitOptions( dsi, eventCounts, normMode, splitSeed, splitMode , mixMode );
   // ======= build event-vector from input, apply preselection ===============
   EventVectorOfClassesOfTreeType tmpEventVector;
   BuildEventVector( dsi, dataInput, tmpEventVector, eventCounts );

   DataSet* ds = MixEvents( dsi, tmpEventVector, eventCounts,
                            splitMode, mixMode, normMode, splitSeed );

   const Bool_t showCollectedOutput = kFALSE;
   if (showCollectedOutput) {
      Int_t maxL = dsi.GetClassNameMaxLength();
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Collected:" << Endl;
      for (UInt_t cl = 0; cl < dsi.GetNClasses(); cl++) {
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "    "
               << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
               << " training entries: " << ds->GetNClassEvents( 0, cl ) << Endl;
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "    "
               << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
               << " testing  entries: " << ds->GetNClassEvents( 1, cl ) << Endl;
      }
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " " << Endl;
   }

   return ds;
}

////////////////////////////////////////////////////////////////////////////////
/// checks a TTreeFormula for problems

Bool_t TMVA::DataSetFactory::CheckTTreeFormula( TTreeFormula* ttf,
                                                const TString& expression,
                                                Bool_t& hasDollar )
{
   Bool_t worked = kTRUE;

   if( ttf->GetNdim() <= 0 )
      Log() << kFATAL << "Expression " << expression.Data()
            << " could not be resolved to a valid formula. " << Endl;
   if( ttf->GetNdata() == 0 ){
      Log() << kWARNING << "Expression: " << expression.Data()
            << " does not provide data for this event. "
            << "This event is not taken into account. --> please check if you use as a variable "
            << "an entry of an array which is not filled for some events "
            << "(e.g. arr[4] when arr has only 3 elements)." << Endl;
      Log() << kWARNING << "If you want to take the event into account you can do something like: "
            << "\"Alt$(arr[4],0)\" where in cases where arr doesn't have a 4th element, "
            << " 0 is taken as an alternative." << Endl;
      worked = kFALSE;
   }
   if( expression.Contains("$") )
      hasDollar = kTRUE;
   else
      {
         for (int i = 0, iEnd = ttf->GetNcodes (); i < iEnd; ++i)
            {
               TLeaf* leaf = ttf->GetLeaf (i);
               if (!leaf->IsOnTerminalBranch())
                  hasDollar = kTRUE;
            }
      }
   return worked;
}


////////////////////////////////////////////////////////////////////////////////
/// While the data gets copied into the local training and testing
/// trees, the input tree can change (for instance when changing from
/// signal to background tree, or using TChains as input) The
/// TTreeFormulas, that hold the input expressions need to be
/// re-associated with the new tree, which is done here

void TMVA::DataSetFactory::ChangeToNewTree( TreeInfo& tinfo, const DataSetInfo & dsi )
{
   TTree *tr = tinfo.GetTree()->GetTree();

   //tr->SetBranchStatus("*",1); // nor needed when using TTReeFormula
   tr->ResetBranchAddresses();

   Bool_t hasDollar = kTRUE;  // Set to false if wants to enable only some branch in the tree

   // 1) the input variable formulas
   Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " create input formulas for tree " << tr->GetName() << Endl;
   std::vector<TTreeFormula*>::const_iterator formIt, formItEnd;
   for (formIt = fInputFormulas.begin(), formItEnd=fInputFormulas.end(); formIt!=formItEnd; ++formIt) if (*formIt) delete *formIt;
   fInputFormulas.clear();
   TTreeFormula* ttf = 0;
   fInputTableFormulas.clear();  // this contains shallow pointer copies

   bool firstArrayVar = kTRUE;
   int firstArrayVarIndex = -1;
   int arraySize = -1;
   for (UInt_t i = 0; i < dsi.GetNVariables(); i++) {

      // create TTreeformula
      if (! dsi.IsVariableFromArray(i) )  {
            ttf = new TTreeFormula(Form("Formula%s", dsi.GetVariableInfo(i).GetInternalName().Data()),
                                   dsi.GetVariableInfo(i).GetExpression().Data(), tr);
            CheckTTreeFormula(ttf, dsi.GetVariableInfo(i).GetExpression(), hasDollar);
            fInputFormulas.emplace_back(ttf);
            fInputTableFormulas.emplace_back(std::make_pair(ttf, (Int_t) 0));
      } else {
         // it is a variable from an array
         if (firstArrayVar) {

            // create a new TFormula
            ttf = new TTreeFormula(Form("Formula%s", dsi.GetVariableInfo(i).GetInternalName().Data()),
                                   dsi.GetVariableInfo(i).GetExpression().Data(), tr);
            CheckTTreeFormula(ttf, dsi.GetVariableInfo(i).GetExpression(), hasDollar);
            fInputFormulas.push_back(ttf);

            arraySize = dsi.GetVarArraySize(dsi.GetVariableInfo(i).GetExpression());
            firstArrayVar = kFALSE;
            firstArrayVarIndex = i;

            Log() << kINFO << "Using variable " << dsi.GetVariableInfo(i).GetInternalName() <<
               " from array expression " << dsi.GetVariableInfo(i).GetExpression() << " of size " << arraySize << Endl;
         }
         fInputTableFormulas.push_back(std::make_pair(ttf, (Int_t) i-firstArrayVarIndex));
         if (int(i)-firstArrayVarIndex == arraySize-1 ) {
            // I am the last element of the array
            firstArrayVar = kTRUE;
            firstArrayVarIndex = -1;
            Log() << kDEBUG << "Using Last variable from array : " << dsi.GetVariableInfo(i).GetInternalName() << Endl;
         }
      }

   }

   //
   // targets
   //
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "transform regression targets" << Endl;
   for (formIt = fTargetFormulas.begin(), formItEnd = fTargetFormulas.end(); formIt!=formItEnd; ++formIt) if (*formIt) delete *formIt;
   fTargetFormulas.clear();
   for (UInt_t i=0; i<dsi.GetNTargets(); i++) {
      ttf = new TTreeFormula( Form( "Formula%s", dsi.GetTargetInfo(i).GetInternalName().Data() ),
                              dsi.GetTargetInfo(i).GetExpression().Data(), tr );
      CheckTTreeFormula( ttf, dsi.GetTargetInfo(i).GetExpression(), hasDollar );
      fTargetFormulas.push_back( ttf );
   }

   //
   // spectators
   //
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "transform spectator variables" << Endl;
   for (formIt = fSpectatorFormulas.begin(), formItEnd = fSpectatorFormulas.end(); formIt!=formItEnd; ++formIt) if (*formIt) delete *formIt;
   fSpectatorFormulas.clear();
   for (UInt_t i=0; i<dsi.GetNSpectators(); i++) {
      ttf = new TTreeFormula( Form( "Formula%s", dsi.GetSpectatorInfo(i).GetInternalName().Data() ),
                              dsi.GetSpectatorInfo(i).GetExpression().Data(), tr );
      CheckTTreeFormula( ttf, dsi.GetSpectatorInfo(i).GetExpression(), hasDollar );
      fSpectatorFormulas.push_back( ttf );
   }

   //
   // the cuts (one per class, if non-existent: formula pointer = 0)
   //
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "transform cuts" << Endl;
   for (formIt = fCutFormulas.begin(), formItEnd = fCutFormulas.end(); formIt!=formItEnd; ++formIt) if (*formIt) delete *formIt;
   fCutFormulas.clear();
   for (UInt_t clIdx=0; clIdx<dsi.GetNClasses(); clIdx++) {
      const TCut& tmpCut = dsi.GetClassInfo(clIdx)->GetCut();
      const TString tmpCutExp(tmpCut.GetTitle());
      ttf = 0;
      if (tmpCutExp!="") {
         ttf = new TTreeFormula( Form("CutClass%i",clIdx), tmpCutExp, tr );
         Bool_t worked = CheckTTreeFormula( ttf, tmpCutExp, hasDollar );
         if( !worked ){
            Log() << kWARNING << "Please check class \"" << dsi.GetClassInfo(clIdx)->GetName()
                  << "\" cut \"" << dsi.GetClassInfo(clIdx)->GetCut() << Endl;
         }
      }
      fCutFormulas.push_back( ttf );
   }

   //
   // the weights (one per class, if non-existent: formula pointer = 0)
   //
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "transform weights" << Endl;
   for (formIt = fWeightFormula.begin(), formItEnd = fWeightFormula.end(); formIt!=formItEnd; ++formIt) if (*formIt) delete *formIt;
   fWeightFormula.clear();
   for (UInt_t clIdx=0; clIdx<dsi.GetNClasses(); clIdx++) {
      const TString tmpWeight = dsi.GetClassInfo(clIdx)->GetWeight();

      if (dsi.GetClassInfo(clIdx)->GetName() != tinfo.GetClassName() ) { // if the tree is of another class
         fWeightFormula.push_back( 0 );
         continue;
      }

      ttf = 0;
      if (tmpWeight!="") {
         ttf = new TTreeFormula( "FormulaWeight", tmpWeight, tr );
         Bool_t worked = CheckTTreeFormula( ttf, tmpWeight, hasDollar );
         if( !worked ){
            Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName()) << "Please check class \"" << dsi.GetClassInfo(clIdx)->GetName()
                  << "\" weight \"" << dsi.GetClassInfo(clIdx)->GetWeight() << Endl;
         }
      }
      else {
         ttf = 0;
      }
      fWeightFormula.push_back( ttf );
   }
   return;
   // all this code below is not needed when using TTReeFormula

   Log() << kDEBUG << Form("Dataset[%s] : ", dsi.GetName()) << "enable branches" << Endl;
   // now enable only branches that are needed in any input formula, target, cut, weight

   if (!hasDollar) {
      tr->SetBranchStatus("*",0);
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "enable branches: input variables" << Endl;
      // input vars
      for (formIt = fInputFormulas.begin(); formIt!=fInputFormulas.end(); ++formIt) {
         ttf = *formIt;
         for (Int_t bi = 0; bi<ttf->GetNcodes(); bi++) {
            tr->SetBranchStatus( ttf->GetLeaf(bi)->GetBranch()->GetName(), 1 );
         }
      }
      // targets
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "enable branches: targets" << Endl;
      for (formIt = fTargetFormulas.begin(); formIt!=fTargetFormulas.end(); ++formIt) {
         ttf = *formIt;
         for (Int_t bi = 0; bi<ttf->GetNcodes(); bi++)
            tr->SetBranchStatus( ttf->GetLeaf(bi)->GetBranch()->GetName(), 1 );
      }
      // spectators
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "enable branches: spectators" << Endl;
      for (formIt = fSpectatorFormulas.begin(); formIt!=fSpectatorFormulas.end(); ++formIt) {
         ttf = *formIt;
         for (Int_t bi = 0; bi<ttf->GetNcodes(); bi++)
            tr->SetBranchStatus( ttf->GetLeaf(bi)->GetBranch()->GetName(), 1 );
      }
      // cuts
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "enable branches: cuts" << Endl;
      for (formIt = fCutFormulas.begin(); formIt!=fCutFormulas.end(); ++formIt) {
         ttf = *formIt;
         if (!ttf) continue;
         for (Int_t bi = 0; bi<ttf->GetNcodes(); bi++)
            tr->SetBranchStatus( ttf->GetLeaf(bi)->GetBranch()->GetName(), 1 );
      }
      // weights
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName()) << "enable branches: weights" << Endl;
      for (formIt = fWeightFormula.begin(); formIt!=fWeightFormula.end(); ++formIt) {
         ttf = *formIt;
         if (!ttf) continue;
         for (Int_t bi = 0; bi<ttf->GetNcodes(); bi++)
            tr->SetBranchStatus( ttf->GetLeaf(bi)->GetBranch()->GetName(), 1 );
      }
   }
   Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "tree initialized" << Endl;
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// compute covariance matrix

void TMVA::DataSetFactory::CalcMinMax( DataSet* ds, TMVA::DataSetInfo& dsi )
{
   const UInt_t nvar  = ds->GetNVariables();
   const UInt_t ntgts = ds->GetNTargets();
   const UInt_t nvis  = ds->GetNSpectators();

   Float_t *min = new Float_t[nvar];
   Float_t *max = new Float_t[nvar];
   Float_t *tgmin = new Float_t[ntgts];
   Float_t *tgmax = new Float_t[ntgts];
   Float_t *vmin  = new Float_t[nvis];
   Float_t *vmax  = new Float_t[nvis];

   for (UInt_t ivar=0; ivar<nvar ; ivar++) {   min[ivar] = FLT_MAX;   max[ivar] = -FLT_MAX; }
   for (UInt_t ivar=0; ivar<ntgts; ivar++) { tgmin[ivar] = FLT_MAX; tgmax[ivar] = -FLT_MAX; }
   for (UInt_t ivar=0; ivar<nvis;  ivar++) {  vmin[ivar] = FLT_MAX;  vmax[ivar] = -FLT_MAX; }

   // perform event loop

   for (Int_t i=0; i<ds->GetNEvents(); i++) {
      const Event * ev = ds->GetEvent(i);
      for (UInt_t ivar=0; ivar<nvar; ivar++) {
         Double_t v = ev->GetValue(ivar);
         if (v<min[ivar]) min[ivar] = v;
         if (v>max[ivar]) max[ivar] = v;
      }
      for (UInt_t itgt=0; itgt<ntgts; itgt++) {
         Double_t v = ev->GetTarget(itgt);
         if (v<tgmin[itgt]) tgmin[itgt] = v;
         if (v>tgmax[itgt]) tgmax[itgt] = v;
      }
      for (UInt_t ivis=0; ivis<nvis; ivis++) {
         Double_t v = ev->GetSpectator(ivis);
         if (v<vmin[ivis]) vmin[ivis] = v;
         if (v>vmax[ivis]) vmax[ivis] = v;
      }
   }

   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      dsi.GetVariableInfo(ivar).SetMin(min[ivar]);
      dsi.GetVariableInfo(ivar).SetMax(max[ivar]);
      if( TMath::Abs(max[ivar]-min[ivar]) <= FLT_MIN )
         Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName()) << "Variable " << dsi.GetVariableInfo(ivar).GetExpression().Data() << " is constant. Please remove the variable." << Endl;
   }
   for (UInt_t ivar=0; ivar<ntgts; ivar++) {
      dsi.GetTargetInfo(ivar).SetMin(tgmin[ivar]);
      dsi.GetTargetInfo(ivar).SetMax(tgmax[ivar]);
      if( TMath::Abs(tgmax[ivar]-tgmin[ivar]) <= FLT_MIN )
         Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName()) << "Target " << dsi.GetTargetInfo(ivar).GetExpression().Data() << " is constant. Please remove the variable." << Endl;
   }
   for (UInt_t ivar=0; ivar<nvis; ivar++) {
      dsi.GetSpectatorInfo(ivar).SetMin(vmin[ivar]);
      dsi.GetSpectatorInfo(ivar).SetMax(vmax[ivar]);
      //       if( TMath::Abs(vmax[ivar]-vmin[ivar]) <= FLT_MIN )
      //          Log() << kWARNING << "Spectator variable " << dsi.GetSpectatorInfo(ivar).GetExpression().Data() << " is constant." << Endl;
   }
   delete [] min;
   delete [] max;
   delete [] tgmin;
   delete [] tgmax;
   delete [] vmin;
   delete [] vmax;
}

////////////////////////////////////////////////////////////////////////////////
/// computes correlation matrix for variables "theVars" in tree;
/// "theType" defines the required event "type"
/// ("type" variable must be present in tree)

TMatrixD* TMVA::DataSetFactory::CalcCorrelationMatrix( DataSet* ds, const UInt_t classNumber )
{
   // first compute variance-covariance
   TMatrixD* mat = CalcCovarianceMatrix( ds, classNumber );

   // now the correlation
   UInt_t nvar = ds->GetNVariables(), ivar, jvar;

   for (ivar=0; ivar<nvar; ivar++) {
      for (jvar=0; jvar<nvar; jvar++) {
         if (ivar != jvar) {
            Double_t d = (*mat)(ivar, ivar)*(*mat)(jvar, jvar);
            if (d > 0) (*mat)(ivar, jvar) /= sqrt(d);
            else {
               Log() << kWARNING << Form("Dataset[%s] : ",DataSetInfo().GetName())<< "<GetCorrelationMatrix> Zero variances for variables "
                     << "(" << ivar << ", " << jvar << ") = " << d
                     << Endl;
               (*mat)(ivar, jvar) = 0;
            }
         }
      }
   }

   for (ivar=0; ivar<nvar; ivar++) (*mat)(ivar, ivar) = 1.0;

   return mat;
}

////////////////////////////////////////////////////////////////////////////////
/// compute covariance matrix

TMatrixD* TMVA::DataSetFactory::CalcCovarianceMatrix( DataSet * ds, const UInt_t classNumber )
{
   UInt_t nvar = ds->GetNVariables();
   UInt_t ivar = 0, jvar = 0;

   TMatrixD* mat = new TMatrixD( nvar, nvar );

   // init matrices
   TVectorD vec(nvar);
   TMatrixD mat2(nvar, nvar);
   for (ivar=0; ivar<nvar; ivar++) {
      vec(ivar) = 0;
      for (jvar=0; jvar<nvar; jvar++) mat2(ivar, jvar) = 0;
   }

   // perform event loop
   Double_t ic = 0;
   for (Int_t i=0; i<ds->GetNEvents(); i++) {

      const Event * ev = ds->GetEvent(i);
      if (ev->GetClass() != classNumber ) continue;

      Double_t weight = ev->GetWeight();
      ic += weight; // count used events

      for (ivar=0; ivar<nvar; ivar++) {

         Double_t xi = ev->GetValue(ivar);
         vec(ivar) += xi*weight;
         mat2(ivar, ivar) += (xi*xi*weight);

         for (jvar=ivar+1; jvar<nvar; jvar++) {
            Double_t xj =  ev->GetValue(jvar);
            mat2(ivar, jvar) += (xi*xj*weight);
         }
      }
   }

   for (ivar=0; ivar<nvar; ivar++)
      for (jvar=ivar+1; jvar<nvar; jvar++)
         mat2(jvar, ivar) = mat2(ivar, jvar); // symmetric matrix


   // variance-covariance
   for (ivar=0; ivar<nvar; ivar++) {
      for (jvar=0; jvar<nvar; jvar++) {
         (*mat)(ivar, jvar) = mat2(ivar, jvar)/ic - vec(ivar)*vec(jvar)/(ic*ic);
      }
   }

   return mat;
}

// --------------------------------------- new versions

////////////////////////////////////////////////////////////////////////////////
/// the dataset splitting

void
TMVA::DataSetFactory::InitOptions( TMVA::DataSetInfo& dsi,
                                   EvtStatsPerClass& nEventRequests,
                                   TString& normMode,
                                   UInt_t&  splitSeed,
                                   TString& splitMode,
                                   TString& mixMode)
{
   Configurable splitSpecs( dsi.GetSplitOptions() );
   splitSpecs.SetConfigName("DataSetFactory");
   splitSpecs.SetConfigDescription( "Configuration options given in the \"PrepareForTrainingAndTesting\" call; these options define the creation of the data sets used for training and expert validation by TMVA" );

   splitMode = "Random";    // the splitting mode
   splitSpecs.DeclareOptionRef( splitMode, "SplitMode",
                                "Method of picking training and testing events (default: random)" );
   splitSpecs.AddPreDefVal(TString("Random"));
   splitSpecs.AddPreDefVal(TString("Alternate"));
   splitSpecs.AddPreDefVal(TString("Block"));

   mixMode = "SameAsSplitMode";    // the splitting mode
   splitSpecs.DeclareOptionRef( mixMode, "MixMode",
                                "Method of mixing events of different classes into one dataset (default: SameAsSplitMode)" );
   splitSpecs.AddPreDefVal(TString("SameAsSplitMode"));
   splitSpecs.AddPreDefVal(TString("Random"));
   splitSpecs.AddPreDefVal(TString("Alternate"));
   splitSpecs.AddPreDefVal(TString("Block"));

   splitSeed = 100;
   splitSpecs.DeclareOptionRef( splitSeed, "SplitSeed",
                                "Seed for random event shuffling" );

   normMode = "EqualNumEvents";  // the weight normalisation modes
   splitSpecs.DeclareOptionRef( normMode, "NormMode",
                                "Overall renormalisation of  event-by-event weights used in the training (NumEvents: average weight of 1 per event, independently for signal and background; EqualNumEvents: average weight of 1 per event for signal, and sum of weights for background equal to sum of weights for signal)" );
   splitSpecs.AddPreDefVal(TString("None"));
   splitSpecs.AddPreDefVal(TString("NumEvents"));
   splitSpecs.AddPreDefVal(TString("EqualNumEvents"));

   splitSpecs.DeclareOptionRef(fScaleWithPreselEff=kFALSE,"ScaleWithPreselEff","Scale the number of requested events by the eff. of the preselection cuts (or not)" );

   // the number of events

   // fill in the numbers
   for (UInt_t cl = 0; cl < dsi.GetNClasses(); cl++) {
      TString clName = dsi.GetClassInfo(cl)->GetName();
      TString titleTrain =  TString().Format("Number of training events of class %s (default: 0 = all)",clName.Data()).Data();
      TString titleTest  =  TString().Format("Number of test events of class %s (default: 0 = all)",clName.Data()).Data();
      TString titleSplit =  TString().Format("Split in training and test events of class %s (default: 0 = deactivated)",clName.Data()).Data();

      splitSpecs.DeclareOptionRef( nEventRequests.at(cl).nTrainingEventsRequested, TString("nTrain_")+clName, titleTrain );
      splitSpecs.DeclareOptionRef( nEventRequests.at(cl).nTestingEventsRequested , TString("nTest_")+clName , titleTest  );
      splitSpecs.DeclareOptionRef( nEventRequests.at(cl).TrainTestSplitRequested , TString("TrainTestSplit_")+clName , titleTest  );
   }

   splitSpecs.DeclareOptionRef( fVerbose, "V", "Verbosity (default: true)" );

   splitSpecs.DeclareOptionRef( fVerboseLevel=TString("Info"), "VerboseLevel", "VerboseLevel (Debug/Verbose/Info)" );
   splitSpecs.AddPreDefVal(TString("Debug"));
   splitSpecs.AddPreDefVal(TString("Verbose"));
   splitSpecs.AddPreDefVal(TString("Info"));

   fCorrelations = kTRUE;
   splitSpecs.DeclareOptionRef(fCorrelations, "Correlations", "Boolean to show correlation output (Default: true)");
   fComputeCorrelations = kTRUE;
   splitSpecs.DeclareOptionRef(fComputeCorrelations, "CalcCorrelations", "Compute correlations and also some variable statistics, e.g. min/max (Default: true )");

   splitSpecs.ParseOptions();
   splitSpecs.CheckForUnusedOptions();

   // output logging verbosity
   if (Verbose()) fLogger->SetMinType( kVERBOSE );
   if (fVerboseLevel.CompareTo("Debug")   ==0) fLogger->SetMinType( kDEBUG );
   if (fVerboseLevel.CompareTo("Verbose") ==0) fLogger->SetMinType( kVERBOSE );
   if (fVerboseLevel.CompareTo("Info")    ==0) fLogger->SetMinType( kINFO );

   // put all to upper case
   splitMode.ToUpper(); mixMode.ToUpper(); normMode.ToUpper();
   // adjust mixmode if same as splitmode option has been set
   Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
    << "\tSplitmode is: \"" << splitMode << "\" the mixmode is: \"" << mixMode << "\"" << Endl;
   if (mixMode=="SAMEASSPLITMODE") mixMode = splitMode;
   else if (mixMode!=splitMode)
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "DataSet splitmode="<<splitMode
            <<" differs from mixmode="<<mixMode<<Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// build empty event vectors
/// distributes events between kTraining/kTesting/kMaxTreeType

void
TMVA::DataSetFactory::BuildEventVector( TMVA::DataSetInfo& dsi,
                                        TMVA::DataInputHandler& dataInput,
                                        EventVectorOfClassesOfTreeType& eventsmap,
                                        EvtStatsPerClass& eventCounts)
{
   const UInt_t nclasses = dsi.GetNClasses();

   eventsmap[ Types::kTraining ]    = EventVectorOfClasses(nclasses);
   eventsmap[ Types::kTesting ]     = EventVectorOfClasses(nclasses);
   eventsmap[ Types::kMaxTreeType ] = EventVectorOfClasses(nclasses);

   // create the type, weight and boostweight branches
   const UInt_t nvars = dsi.GetNVariables();
   const UInt_t ntgts = dsi.GetNTargets();
   const UInt_t nvis  = dsi.GetNSpectators();

   for (size_t i=0; i<nclasses; i++) {
      eventCounts[i].varAvLength = new Float_t[nvars];
      for (UInt_t ivar=0; ivar<nvars; ivar++)
         eventCounts[i].varAvLength[ivar] = 0;
   }

   //Bool_t haveArrayVariable = kFALSE;
   //Bool_t *varIsArray = new Bool_t[nvars];

   // If there are NaNs in the tree:
   // => warn if used variables/cuts/weights contain nan (no problem if event is cut out)
   // => fatal if cut value is nan or (event not cut out and nans somewhere)
   // Count & collect all these warnings/errors and output them at the end.
   std::map<TString, int> nanInfWarnings;
   std::map<TString, int> nanInfErrors;

   // if we work with chains we need to remember the current tree if
   // the chain jumps to a new tree we have to reset the formulas
   for (UInt_t cl=0; cl<nclasses; cl++) {

     //Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Create training and testing trees -- looping over class \"" << dsi.GetClassInfo(cl)->GetName() << "\" ..." << Endl;

      EventStats& classEventCounts = eventCounts[cl];

      // info output for weights
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "\tWeight expression for class \'" << dsi.GetClassInfo(cl)->GetName() << "\': \""
            << dsi.GetClassInfo(cl)->GetWeight() << "\"" << Endl;

      // used for chains only
      TString currentFileName("");

      std::vector<TreeInfo>::const_iterator treeIt(dataInput.begin(dsi.GetClassInfo(cl)->GetName()));
      for (;treeIt!=dataInput.end(dsi.GetClassInfo(cl)->GetName()); ++treeIt) {

         // read first the variables
         std::vector<Float_t> vars(nvars);
         std::vector<Float_t> tgts(ntgts);
         std::vector<Float_t> vis(nvis);
         TreeInfo currentInfo = *treeIt;

         Log() << kINFO << "Building event vectors for type " << currentInfo.GetTreeType() << " " << currentInfo.GetClassName() <<  Endl;

         EventVector& event_v = eventsmap[currentInfo.GetTreeType()].at(cl);

         Bool_t isChain = (TString("TChain") == currentInfo.GetTree()->ClassName());
         currentInfo.GetTree()->LoadTree(0);
         // create the TTReeFormula to evalute later on on each single event
         ChangeToNewTree( currentInfo, dsi );

         // count number of events in tree before cut
         classEventCounts.nInitialEvents += currentInfo.GetTree()->GetEntries();

         // flag to control a warning message when size of array in disk are bigger than what requested
         Bool_t foundLargerArraySize = kFALSE;

         // loop over events in ntuple
         const UInt_t nEvts = currentInfo.GetTree()->GetEntries();
         for (Long64_t evtIdx = 0; evtIdx < nEvts; evtIdx++) {
            currentInfo.GetTree()->LoadTree(evtIdx);

            // may need to reload tree in case of chains
            if (isChain) {
               if (currentInfo.GetTree()->GetTree()->GetDirectory()->GetFile()->GetName() != currentFileName) {
                  currentFileName = currentInfo.GetTree()->GetTree()->GetDirectory()->GetFile()->GetName();
                  ChangeToNewTree( currentInfo, dsi );
               }
            }
            currentInfo.GetTree()->GetEntry(evtIdx);
            Int_t sizeOfArrays = 1;
            Int_t prevArrExpr = 0;
            Bool_t haveAllArrayData = kFALSE;

            // ======= evaluate all formulas =================

            // first we check if some of the formulas are arrays
            // This is the case when all inputs (variables, targets and spectetors are array and a TMVA event is not
            // an event of the tree but an event + array index). In this case we set the flag haveAllArrayData = true
            // Otherwise we support for arrays of variables where each
            // element of the array corresponds to a different variable like in the case of image
            // In that case the VAriableInfo has a bit, IsVariableFromArray that is set and we have a single formula for the array
            // fInputFormulaTable contains a map of the formula and the variable index to evaluate the formula
            for (UInt_t ivar = 0; ivar < nvars; ivar++) {
               // distinguish case where variable is not from an array
               if (dsi.IsVariableFromArray(ivar)) continue;
               auto inputFormula = fInputTableFormulas[ivar].first;

               Int_t ndata = inputFormula->GetNdata();

               classEventCounts.varAvLength[ivar] += ndata;
               if (ndata == 1) continue;
               haveAllArrayData = kTRUE;
               //varIsArray[ivar] = kTRUE;
               //std::cout << "Found array !!!" << std::endl;
               if (sizeOfArrays == 1) {
                  sizeOfArrays = ndata;
                  prevArrExpr = ivar;
               }
               else if (sizeOfArrays!=ndata) {
                  Log() << kERROR << Form("Dataset[%s] : ",dsi.GetName())<< "ERROR while preparing training and testing trees:" << Endl;
                  Log() << Form("Dataset[%s] : ",dsi.GetName())<< "   multiple array-type expressions of different length were encountered" << Endl;
                  Log() << Form("Dataset[%s] : ",dsi.GetName())<< "   location of error: event " << evtIdx
                        << " in tree " << currentInfo.GetTree()->GetName()
                        << " of file " << currentInfo.GetTree()->GetCurrentFile()->GetName() << Endl;
                  Log() << Form("Dataset[%s] : ",dsi.GetName())<< "   expression " << inputFormula->GetTitle() << " has "
                        << Form("Dataset[%s] : ",dsi.GetName()) << ndata << " entries, while" << Endl;
                  Log() << Form("Dataset[%s] : ",dsi.GetName())<< "   expression " << fInputTableFormulas[prevArrExpr].first->GetTitle() << " has "
                        << Form("Dataset[%s] : ",dsi.GetName())<< fInputTableFormulas[prevArrExpr].first->GetNdata() << " entries" << Endl;
                  Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName())<< "Need to abort" << Endl;
               }
            }

            // now we read the information
            for (Int_t idata = 0;  idata<sizeOfArrays; idata++) {
               Bool_t contains_NaN_or_inf = kFALSE;

               auto checkNanInf = [&](std::map<TString, int> &msgMap, Float_t value, const char *what, const char *formulaTitle) {
                  if (TMath::IsNaN(value)) {
                     contains_NaN_or_inf = kTRUE;
                     ++msgMap[TString::Format("Dataset[%s] : %s expression resolves to indeterminate value (NaN): %s", dsi.GetName(), what, formulaTitle)];
                  } else if (!TMath::Finite(value)) {
                     contains_NaN_or_inf = kTRUE;
                     ++msgMap[TString::Format("Dataset[%s] : %s expression resolves to infinite value (+inf or -inf): %s", dsi.GetName(), what, formulaTitle)];
                  }
               };

               TTreeFormula* formula = 0;

               // the cut expression
               Double_t cutVal = 1.;
               formula = fCutFormulas[cl];
               if (formula) {
                  Int_t ndata = formula->GetNdata();
                  cutVal = (ndata==1 ?
                            formula->EvalInstance(0) :
                            formula->EvalInstance(idata));
                  checkNanInf(nanInfErrors, cutVal, "Cut", formula->GetTitle());
               }

               // if event is cut out, add to warnings, else add to errors.
               auto &nanMessages = cutVal < 0.5 ? nanInfWarnings : nanInfErrors;

               // the input variable
               for (UInt_t ivar=0; ivar<nvars; ivar++) {
                  auto formulaMap = fInputTableFormulas[ivar];
                  formula = formulaMap.first;
                  int inputVarIndex = formulaMap.second;
                  // check fomula ndata size (in case of arrays variable)
                  // enough to check for ivarindex = 0 then formula is the same
                  // this check might take some time. Maybe do only in debug mode
                  if (inputVarIndex == 0 && dsi.IsVariableFromArray(ivar)) {
                     Int_t ndata = formula->GetNdata();
                     Int_t arraySize = dsi.GetVarArraySize(dsi.GetVariableInfo(ivar).GetExpression());
                     if (ndata < arraySize) {
                        Log() << kFATAL << "Size of array " << dsi.GetVariableInfo(ivar).GetExpression()
                              << " in the current tree " << currentInfo.GetTree()->GetName() << " for the event " << evtIdx
                              << " is " << ndata << " instead of " << arraySize << Endl;
                     } else if (ndata > arraySize && !foundLargerArraySize) {
                        Log() << kWARNING << "Size of array " << dsi.GetVariableInfo(ivar).GetExpression()
                              << " in the current tree " << currentInfo.GetTree()->GetName() << " for the event "
                              << evtIdx << " is " << ndata << ", larger than " << arraySize << Endl;
                        Log() << kWARNING << "Some data will then be ignored. This WARNING is printed only once, "
                              << " check in case for the other variables and events " << Endl;
                           // note that following warnings will be suppressed
                        foundLargerArraySize = kTRUE;
                     }
                  }
                  formula->SetQuickLoad(true); // is this needed ???

                  vars[ivar] =  ( !haveAllArrayData ?
                                 formula->EvalInstance(inputVarIndex) :
                                 formula->EvalInstance(idata));
                  checkNanInf(nanMessages, vars[ivar], "Input", formula->GetTitle());
               }

               // the targets
               for (UInt_t itrgt=0; itrgt<ntgts; itrgt++) {
                  formula = fTargetFormulas[itrgt];
                  Int_t ndata = formula->GetNdata();
                  tgts[itrgt] = (ndata == 1 ?
                                 formula->EvalInstance(0) :
                                 formula->EvalInstance(idata));
                  checkNanInf(nanMessages, tgts[itrgt], "Target", formula->GetTitle());
               }

               // the spectators
               for (UInt_t itVis=0; itVis<nvis; itVis++) {
                  formula = fSpectatorFormulas[itVis];
                  Int_t ndata = formula->GetNdata();
                  vis[itVis] = (ndata == 1 ?
                                formula->EvalInstance(0) :
                                formula->EvalInstance(idata));
                  checkNanInf(nanMessages, vis[itVis], "Spectator", formula->GetTitle());
               }


               // the weight
               Float_t weight = currentInfo.GetWeight(); // multiply by tree weight
               formula = fWeightFormula[cl];
               if (formula!=0) {
                  Int_t ndata = formula->GetNdata();
                  weight *= (ndata == 1 ?
                             formula->EvalInstance() :
                             formula->EvalInstance(idata));
                  checkNanInf(nanMessages, weight, "Weight", formula->GetTitle());
               }

               // Count the events before rejection due to cut or NaN
               // value (weighted and unweighted)
               classEventCounts.nEvBeforeCut++;
               if (!TMath::IsNaN(weight))
                  classEventCounts.nWeEvBeforeCut += weight;

               // apply the cut, skip rest if cut is not fulfilled
               if (cutVal<0.5) continue;

               // global flag if negative weights exist -> can be used
               // by classifiers who may require special data
               // treatment (also print warning)
               if (weight < 0) classEventCounts.nNegWeights++;

               // now read the event-values (variables and regression targets)

               if (contains_NaN_or_inf) {
                  Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName())<< "NaN or +-inf in Event " << evtIdx << Endl;
                  if (sizeOfArrays>1) Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName())<< " rejected" << Endl;
                  continue;
               }

               // Count the events after rejection due to cut or NaN value
               // (weighted and unweighted)
               classEventCounts.nEvAfterCut++;
               classEventCounts.nWeEvAfterCut += weight;

               // event accepted, fill temporary ntuple
               event_v.push_back(new Event(vars, tgts , vis, cl , weight));
            }
         }
         currentInfo.GetTree()->ResetBranchAddresses();
      }
   }

   if (!nanInfWarnings.empty()) {
      Log() << kWARNING << "Found events with NaN and/or +-inf values" << Endl;
      for (const auto &warning : nanInfWarnings) {
         auto &log = Log() << kWARNING << warning.first;
         if (warning.second > 1) log << " (" << warning.second << " times)";
         log << Endl;
      }
      Log() << kWARNING << "These NaN and/or +-infs were all removed by the specified cut, continuing." << Endl;
      Log() << Endl;
   }

   if (!nanInfErrors.empty()) {
      Log() << kWARNING << "Found events with NaN and/or +-inf values (not removed by cut)" << Endl;
      for (const auto &error : nanInfErrors) {
         auto &log = Log() << kWARNING << error.first;
         if (error.second > 1) log << " (" << error.second << " times)";
         log << Endl;
      }
      Log() << kFATAL << "How am I supposed to train a NaN or +-inf?!" << Endl;
   }

   // for output format, get the maximum class name length
   Int_t maxL = dsi.GetClassNameMaxLength();

   Log() << kHEADER << Form("[%s] : ",dsi.GetName()) << "Number of events in input trees" << Endl;
   Log() << kDEBUG << "(after possible flattening of arrays):" << Endl;


   for (UInt_t cl = 0; cl < dsi.GetNClasses(); cl++) {
      Log() << kDEBUG //<< Form("[%s] : ",dsi.GetName())
             << "    "
            << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
            << "      -- number of events       : "
            << std::setw(5) << eventCounts[cl].nEvBeforeCut
            << "  / sum of weights: " << std::setw(5) << eventCounts[cl].nWeEvBeforeCut << Endl;
   }

   for (UInt_t cl = 0; cl < dsi.GetNClasses(); cl++) {
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "    " << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
            <<" tree -- total number of entries: "
            << std::setw(5) << dataInput.GetEntries(dsi.GetClassInfo(cl)->GetName()) << Endl;
   }

   if (fScaleWithPreselEff)
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "\tPreselection: (will affect number of requested training and testing events)" << Endl;
   else
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "\tPreselection: (will NOT affect number of requested training and testing events)" << Endl;

   if (dsi.HasCuts()) {
      for (UInt_t cl = 0; cl< dsi.GetNClasses(); cl++) {
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "    " << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
               << " requirement: \"" << dsi.GetClassInfo(cl)->GetCut() << "\"" << Endl;
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "    "
               << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
               << "      -- number of events passed: "
               << std::setw(5) << eventCounts[cl].nEvAfterCut
               << "  / sum of weights: " << std::setw(5) << eventCounts[cl].nWeEvAfterCut << Endl;
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "    "
               << setiosflags(ios::left) << std::setw(maxL) << dsi.GetClassInfo(cl)->GetName()
               << "      -- efficiency             : "
               << std::setw(6) << eventCounts[cl].nWeEvAfterCut/eventCounts[cl].nWeEvBeforeCut << Endl;
      }
   }
   else Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
         << "    No preselection cuts applied on event classes" << Endl;

   //delete[] varIsArray;

}

////////////////////////////////////////////////////////////////////////////////
/// Select and distribute unassigned events to kTraining and kTesting

TMVA::DataSet*
TMVA::DataSetFactory::MixEvents( DataSetInfo& dsi,
                                 EventVectorOfClassesOfTreeType& tmpEventVector,
                                 EvtStatsPerClass& eventCounts,
                                 const TString& splitMode,
                                 const TString& mixMode,
                                 const TString& normMode,
                                 UInt_t splitSeed)
{
   TMVA::RandomGenerator<TRandom3> rndm(splitSeed);

   // ==== splitting of undefined events to kTraining and kTesting

   // if splitMode contains "RANDOM", then shuffle the undefined events
   if (splitMode.Contains( "RANDOM" ) /*&& !emptyUndefined*/ ) {
      // random shuffle the undefined events of each class
      for( UInt_t cls = 0; cls < dsi.GetNClasses(); ++cls ){
         EventVector& unspecifiedEvents = tmpEventVector[Types::kMaxTreeType].at(cls);
         if( ! unspecifiedEvents.empty() ) {
            Log() << kDEBUG << "randomly shuffling "
                  << unspecifiedEvents.size()
                  << " events of class " << cls
                  << " which are not yet associated to testing or training" << Endl;
            std::shuffle(unspecifiedEvents.begin(), unspecifiedEvents.end(), rndm);
         }
      }
   }

   // check for each class the number of training and testing events, the requested number and the available number
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "SPLITTING ========" << Endl;
   for( UInt_t cls = 0; cls < dsi.GetNClasses(); ++cls ){
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "---- class " << cls << Endl;
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "check number of training/testing events, requested and available number of events and for class " << cls << Endl;

      // check if enough or too many events are already in the training/testing eventvectors of the class cls
      EventVector& eventVectorTraining  = tmpEventVector[ Types::kTraining    ].at(cls);
      EventVector& eventVectorTesting   = tmpEventVector[ Types::kTesting     ].at(cls);
      EventVector& eventVectorUndefined = tmpEventVector[ Types::kMaxTreeType ].at(cls);

      Int_t availableTraining  = eventVectorTraining.size();
      Int_t availableTesting   = eventVectorTesting.size();
      Int_t availableUndefined = eventVectorUndefined.size();

      Float_t presel_scale;
      if (fScaleWithPreselEff) {
         presel_scale = eventCounts[cls].cutScaling();
         if (presel_scale < 1)
            Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " you have opted for scaling the number of requested training/testing events\n to be scaled by the preselection efficiency"<< Endl;
      }else{
         presel_scale = 1.; // this scaling was too confusing to most people, including me! Sorry... (Helge)
         if (eventCounts[cls].cutScaling() < 1)
            Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " you have opted for interpreting the requested number of training/testing events\n to be the number of events AFTER your preselection cuts" <<  Endl;

      }

      // If TrainTestSplit_<class> is set, set number of requested training events to split*num_all_events
      // Requested number of testing events is set to zero and therefore takes all other events
      // The option TrainTestSplit_<class> overrides nTrain_<class> or nTest_<class>
      if(eventCounts[cls].TrainTestSplitRequested < 1.0 && eventCounts[cls].TrainTestSplitRequested > 0.0){
         eventCounts[cls].nTrainingEventsRequested = Int_t(eventCounts[cls].TrainTestSplitRequested*(availableTraining+availableTesting+availableUndefined));
         eventCounts[cls].nTestingEventsRequested = Int_t(0);
      }
      else if(eventCounts[cls].TrainTestSplitRequested != 0.0) Log() << kFATAL << Form("The option TrainTestSplit_<class> has to be in range (0, 1] but is set to %f.",eventCounts[cls].TrainTestSplitRequested) << Endl;
      Int_t requestedTraining = Int_t(eventCounts[cls].nTrainingEventsRequested * presel_scale);
      Int_t requestedTesting  = Int_t(eventCounts[cls].nTestingEventsRequested  * presel_scale);

      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "events in training trees    : " << availableTraining  << Endl;
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "events in testing trees     : " << availableTesting   << Endl;
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "events in unspecified trees : " << availableUndefined << Endl;
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "requested for training      : " << requestedTraining << Endl;;

      if(presel_scale<1)
         Log() << " ( " << eventCounts[cls].nTrainingEventsRequested
               << " * " << presel_scale << " preselection efficiency)" << Endl;
      else
         Log() << Endl;
      Log() << kDEBUG << "requested for testing       : " << requestedTesting;
      if(presel_scale<1)
         Log() << " ( " << eventCounts[cls].nTestingEventsRequested
               << " * " << presel_scale << " preselection efficiency)" << Endl;
      else
         Log() << Endl;

      // nomenclature r = available training
      //              s = available testing
      //              u = available undefined
      //              R = requested training
      //              S = requested testing
      //              nR = to be used to select training events
      //              nS = to be used to select test events
      //              we have the constraint: nR + nS < r+s+u,
      //                 since we can not use more events than we have
      //              free events: Nfree = u-Thet(R-r)-Thet(S-s)
      //              nomenclature: Thet(x) = x,  if x>0 else 0
      //              nR = max(R,r) + 0.5 * Nfree
      //              nS = max(S,s) + 0.5 * Nfree
      //              nR +nS = R+S + u-R+r-S+s = u+r+s= ok! for R>r
      //              nR +nS = r+S + u-S+s = u+r+s= ok! for r>R

      // three different cases might occur here
      //
      // Case a
      // requestedTraining and requestedTesting >0
      // free events: Nfree = u-Thet(R-r)-Thet(S-s)
      //              nR = Max(R,r) + 0.5 * Nfree
      //              nS = Max(S,s) + 0.5 * Nfree
      //
      // Case b
      // exactly one of requestedTraining or requestedTesting >0
      // assume training R >0
      //    nR  = max(R,r)
      //    nS  = s+u+r-nR
      //    and  s=nS
      //
      // Case c
      // requestedTraining=0, requestedTesting=0
      // Nfree = u-|r-s|
      // if NFree >=0
      //    R = Max(r,s) + 0.5 * Nfree = S
      // else if r>s
      //    R = r; S=s+u
      // else
      //    R = r+u; S=s
      //
      // Next steps:
      // Determination of Event numbers R,S, nR, nS
      // distribute undefined events according to nR, nS
      // finally determine actual sub samples from nR and nS to be used in training / testing
      //

      Int_t useForTesting(0),useForTraining(0);
      Int_t allAvailable(availableUndefined + availableTraining + availableTesting);

      if( (requestedTraining == 0) && (requestedTesting == 0)){

         // Case C: balance the number of training and testing events

         if ( availableUndefined >= TMath::Abs(availableTraining - availableTesting) ) {
            // enough unspecified are available to equal training and testing
            useForTraining = useForTesting = allAvailable/2;
         } else {
            // all unspecified are assigned to the smaller of training / testing
            useForTraining = availableTraining;
            useForTesting  = availableTesting;
            if (availableTraining < availableTesting)
               useForTraining += availableUndefined;
            else
               useForTesting += availableUndefined;
         }
         requestedTraining = useForTraining;
         requestedTesting  = useForTesting;
      }

      else if (requestedTesting == 0){
         // case B
         useForTraining = TMath::Max(requestedTraining,availableTraining);
         if (allAvailable <  useForTraining) {
            Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName())<< "More events requested for training ("
                  << requestedTraining << ") than available ("
                  << allAvailable << ")!" << Endl;
         }
         useForTesting  = allAvailable - useForTraining; // the rest
         requestedTesting = useForTesting;
      }

      else if (requestedTraining == 0){ // case B)
         useForTesting = TMath::Max(requestedTesting,availableTesting);
         if (allAvailable <  useForTesting) {
            Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName())<< "More events requested for testing ("
                  << requestedTesting << ") than available ("
                  << allAvailable << ")!" << Endl;
         }
         useForTraining= allAvailable - useForTesting; // the rest
         requestedTraining = useForTraining;
      }

      else {
         // Case A
         // requestedTraining R and requestedTesting S >0
         // free events: Nfree = u-Thet(R-r)-Thet(S-s)
         //              nR = Max(R,r) + 0.5 * Nfree
         //              nS = Max(S,s) + 0.5 * Nfree
         Int_t stillNeedForTraining = TMath::Max(requestedTraining-availableTraining,0);
         Int_t stillNeedForTesting = TMath::Max(requestedTesting-availableTesting,0);

         int NFree = availableUndefined - stillNeedForTraining - stillNeedForTesting;
         if (NFree <0) NFree = 0;
         useForTraining = TMath::Max(requestedTraining,availableTraining) + NFree/2;
         useForTesting= allAvailable - useForTraining; // the rest
      }

      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "determined event sample size to select training sample from="<<useForTraining<<Endl;
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "determined event sample size to select test sample from="<<useForTesting<<Endl;



      // associate undefined events
      if( splitMode == "ALTERNATE" ){
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "split 'ALTERNATE'" << Endl;
         Int_t nTraining = availableTraining;
         Int_t nTesting  = availableTesting;
         for( EventVector::iterator it = eventVectorUndefined.begin(), itEnd = eventVectorUndefined.end(); it != itEnd; ){
            ++nTraining;
            if( nTraining <= requestedTraining ){
               eventVectorTraining.insert( eventVectorTraining.end(), (*it) );
               ++it;
            }
            if( it != itEnd ){
               ++nTesting;
               eventVectorTesting.insert( eventVectorTesting.end(), (*it) );
               ++it;
            }
         }
      } else {
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "split '" << splitMode << "'" << Endl;

         // test if enough events are available
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "availableundefined : " << availableUndefined << Endl;
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "useForTraining     : " << useForTraining << Endl;
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "useForTesting      : " << useForTesting  << Endl;
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "availableTraining      : " << availableTraining  << Endl;
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "availableTesting       : " << availableTesting  << Endl;

         if( availableUndefined<(useForTraining-availableTraining) ||
             availableUndefined<(useForTesting -availableTesting ) ||
             availableUndefined<(useForTraining+useForTesting-availableTraining-availableTesting ) ){
            Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName())<< "More events requested than available!" << Endl;
         }

         // select the events
         if (useForTraining>availableTraining){
            eventVectorTraining.insert(  eventVectorTraining.end() , eventVectorUndefined.begin(), eventVectorUndefined.begin()+ useForTraining- availableTraining );
            eventVectorUndefined.erase( eventVectorUndefined.begin(), eventVectorUndefined.begin() + useForTraining- availableTraining);
         }
         if (useForTesting>availableTesting){
            eventVectorTesting.insert(  eventVectorTesting.end() , eventVectorUndefined.begin(), eventVectorUndefined.begin()+ useForTesting- availableTesting );
         }
      }
      eventVectorUndefined.clear();

      // finally shorten the event vectors to the requested size by removing random events
      if (splitMode.Contains( "RANDOM" )){
         UInt_t sizeTraining  = eventVectorTraining.size();
         if( sizeTraining > UInt_t(requestedTraining) ){
            std::vector<UInt_t> indicesTraining( sizeTraining );
            // make indices
            std::generate( indicesTraining.begin(), indicesTraining.end(), TMVA::Increment<UInt_t>(0) );
            // shuffle indices
            std::shuffle(indicesTraining.begin(), indicesTraining.end(), rndm);
            // erase indices of not needed events
            indicesTraining.erase( indicesTraining.begin()+sizeTraining-UInt_t(requestedTraining), indicesTraining.end() );
            // delete all events with the given indices
            for( std::vector<UInt_t>::iterator it = indicesTraining.begin(), itEnd = indicesTraining.end(); it != itEnd; ++it ){
               delete eventVectorTraining.at( (*it) ); // delete event
               eventVectorTraining.at( (*it) ) = NULL; // set pointer to NULL
            }
            // now remove and erase all events with pointer==NULL
            eventVectorTraining.erase( std::remove( eventVectorTraining.begin(), eventVectorTraining.end(), (void*)NULL ), eventVectorTraining.end() );
         }

         UInt_t sizeTesting   = eventVectorTesting.size();
         if( sizeTesting > UInt_t(requestedTesting) ){
            std::vector<UInt_t> indicesTesting( sizeTesting );
            // make indices
            std::generate( indicesTesting.begin(), indicesTesting.end(), TMVA::Increment<UInt_t>(0) );
            // shuffle indices
            std::shuffle(indicesTesting.begin(), indicesTesting.end(), rndm);
            // erase indices of not needed events
            indicesTesting.erase( indicesTesting.begin()+sizeTesting-UInt_t(requestedTesting), indicesTesting.end() );
            // delete all events with the given indices
            for( std::vector<UInt_t>::iterator it = indicesTesting.begin(), itEnd = indicesTesting.end(); it != itEnd; ++it ){
               delete eventVectorTesting.at( (*it) ); // delete event
               eventVectorTesting.at( (*it) ) = NULL; // set pointer to NULL
            }
            // now remove and erase all events with pointer==NULL
            eventVectorTesting.erase( std::remove( eventVectorTesting.begin(), eventVectorTesting.end(), (void*)NULL ), eventVectorTesting.end() );
         }
      }
      else { // erase at end
         if( eventVectorTraining.size() < UInt_t(requestedTraining) )
            Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName())<< "DataSetFactory/requested number of training samples larger than size of eventVectorTraining.\n"
                  << "There is probably an issue. Please contact the TMVA developers." << Endl;
         std::for_each( eventVectorTraining.begin()+requestedTraining, eventVectorTraining.end(), DeleteFunctor<Event>() );
         eventVectorTraining.erase(eventVectorTraining.begin()+requestedTraining,eventVectorTraining.end());

         if( eventVectorTesting.size() < UInt_t(requestedTesting) )
            Log() << kWARNING << Form("Dataset[%s] : ",dsi.GetName())<< "DataSetFactory/requested number of testing samples larger than size of eventVectorTesting.\n"
                  << "There is probably an issue. Please contact the TMVA developers." << Endl;
         std::for_each( eventVectorTesting.begin()+requestedTesting, eventVectorTesting.end(), DeleteFunctor<Event>() );
         eventVectorTesting.erase(eventVectorTesting.begin()+requestedTesting,eventVectorTesting.end());
      }
   }

   TMVA::DataSetFactory::RenormEvents( dsi, tmpEventVector, eventCounts, normMode );

   Int_t trainingSize = 0;
   Int_t testingSize  = 0;

   // sum up number of training and testing events
   for( UInt_t cls = 0; cls < dsi.GetNClasses(); ++cls ){
      trainingSize += tmpEventVector[Types::kTraining].at(cls).size();
      testingSize  += tmpEventVector[Types::kTesting].at(cls).size();
   }

   // --- collect all training (testing) events into the training (testing) eventvector

   // create event vectors reserve enough space
   EventVector* trainingEventVector = new EventVector();
   EventVector* testingEventVector  = new EventVector();

   trainingEventVector->reserve( trainingSize );
   testingEventVector->reserve( testingSize );


   // collect the events

   // mixing of kTraining and kTesting data sets
   Log() << kDEBUG << " MIXING ============= " << Endl;

   if( mixMode == "ALTERNATE" ){
      // Inform user if he tries to use alternate mixmode for
      // event classes with different number of events, this works but the alternation stops at the last event of the smaller class
      for( UInt_t cls = 1; cls < dsi.GetNClasses(); ++cls ){
         if (tmpEventVector[Types::kTraining].at(cls).size() != tmpEventVector[Types::kTraining].at(0).size()){
            Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Training sample: You are trying to mix events in alternate mode although the classes have different event numbers. This works but the alternation stops at the last event of the smaller class."<<Endl;
         }
         if (tmpEventVector[Types::kTesting].at(cls).size() != tmpEventVector[Types::kTesting].at(0).size()){
            Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Testing sample: You are trying to mix events in alternate mode although the classes have different event numbers. This works but the alternation stops at the last event of the smaller class."<<Endl;
         }
      }
      typedef EventVector::iterator EvtVecIt;
      EvtVecIt itEvent, itEventEnd;

      // insert first class
      Log() << kDEBUG << "insert class 0 into training and test vector" << Endl;
      trainingEventVector->insert( trainingEventVector->end(), tmpEventVector[Types::kTraining].at(0).begin(), tmpEventVector[Types::kTraining].at(0).end() );
      testingEventVector->insert( testingEventVector->end(),   tmpEventVector[Types::kTesting].at(0).begin(),  tmpEventVector[Types::kTesting].at(0).end() );

      // insert other classes
      EvtVecIt itTarget;
      for( UInt_t cls = 1; cls < dsi.GetNClasses(); ++cls ){
         Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "insert class " << cls << Endl;
         // training vector
         itTarget = trainingEventVector->begin() - 1; // start one before begin
         // loop over source
         for( itEvent = tmpEventVector[Types::kTraining].at(cls).begin(), itEventEnd = tmpEventVector[Types::kTraining].at(cls).end(); itEvent != itEventEnd; ++itEvent ){
            //            if( std::distance( itTarget, trainingEventVector->end()) < Int_t(cls+1) ) {
            if( (trainingEventVector->end() - itTarget) < Int_t(cls+1) ) {
               itTarget = trainingEventVector->end();
               trainingEventVector->insert( itTarget, itEvent, itEventEnd ); // fill in the rest without mixing
               break;
            }else{
               itTarget += cls+1;
               trainingEventVector->insert( itTarget, (*itEvent) ); // fill event
            }
         }
         // testing vector
         itTarget = testingEventVector->begin() - 1;
         // loop over source
         for( itEvent = tmpEventVector[Types::kTesting].at(cls).begin(), itEventEnd = tmpEventVector[Types::kTesting].at(cls).end(); itEvent != itEventEnd; ++itEvent ){
            //             if( std::distance( itTarget, testingEventVector->end()) < Int_t(cls+1) ) {
            if( ( testingEventVector->end() - itTarget ) < Int_t(cls+1) ) {
               itTarget = testingEventVector->end();
               testingEventVector->insert( itTarget, itEvent, itEventEnd ); // fill in the rest without mixing
               break;
            }else{
               itTarget += cls+1;
               testingEventVector->insert( itTarget, (*itEvent) ); // fill event
            }
         }
      }
   }else{
      for( UInt_t cls = 0; cls < dsi.GetNClasses(); ++cls ){
         trainingEventVector->insert( trainingEventVector->end(), tmpEventVector[Types::kTraining].at(cls).begin(), tmpEventVector[Types::kTraining].at(cls).end() );
         testingEventVector->insert ( testingEventVector->end(),  tmpEventVector[Types::kTesting].at(cls).begin(),  tmpEventVector[Types::kTesting].at(cls).end()  );
      }
   }
   // delete the tmpEventVector (but not the events therein)
   tmpEventVector[Types::kTraining].clear();
   tmpEventVector[Types::kTesting].clear();

   tmpEventVector[Types::kMaxTreeType].clear();

   if (mixMode == "RANDOM") {
      Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "shuffling events"<<Endl;

      std::shuffle(trainingEventVector->begin(), trainingEventVector->end(), rndm);
      std::shuffle(testingEventVector->begin(),  testingEventVector->end(),  rndm);
   }

   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "trainingEventVector " << trainingEventVector->size() << Endl;
   Log() << kDEBUG << Form("Dataset[%s] : ",dsi.GetName())<< "testingEventVector  " << testingEventVector->size() << Endl;

   // create dataset
   DataSet* ds = new DataSet(dsi);

   // Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Create internal training tree" << Endl;
   ds->SetEventCollection(trainingEventVector, Types::kTraining );
   // Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Create internal testing tree" << Endl;
   ds->SetEventCollection(testingEventVector,  Types::kTesting  );


   if (ds->GetNTrainingEvents() < 1){
      Log() << kFATAL << "Dataset " << std::string(dsi.GetName()) << " does not have any training events, I better stop here and let you fix that one first " << Endl;
   }

   if (ds->GetNTestEvents() < 1) {
      Log() << kERROR << "Dataset " << std::string(dsi.GetName()) << " does not have any testing events, guess that will cause problems later..but for now, I continue " << Endl;
   }

   delete trainingEventVector;
   delete testingEventVector;
   return ds;

}

////////////////////////////////////////////////////////////////////////////////
/// renormalisation of the TRAINING event weights
///  - none       (kind of obvious) .. use the weights as supplied by the
///                user..  (we store however the relative weight for later use)
///  - numEvents
///  - equalNumEvents reweight the training events such that the sum of all
///                   backgr. (class > 0) weights equal that of the signal (class 0)

void
TMVA::DataSetFactory::RenormEvents( TMVA::DataSetInfo& dsi,
                                    EventVectorOfClassesOfTreeType& tmpEventVector,
                                    const EvtStatsPerClass& eventCounts,
                                    const TString& normMode )
{


   // print rescaling info
   // ---------------------------------
   // compute sizes and sums of weights
   Int_t trainingSize = 0;
   Int_t testingSize  = 0;

   ValuePerClass trainingSumWeightsPerClass( dsi.GetNClasses() );
   ValuePerClass testingSumWeightsPerClass( dsi.GetNClasses() );

   NumberPerClass trainingSizePerClass( dsi.GetNClasses() );
   NumberPerClass testingSizePerClass( dsi.GetNClasses() );

   Double_t trainingSumSignalWeights = 0;
   Double_t trainingSumBackgrWeights = 0; // Backgr. includes all classes that are not signal
   Double_t testingSumSignalWeights  = 0;
   Double_t testingSumBackgrWeights  = 0; // Backgr. includes all classes that are not signal



   for( UInt_t cls = 0, clsEnd = dsi.GetNClasses(); cls < clsEnd; ++cls ){
      trainingSizePerClass.at(cls) = tmpEventVector[Types::kTraining].at(cls).size();
      testingSizePerClass.at(cls)  = tmpEventVector[Types::kTesting].at(cls).size();

      trainingSize += trainingSizePerClass.back();
      testingSize  += testingSizePerClass.back();

      // the functional solution
      // sum up the weights in Double_t although the individual weights are Float_t to prevent rounding issues in addition of floating points
      //
      // accumulate --> does what the name says
      //     begin() and end() denote the range of the vector to be accumulated
      //     Double_t(0) tells accumulate the type and the starting value
      //     compose_binary creates a BinaryFunction of ...
      //         std::plus<Double_t>() knows how to sum up two doubles
      //         null<Double_t>() leaves the first argument (the running sum) unchanged and returns it
      //
      // all together sums up all the event-weights of the events in the vector and returns it
      trainingSumWeightsPerClass.at(cls) =
         std::accumulate(tmpEventVector[Types::kTraining].at(cls).begin(),
                         tmpEventVector[Types::kTraining].at(cls).end(),
                         Double_t(0), [](Double_t w, const TMVA::Event *E) { return w + E->GetOriginalWeight(); });

      testingSumWeightsPerClass.at(cls) =
         std::accumulate(tmpEventVector[Types::kTesting].at(cls).begin(),
                         tmpEventVector[Types::kTesting].at(cls).end(),
                         Double_t(0), [](Double_t w, const TMVA::Event *E) { return w + E->GetOriginalWeight(); });

      if ( cls == dsi.GetSignalClassIndex()){
         trainingSumSignalWeights += trainingSumWeightsPerClass.at(cls);
         testingSumSignalWeights  += testingSumWeightsPerClass.at(cls);
      }else{
         trainingSumBackgrWeights += trainingSumWeightsPerClass.at(cls);
         testingSumBackgrWeights  += testingSumWeightsPerClass.at(cls);
      }
   }

   // ---------------------------------
   // compute renormalization factors

   ValuePerClass renormFactor( dsi.GetNClasses() );


   // for information purposes
   dsi.SetNormalization( normMode );
   // !! these will be overwritten later by the 'rescaled' ones if
   //    NormMode != None  !!!
   dsi.SetTrainingSumSignalWeights(trainingSumSignalWeights);
   dsi.SetTrainingSumBackgrWeights(trainingSumBackgrWeights);
   dsi.SetTestingSumSignalWeights(testingSumSignalWeights);
   dsi.SetTestingSumBackgrWeights(testingSumBackgrWeights);


   if (normMode == "NONE") {
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "No weight renormalisation applied: use original global and event weights" << Endl;
      return;
   }
   //changed by Helge 27.5.2013     What on earth was done here before? I still remember the idea behind this which apparently was
   //NOT understood by the 'programmer' :)  .. the idea was to have SAME amount of effective TRAINING data for signal and background.
   // Testing events are totally irrelevant for this and might actually skew the whole normalisation!!
   else if (normMode == "NUMEVENTS") {
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "\tWeight renormalisation mode: \"NumEvents\": renormalises all event classes " << Endl;
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << " such that the effective (weighted) number of events in each class equals the respective " << Endl;
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << " number of events (entries) that you demanded in PrepareTrainingAndTestTree(\"\",\"nTrain_Signal=.. )" << Endl;
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << " ... i.e. such that Sum[i=1..N_j]{w_i} = N_j, j=0,1,2..." << Endl;
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << " ... (note that N_j is the sum of TRAINING events (nTrain_j...with j=Signal,Background.." << Endl;
      Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << " ..... Testing events are not renormalised nor included in the renormalisation factor! )"<< Endl;

      for( UInt_t cls = 0, clsEnd = dsi.GetNClasses(); cls < clsEnd; ++cls ){
         //         renormFactor.at(cls) = ( (trainingSizePerClass.at(cls) + testingSizePerClass.at(cls))/
         //                                  (trainingSumWeightsPerClass.at(cls) + testingSumWeightsPerClass.at(cls)) );
         //changed by Helge 27.5.2013
         renormFactor.at(cls) = ((Float_t)trainingSizePerClass.at(cls) )/
            (trainingSumWeightsPerClass.at(cls)) ;
      }
   }
   else if (normMode == "EQUALNUMEVENTS") {
      //changed by Helge 27.5.2013     What on earth was done here before? I still remember the idea behind this which apparently was
      //NOT understood by the 'programmer' :)  .. the idea was to have SAME amount of effective TRAINING data for signal and background.
      //done here was something like having each data source normalized to its number of entries and this even for training+testing together.
      // what should this have been good for ???

      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "Weight renormalisation mode: \"EqualNumEvents\": renormalises all event classes ..." << Endl;
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " such that the effective (weighted) number of events in each class is the same " << Endl;
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " (and equals the number of events (entries) given for class=0 )" << Endl;
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "... i.e. such that Sum[i=1..N_j]{w_i} = N_classA, j=classA, classB, ..." << Endl;
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << "... (note that N_j is the sum of TRAINING events" << Endl;
      Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << " ..... Testing events are not renormalised nor included in the renormalisation factor!)" << Endl;

      // normalize to size of first class
      UInt_t referenceClass = 0;
      for (UInt_t cls = 0, clsEnd = dsi.GetNClasses(); cls < clsEnd; ++cls ) {
         renormFactor.at(cls) = Float_t(trainingSizePerClass.at(referenceClass))/
            (trainingSumWeightsPerClass.at(cls));
      }
   }
   else {
      Log() << kFATAL << Form("Dataset[%s] : ",dsi.GetName())<< "<PrepareForTrainingAndTesting> Unknown NormMode: " << normMode << Endl;
   }

   // ---------------------------------
   // now apply the normalization factors
   Int_t maxL = dsi.GetClassNameMaxLength();
   for (UInt_t cls = 0, clsEnd = dsi.GetNClasses(); cls<clsEnd; ++cls) {
     Log() << kDEBUG //<< Form("Dataset[%s] : ",dsi.GetName())
       << "--> Rescale " << setiosflags(ios::left) << std::setw(maxL)
            << dsi.GetClassInfo(cls)->GetName() << " event weights by factor: " << renormFactor.at(cls) << Endl;
      for (EventVector::iterator it = tmpEventVector[Types::kTraining].at(cls).begin(),
              itEnd = tmpEventVector[Types::kTraining].at(cls).end(); it != itEnd; ++it){
         (*it)->SetWeight ((*it)->GetWeight() * renormFactor.at(cls));
      }

   }


   // print out the result
   // (same code as before --> this can be done nicer )
   //

   Log() << kINFO //<< Form("Dataset[%s] : ",dsi.GetName())
    << "Number of training and testing events" << Endl;
   Log() << kDEBUG << "\tafter rescaling:" << Endl;
   Log() << kINFO //<< Form("Dataset[%s] : ",dsi.GetName())
    << "---------------------------------------------------------------------------" << Endl;

   trainingSumSignalWeights = 0;
   trainingSumBackgrWeights = 0; // Backgr. includes all classes that are not signal
   testingSumSignalWeights  = 0;
   testingSumBackgrWeights  = 0; // Backgr. includes all classes that are not signal

   for( UInt_t cls = 0, clsEnd = dsi.GetNClasses(); cls < clsEnd; ++cls ){
      trainingSumWeightsPerClass.at(cls) =
         std::accumulate(tmpEventVector[Types::kTraining].at(cls).begin(),
                         tmpEventVector[Types::kTraining].at(cls).end(),
                         Double_t(0), [](Double_t w, const TMVA::Event *E) { return w + E->GetOriginalWeight(); });

      testingSumWeightsPerClass.at(cls) =
         std::accumulate(tmpEventVector[Types::kTesting].at(cls).begin(),
                         tmpEventVector[Types::kTesting].at(cls).end(),
                         Double_t(0), [](Double_t w, const TMVA::Event *E) { return w + E->GetOriginalWeight(); });

      if ( cls == dsi.GetSignalClassIndex()){
         trainingSumSignalWeights += trainingSumWeightsPerClass.at(cls);
         testingSumSignalWeights  += testingSumWeightsPerClass.at(cls);
      }else{
         trainingSumBackgrWeights += trainingSumWeightsPerClass.at(cls);
         testingSumBackgrWeights  += testingSumWeightsPerClass.at(cls);
      }

      // output statistics

      Log() << kINFO //<< Form("Dataset[%s] : ",dsi.GetName())
       << setiosflags(ios::left) << std::setw(maxL)
            << dsi.GetClassInfo(cls)->GetName() << " -- "
            << "training events            : " << trainingSizePerClass.at(cls) << Endl;
      Log() << kDEBUG << "\t(sum of weights: " << trainingSumWeightsPerClass.at(cls) << ")"
            <<  " - requested were " << eventCounts[cls].nTrainingEventsRequested << " events" << Endl;
      Log() << kINFO //<< Form("Dataset[%s] : ",dsi.GetName())
       << setiosflags(ios::left) << std::setw(maxL)
            << dsi.GetClassInfo(cls)->GetName() << " -- "
            << "testing events             : " << testingSizePerClass.at(cls) << Endl;
      Log() << kDEBUG << "\t(sum of weights: " << testingSumWeightsPerClass.at(cls) << ")"
            <<  " - requested were " << eventCounts[cls].nTestingEventsRequested << " events" << Endl;
      Log() << kINFO //<< Form("Dataset[%s] : ",dsi.GetName())
       << setiosflags(ios::left) << std::setw(maxL)
            << dsi.GetClassInfo(cls)->GetName() << " -- "
            << "training and testing events: "
            << (trainingSizePerClass.at(cls)+testingSizePerClass.at(cls)) << Endl;
      Log() << kDEBUG << "\t(sum of weights: "
            << (trainingSumWeightsPerClass.at(cls)+testingSumWeightsPerClass.at(cls)) << ")" << Endl;
      if(eventCounts[cls].nEvAfterCut<eventCounts[cls].nEvBeforeCut) {
         Log() << kINFO << Form("Dataset[%s] : ",dsi.GetName()) << setiosflags(ios::left) << std::setw(maxL)
               << dsi.GetClassInfo(cls)->GetName() << " -- "
               << "due to the preselection a scaling factor has been applied to the numbers of requested events: "
               << eventCounts[cls].cutScaling() << Endl;
      }
   }
   Log() << kINFO << Endl;

   // for information purposes
   dsi.SetTrainingSumSignalWeights(trainingSumSignalWeights);
   dsi.SetTrainingSumBackgrWeights(trainingSumBackgrWeights);
   dsi.SetTestingSumSignalWeights(testingSumSignalWeights);
   dsi.SetTestingSumBackgrWeights(testingSumBackgrWeights);


}
