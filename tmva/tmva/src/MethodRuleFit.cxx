// @(#)root/tmva $Id$
// Author: Fredrik Tegenfeldt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodRuleFit                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch>  - Iowa State U., USA     *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U.                                                             *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodRuleFit
\ingroup TMVA
J Friedman's RuleFit method
*/

#include "TMVA/MethodRuleFit.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/CrossEntropy.h"
#include "TMVA/DataSet.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/GiniIndex.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MisClassificationError.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Ranking.h"
#include "TMVA/RuleFitAPI.h"
#include "TMVA/SdivSqrtSplusB.h"
#include "TMVA/SeparationBase.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TRandom3.h"
#include "TMatrix.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <list>
#include <random>

using std::min;

REGISTER_METHOD(RuleFit)

ClassImp(TMVA::MethodRuleFit);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodRuleFit::MethodRuleFit( const TString& jobName,
                                       const TString& methodTitle,
                                       DataSetInfo& theData,
                                       const TString& theOption) :
   MethodBase( jobName, Types::kRuleFit, methodTitle, theData, theOption)
   , fSignalFraction(0)
   , fNTImportance(0)
   , fNTCoefficient(0)
   , fNTSupport(0)
   , fNTNcuts(0)
   , fNTNvars(0)
   , fNTPtag(0)
   , fNTPss(0)
   , fNTPsb(0)
   , fNTPbs(0)
   , fNTPbb(0)
   , fNTSSB(0)
   , fNTType(0)
   , fUseRuleFitJF(kFALSE)
   , fRFNrules(0)
   , fRFNendnodes(0)
   , fNTrees(0)
   , fTreeEveFrac(0)
   , fSepType(0)
   , fMinFracNEve(0)
   , fMaxFracNEve(0)
   , fNCuts(0)
   , fPruneMethod(TMVA::DecisionTree::kCostComplexityPruning)
   , fPruneStrength(0)
   , fUseBoost(kFALSE)
   , fGDPathEveFrac(0)
   , fGDValidEveFrac(0)
   , fGDTau(0)
   , fGDTauPrec(0)
   , fGDTauMin(0)
   , fGDTauMax(0)
   , fGDTauScan(0)
   , fGDPathStep(0)
   , fGDNPathSteps(0)
   , fGDErrScale(0)
   , fMinimp(0)
   , fRuleMinDist(0)
   , fLinQuantile(0)
{
      fMonitorNtuple = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodRuleFit::MethodRuleFit( DataSetInfo& theData,
                                    const TString& theWeightFile) :
   MethodBase( Types::kRuleFit, theData, theWeightFile)
   , fSignalFraction(0)
   , fNTImportance(0)
   , fNTCoefficient(0)
   , fNTSupport(0)
   , fNTNcuts(0)
   , fNTNvars(0)
   , fNTPtag(0)
   , fNTPss(0)
   , fNTPsb(0)
   , fNTPbs(0)
   , fNTPbb(0)
   , fNTSSB(0)
   , fNTType(0)
   , fUseRuleFitJF(kFALSE)
   , fRFNrules(0)
   , fRFNendnodes(0)
   , fNTrees(0)
   , fTreeEveFrac(0)
   , fSepType(0)
   , fMinFracNEve(0)
   , fMaxFracNEve(0)
   , fNCuts(0)
   , fPruneMethod(TMVA::DecisionTree::kCostComplexityPruning)
   , fPruneStrength(0)
   , fUseBoost(kFALSE)
   , fGDPathEveFrac(0)
   , fGDValidEveFrac(0)
   , fGDTau(0)
   , fGDTauPrec(0)
   , fGDTauMin(0)
   , fGDTauMax(0)
   , fGDTauScan(0)
   , fGDPathStep(0)
   , fGDNPathSteps(0)
   , fGDErrScale(0)
   , fMinimp(0)
   , fRuleMinDist(0)
   , fLinQuantile(0)
{
      fMonitorNtuple = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodRuleFit::~MethodRuleFit( void )
{
   for (UInt_t i=0; i<fEventSample.size(); i++) delete fEventSample[i];
   for (UInt_t i=0; i<fForest.size(); i++)      delete fForest[i];
}

////////////////////////////////////////////////////////////////////////////////
/// RuleFit can handle classification with 2 classes

Bool_t TMVA::MethodRuleFit::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
/// know options.
///
/// #### general
///
///  - RuleFitModule  <string>
///    available values are:
///    - RFTMVA      - use TMVA implementation
///    - RFFriedman  - use Friedmans original implementation
///
/// #### Path search (fitting)
///
///  - GDTau          <float>      gradient-directed path: fit threshold, default
///  - GDTauPrec      <float>      gradient-directed path: precision of estimated tau
///  - GDStep         <float>      gradient-directed path: step size
///  - GDNSteps       <float>      gradient-directed path: number of steps
///  - GDErrScale     <float>      stop scan when error>scale*errmin
///
/// #### Tree generation
///
///  - fEventsMin     <float>      minimum fraction of events in a splittable node
///  - fEventsMax     <float>      maximum fraction of events in a splittable node
///  - nTrees         <float>      number of trees in forest.
///  - ForestType     <string>
///    available values are:
///    - Random    - create forest using random subsample and only random variables subset at each node
///    - AdaBoost  - create forest with boosted events
///
/// #### Model creation
///
///  - RuleMinDist    <float>      min distance allowed between rules
///  - MinImp         <float>      minimum rule importance accepted
///  - Model          <string>     model to be used
///    available values are:
///    - ModRuleLinear <default>
///    - ModRule
///    - ModLinear
///
/// #### Friedmans module
///
///  - RFWorkDir      <string>     directory where Friedmans module (rf_go.exe) is installed
///  - RFNrules       <int>        maximum number of rules allowed
///  - RFNendnodes    <int>        average number of end nodes in the forest of trees

void TMVA::MethodRuleFit::DeclareOptions()
{
   DeclareOptionRef(fGDTau=-1,             "GDTau",          "Gradient-directed (GD) path: default fit cut-off");
   DeclareOptionRef(fGDTauPrec=0.01,       "GDTauPrec",      "GD path: precision of tau");
   DeclareOptionRef(fGDPathStep=0.01,      "GDStep",         "GD path: step size");
   DeclareOptionRef(fGDNPathSteps=10000,   "GDNSteps",       "GD path: number of steps");
   DeclareOptionRef(fGDErrScale=1.1,       "GDErrScale",     "Stop scan when error > scale*errmin");
   DeclareOptionRef(fLinQuantile,           "LinQuantile",  "Quantile of linear terms (removes outliers)");
   DeclareOptionRef(fGDPathEveFrac=0.5,    "GDPathEveFrac",  "Fraction of events used for the path search");
   DeclareOptionRef(fGDValidEveFrac=0.5,   "GDValidEveFrac", "Fraction of events used for the validation");
   // tree options
   DeclareOptionRef(fMinFracNEve=0.1,      "fEventsMin",     "Minimum fraction of events in a splittable node");
   DeclareOptionRef(fMaxFracNEve=0.9,      "fEventsMax",     "Maximum fraction of events in a splittable node");
   DeclareOptionRef(fNTrees=20,            "nTrees",         "Number of trees in forest.");

   DeclareOptionRef(fForestTypeS="AdaBoost",  "ForestType",   "Method to use for forest generation (AdaBoost or RandomForest)");
   AddPreDefVal(TString("AdaBoost"));
   AddPreDefVal(TString("Random"));
   // rule cleanup options
   DeclareOptionRef(fRuleMinDist=0.001,    "RuleMinDist",    "Minimum distance between rules");
   DeclareOptionRef(fMinimp=0.01,          "MinImp",         "Minimum rule importance accepted");
   // rule model option
   DeclareOptionRef(fModelTypeS="ModRuleLinear", "Model",    "Model to be used");
   AddPreDefVal(TString("ModRule"));
   AddPreDefVal(TString("ModRuleLinear"));
   AddPreDefVal(TString("ModLinear"));
   DeclareOptionRef(fRuleFitModuleS="RFTMVA",  "RuleFitModule","Which RuleFit module to use");
   AddPreDefVal(TString("RFTMVA"));
   AddPreDefVal(TString("RFFriedman"));

   DeclareOptionRef(fRFWorkDir="./rulefit", "RFWorkDir",    "Friedman\'s RuleFit module (RFF): working dir");
   DeclareOptionRef(fRFNrules=2000,         "RFNrules",     "RFF: Mximum number of rules");
   DeclareOptionRef(fRFNendnodes=4,         "RFNendnodes",  "RFF: Average number of end nodes");
}

////////////////////////////////////////////////////////////////////////////////
/// process the options specified by the user

void TMVA::MethodRuleFit::ProcessOptions()
{
   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not yet available for method: "
            << GetMethodTypeName()
            << " --> please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }

   fRuleFitModuleS.ToLower();
   if      (fRuleFitModuleS == "rftmva")     fUseRuleFitJF = kFALSE;
   else if (fRuleFitModuleS == "rffriedman") fUseRuleFitJF = kTRUE;
   else                                      fUseRuleFitJF = kTRUE;

   fSepTypeS.ToLower();
   if      (fSepTypeS == "misclassificationerror") fSepType = new MisClassificationError();
   else if (fSepTypeS == "giniindex")              fSepType = new GiniIndex();
   else if (fSepTypeS == "crossentropy")           fSepType = new CrossEntropy();
   else                                            fSepType = new SdivSqrtSplusB();

   fModelTypeS.ToLower();
   if      (fModelTypeS == "modlinear" ) fRuleFit.SetModelLinear();
   else if (fModelTypeS == "modrule" )   fRuleFit.SetModelRules();
   else                                  fRuleFit.SetModelFull();

   fPruneMethodS.ToLower();
   if      (fPruneMethodS == "expectederror" )   fPruneMethod  = DecisionTree::kExpectedErrorPruning;
   else if (fPruneMethodS == "costcomplexity" )  fPruneMethod  = DecisionTree::kCostComplexityPruning;
   else                                          fPruneMethod  = DecisionTree::kNoPruning;

   fForestTypeS.ToLower();
   if      (fForestTypeS == "random" )   fUseBoost = kFALSE;
   else if (fForestTypeS == "adaboost" ) fUseBoost = kTRUE;
   else                                  fUseBoost = kTRUE;
   //
   // if creating the forest by boosting the events
   // the full training sample is used per tree
   // -> only true for the TMVA version of RuleFit.
   if (fUseBoost && (!fUseRuleFitJF)) fTreeEveFrac = 1.0;

   // check event fraction for tree generation
   // if <0 set to automatic number
   if (fTreeEveFrac<=0) {
      Int_t nevents = Data()->GetNTrainingEvents();
      Double_t n = static_cast<Double_t>(nevents);
      fTreeEveFrac = min( 0.5, (100.0 +6.0*sqrt(n))/n);
   }
   // verify ranges of options
   VerifyRange(Log(), "nTrees",        fNTrees,0,100000,20);
   VerifyRange(Log(), "MinImp",        fMinimp,0.0,1.0,0.0);
   VerifyRange(Log(), "GDTauPrec",     fGDTauPrec,1e-5,5e-1);
   VerifyRange(Log(), "GDTauMin",      fGDTauMin,0.0,1.0);
   VerifyRange(Log(), "GDTauMax",      fGDTauMax,fGDTauMin,1.0);
   VerifyRange(Log(), "GDPathStep",    fGDPathStep,0.0,100.0,0.01);
   VerifyRange(Log(), "GDErrScale",    fGDErrScale,1.0,100.0,1.1);
   VerifyRange(Log(), "GDPathEveFrac", fGDPathEveFrac,0.01,0.9,0.5);
   VerifyRange(Log(), "GDValidEveFrac",fGDValidEveFrac,0.01,1.0-fGDPathEveFrac,1.0-fGDPathEveFrac);
   VerifyRange(Log(), "fEventsMin",    fMinFracNEve,0.0,1.0);
   VerifyRange(Log(), "fEventsMax",    fMaxFracNEve,fMinFracNEve,1.0);

   fRuleFit.GetRuleEnsemblePtr()->SetLinQuantile(fLinQuantile);
   fRuleFit.GetRuleFitParamsPtr()->SetGDTauRange(fGDTauMin,fGDTauMax);
   fRuleFit.GetRuleFitParamsPtr()->SetGDTau(fGDTau);
   fRuleFit.GetRuleFitParamsPtr()->SetGDTauPrec(fGDTauPrec);
   fRuleFit.GetRuleFitParamsPtr()->SetGDTauScan(fGDTauScan);
   fRuleFit.GetRuleFitParamsPtr()->SetGDPathStep(fGDPathStep);
   fRuleFit.GetRuleFitParamsPtr()->SetGDNPathSteps(fGDNPathSteps);
   fRuleFit.GetRuleFitParamsPtr()->SetGDErrScale(fGDErrScale);
   fRuleFit.SetImportanceCut(fMinimp);
   fRuleFit.SetRuleMinDist(fRuleMinDist);


   // check if Friedmans module is used.
   // print a message concerning the options.
   if (fUseRuleFitJF) {
      Log() << kINFO << "" << Endl;
      Log() << kINFO << "--------------------------------------" <<Endl;
      Log() << kINFO << "Friedmans RuleFit module is selected." << Endl;
      Log() << kINFO << "Only the following options are used:" << Endl;
      Log() << kINFO <<  Endl;
      Log() << kINFO << gTools().Color("bold") << "   Model"        << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   RFWorkDir"    << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   RFNrules"     << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   RFNendnodes"  << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   GDNPathSteps" << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   GDPathStep"   << gTools().Color("reset") << Endl;
      Log() << kINFO << gTools().Color("bold") << "   GDErrScale"   << gTools().Color("reset") << Endl;
      Log() << kINFO << "--------------------------------------" <<Endl;
      Log() << kINFO << Endl;
   }

   // Select what weight to use in the 'importance' rule visualisation plots.
   // Note that if UseCoefficientsVisHists() is selected, the following weight is used:
   //    w = rule coefficient * rule support
   // The support is a positive number which is 0 if no events are accepted by the rule.
   // Normally the importance gives more useful information.
   //
   //fRuleFit.UseCoefficientsVisHists();
   fRuleFit.UseImportanceVisHists();

   fRuleFit.SetMsgType( Log().GetMinType() );

   if (HasTrainingTree()) InitEventSample();

}

////////////////////////////////////////////////////////////////////////////////
/// initialize the monitoring ntuple

void TMVA::MethodRuleFit::InitMonitorNtuple()
{
   BaseDir()->cd();
   fMonitorNtuple= new TTree("MonitorNtuple_RuleFit","RuleFit variables");
   fMonitorNtuple->Branch("importance",&fNTImportance,"importance/D");
   fMonitorNtuple->Branch("support",&fNTSupport,"support/D");
   fMonitorNtuple->Branch("coefficient",&fNTCoefficient,"coefficient/D");
   fMonitorNtuple->Branch("ncuts",&fNTNcuts,"ncuts/I");
   fMonitorNtuple->Branch("nvars",&fNTNvars,"nvars/I");
   fMonitorNtuple->Branch("type",&fNTType,"type/I");
   fMonitorNtuple->Branch("ptag",&fNTPtag,"ptag/D");
   fMonitorNtuple->Branch("pss",&fNTPss,"pss/D");
   fMonitorNtuple->Branch("psb",&fNTPsb,"psb/D");
   fMonitorNtuple->Branch("pbs",&fNTPbs,"pbs/D");
   fMonitorNtuple->Branch("pbb",&fNTPbb,"pbb/D");
   fMonitorNtuple->Branch("soversb",&fNTSSB,"soversb/D");
}

////////////////////////////////////////////////////////////////////////////////
/// default initialization

void TMVA::MethodRuleFit::Init()
{
   // the minimum requirement to declare an event signal-like
   SetSignalReferenceCut( 0.0 );

   // set variables that used to be options
   // any modifications are then made in ProcessOptions()
   fLinQuantile   = 0.025;       // Quantile of linear terms (remove outliers)
   fTreeEveFrac   = -1.0;        // Fraction of events used to train each tree
   fNCuts         = 20;          // Number of steps during node cut optimisation
   fSepTypeS      = "GiniIndex"; // Separation criterion for node splitting; see BDT
   fPruneMethodS  = "NONE";      // Pruning method; see BDT
   fPruneStrength = 3.5;         // Pruning strength; see BDT
   fGDTauMin      = 0.0;         // Gradient-directed path: min fit threshold (tau)
   fGDTauMax      = 1.0;         // Gradient-directed path: max fit threshold (tau)
   fGDTauScan     = 1000;        // Gradient-directed path: number of points scanning for best tau

}

////////////////////////////////////////////////////////////////////////////////
/// write all Events from the Tree into a vector of Events, that are
/// more easily manipulated.
/// This method should never be called without existing trainingTree, as it
/// the vector of events from the ROOT training tree

void TMVA::MethodRuleFit::InitEventSample( void )
{
   if (Data()->GetNEvents()==0) Log() << kFATAL << "<Init> Data().TrainingTree() is zero pointer" << Endl;

   Int_t nevents = Data()->GetNEvents();
   for (Int_t ievt=0; ievt<nevents; ievt++){
      const Event * ev = GetEvent(ievt);
      fEventSample.push_back( new Event(*ev));
   }
   if (fTreeEveFrac<=0) {
      Double_t n = static_cast<Double_t>(nevents);
      fTreeEveFrac = min( 0.5, (100.0 +6.0*sqrt(n))/n);
   }
   if (fTreeEveFrac>1.0) fTreeEveFrac=1.0;
   //
   std::shuffle(fEventSample.begin(), fEventSample.end(), std::default_random_engine{});
   //
   Log() << kDEBUG << "Set sub-sample fraction to " << fTreeEveFrac << Endl;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodRuleFit::Train( void )
{
   TMVA::DecisionTreeNode::SetIsTraining(true);
   // training of rules

  if(!IsSilentFile()) InitMonitorNtuple();

   // fill the STL Vector with the event sample
   this->InitEventSample();

   if (fUseRuleFitJF) {
      TrainJFRuleFit();
   }
   else {
      TrainTMVARuleFit();
   }
   fRuleFit.GetRuleEnsemblePtr()->ClearRuleMap();
   TMVA::DecisionTreeNode::SetIsTraining(false);
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// training of rules using TMVA implementation

void TMVA::MethodRuleFit::TrainTMVARuleFit( void )
{
   if (IsNormalised()) Log() << kFATAL << "\"Normalise\" option cannot be used with RuleFit; "
                             << "please remove the option from the configuration string, or "
                             << "use \"!Normalise\""
                             << Endl;

   // timer
   Timer timer( 1, GetName() );

   // test tree nmin cut -> for debug purposes
   // the routine will generate trees with stopping cut on N(eve) given by
   // a fraction between [20,N(eve)-1].
   //
   //   MakeForestRnd();
   //   exit(1);
   //

   // Init RuleFit object and create rule ensemble
   // + make forest & rules
   fRuleFit.Initialize( this );

   // Make forest of decision trees
   //   if (fRuleFit.GetRuleEnsemble().DoRules()) fRuleFit.MakeForest();

   // Fit the rules
   Log() << kDEBUG << "Fitting rule coefficients ..." << Endl;
   fRuleFit.FitCoefficients();

   // Calculate importance
   Log() << kDEBUG << "Computing rule and variable importance" << Endl;
   fRuleFit.CalcImportance();

   // Output results and fill monitor ntuple
   fRuleFit.GetRuleEnsemblePtr()->Print();
   //
   if(!IsSilentFile())
    {
        Log() << kDEBUG << "Filling rule ntuple" << Endl;
        UInt_t nrules = fRuleFit.GetRuleEnsemble().GetRulesConst().size();
        const Rule *rule;
        for (UInt_t i=0; i<nrules; i++ ) {
            rule            = fRuleFit.GetRuleEnsemble().GetRulesConst(i);
            fNTImportance   = rule->GetRelImportance();
            fNTSupport      = rule->GetSupport();
            fNTCoefficient  = rule->GetCoefficient();
            fNTType         = (rule->IsSignalRule() ? 1:-1 );
            fNTNvars        = rule->GetRuleCut()->GetNvars();
            fNTNcuts        = rule->GetRuleCut()->GetNcuts();
            fNTPtag         = fRuleFit.GetRuleEnsemble().GetRulePTag(i); // should be identical with support
            fNTPss          = fRuleFit.GetRuleEnsemble().GetRulePSS(i);
            fNTPsb          = fRuleFit.GetRuleEnsemble().GetRulePSB(i);
            fNTPbs          = fRuleFit.GetRuleEnsemble().GetRulePBS(i);
            fNTPbb          = fRuleFit.GetRuleEnsemble().GetRulePBB(i);
            fNTSSB          = rule->GetSSB();
            fMonitorNtuple->Fill();
        }

        fRuleFit.MakeVisHists();
        fRuleFit.MakeDebugHists();
   }
   Log() << kDEBUG << "Training done" << Endl;

}

////////////////////////////////////////////////////////////////////////////////
/// training of rules using Jerome Friedmans implementation

void TMVA::MethodRuleFit::TrainJFRuleFit( void )
{
   fRuleFit.InitPtrs( this );
   Data()->SetCurrentType(Types::kTraining);
   UInt_t nevents = Data()->GetNTrainingEvents();
   std::vector<const TMVA::Event*> tmp;
   for (Long64_t ievt=0; ievt<nevents; ievt++) {
      const Event *event = GetEvent(ievt);
      tmp.push_back(event);
   }
   fRuleFit.SetTrainingEvents( tmp );

   RuleFitAPI *rfAPI = new RuleFitAPI( this, &fRuleFit, Log().GetMinType() );

   rfAPI->WelcomeMessage();

   // timer
   Timer timer( 1, GetName() );

   Log() << kINFO << "Training ..." << Endl;
   rfAPI->TrainRuleFit();

   Log() << kDEBUG << "reading model summary from rf_go.exe output" << Endl;
   rfAPI->ReadModelSum();

   //   fRuleFit.GetRuleEnsemblePtr()->MakeRuleMap();

   Log() << kDEBUG << "calculating rule and variable importance" << Endl;
   fRuleFit.CalcImportance();

   // Output results and fill monitor ntuple
   fRuleFit.GetRuleEnsemblePtr()->Print();
   //
   if(!IsSilentFile())fRuleFit.MakeVisHists();

   delete rfAPI;

   Log() << kDEBUG << "done training" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// computes ranking of input variables

const TMVA::Ranking* TMVA::MethodRuleFit::CreateRanking()
{
   // create the ranking object
   fRanking = new Ranking( GetName(), "Importance" );

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fRanking->AddRank( Rank( GetInputLabel(ivar), fRuleFit.GetRuleEnsemble().GetVarImportance(ivar) ) );
   }

   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////
/// add the rules to XML node

void TMVA::MethodRuleFit::AddWeightsXMLTo( void* parent ) const
{
   fRuleFit.GetRuleEnsemble().AddXMLTo( parent );
}

////////////////////////////////////////////////////////////////////////////////
/// read rules from an std::istream

void TMVA::MethodRuleFit::ReadWeightsFromStream( std::istream & istr )
{
   fRuleFit.GetRuleEnsemblePtr()->ReadRaw( istr );
}

////////////////////////////////////////////////////////////////////////////////
/// read rules from XML node

void TMVA::MethodRuleFit::ReadWeightsFromXML( void* wghtnode )
{
   fRuleFit.GetRuleEnsemblePtr()->ReadFromXML( wghtnode );
}

////////////////////////////////////////////////////////////////////////////////
/// returns MVA value for given event

Double_t TMVA::MethodRuleFit::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   // cannot determine error
   NoErrorCalc(err, errUpper);

   return fRuleFit.EvalEvent( *GetEvent() );
}

////////////////////////////////////////////////////////////////////////////////
/// write special monitoring histograms to file (here ntuple)

void  TMVA::MethodRuleFit::WriteMonitoringHistosToFile( void ) const
{
   BaseDir()->cd();
   Log() << kINFO << "Write monitoring ntuple to file: " << BaseDir()->GetPath() << Endl;
   fMonitorNtuple->Write();
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodRuleFit::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   Int_t dp = fout.precision();
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "};" << std::endl;
   fout << "void   " << className << "::Initialize(){}" << std::endl;
   fout << "void   " << className << "::Clear(){}" << std::endl;
   fout << "double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const {" << std::endl;
   fout << "   double rval=" << std::setprecision(10) << fRuleFit.GetRuleEnsemble().GetOffset() << ";" << std::endl;
   MakeClassRuleCuts(fout);
   MakeClassLinear(fout);
   fout << "   return rval;" << std::endl;
   fout << "}" << std::endl;
   fout << std::setprecision(dp);
}

////////////////////////////////////////////////////////////////////////////////
/// print out the rule cuts

void TMVA::MethodRuleFit::MakeClassRuleCuts( std::ostream& fout ) const
{
   Int_t dp = fout.precision();
   if (!fRuleFit.GetRuleEnsemble().DoRules()) {
      fout << "   //" << std::endl;
      fout << "   // ==> MODEL CONTAINS NO RULES <==" << std::endl;
      fout << "   //" << std::endl;
      return;
   }
   const RuleEnsemble *rens = &(fRuleFit.GetRuleEnsemble());
   const std::vector< Rule* > *rules = &(rens->GetRulesConst());
   const RuleCut *ruleCut;
   //
   std::list< std::pair<Double_t,Int_t> > sortedRules;
   for (UInt_t ir=0; ir<rules->size(); ir++) {
      sortedRules.push_back( std::pair<Double_t,Int_t>( (*rules)[ir]->GetImportance()/rens->GetImportanceRef(),ir ) );
   }
   sortedRules.sort();
   //
   fout << "   //" << std::endl;
   fout << "   // here follows all rules ordered in importance (most important first)" << std::endl;
   fout << "   // at the end of each line, the relative importance of the rule is given" << std::endl;
   fout << "   //" << std::endl;
   //
   for ( std::list< std::pair<double,int> >::reverse_iterator itpair = sortedRules.rbegin();
         itpair != sortedRules.rend(); ++itpair ) {
      UInt_t ir     = itpair->second;
      Double_t impr = itpair->first;
      ruleCut = (*rules)[ir]->GetRuleCut();
      if (impr<rens->GetImportanceCut()) fout << "   //" << std::endl;
      fout << "   if (" << std::flush;
      for (UInt_t ic=0; ic<ruleCut->GetNvars(); ic++) {
         Double_t sel    = ruleCut->GetSelector(ic);
         Double_t valmin = ruleCut->GetCutMin(ic);
         Double_t valmax = ruleCut->GetCutMax(ic);
         Bool_t   domin  = ruleCut->GetCutDoMin(ic);
         Bool_t   domax  = ruleCut->GetCutDoMax(ic);
         //
         if (ic>0) fout << "&&" << std::flush;
         if (domin) {
            fout << "(" << std::setprecision(10) << valmin << std::flush;
            fout << "<inputValues[" << sel << "])" << std::flush;
         }
         if (domax) {
            if (domin) fout << "&&" << std::flush;
            fout << "(inputValues[" << sel << "]" << std::flush;
            fout << "<" << std::setprecision(10) << valmax << ")" <<std::flush;
         }
      }
      fout << ") rval+=" << std::setprecision(10) << (*rules)[ir]->GetCoefficient() << ";" << std::flush;
      fout << "   // importance = " << Form("%3.3f",impr) << std::endl;
   }
   fout << std::setprecision(dp);
}

////////////////////////////////////////////////////////////////////////////////
/// print out the linear terms

void TMVA::MethodRuleFit::MakeClassLinear( std::ostream& fout ) const
{
   if (!fRuleFit.GetRuleEnsemble().DoLinear()) {
      fout << "   //" << std::endl;
      fout << "   // ==> MODEL CONTAINS NO LINEAR TERMS <==" << std::endl;
      fout << "   //" << std::endl;
      return;
   }
   fout << "   //" << std::endl;
   fout << "   // here follows all linear terms" << std::endl;
   fout << "   // at the end of each line, the relative importance of the term is given" << std::endl;
   fout << "   //" << std::endl;
   const RuleEnsemble *rens = &(fRuleFit.GetRuleEnsemble());
   UInt_t nlin = rens->GetNLinear();
   for (UInt_t il=0; il<nlin; il++) {
      if (rens->IsLinTermOK(il)) {
         Double_t norm = rens->GetLinNorm(il);
         Double_t imp  = rens->GetLinImportance(il)/rens->GetImportanceRef();
         fout << "   rval+="
            //           << std::setprecision(10) << rens->GetLinCoefficients(il)*norm << "*std::min(" << setprecision(10) << rens->GetLinDP(il)
            //           << ", std::max( inputValues[" << il << "]," << std::setprecision(10) << rens->GetLinDM(il) << "));"
              << std::setprecision(10) << rens->GetLinCoefficients(il)*norm
              << "*std::min( double(" << std::setprecision(10) << rens->GetLinDP(il)
              << "), std::max( double(inputValues[" << il << "]), double(" << std::setprecision(10) << rens->GetLinDM(il) << ")));"
              << std::flush;
         fout << "   // importance = " << Form("%3.3f",imp) << std::endl;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodRuleFit::GetHelpMessage() const
{
   TString col    = gConfig().WriteOptionsReference() ? TString() : gTools().Color("bold");
   TString colres = gConfig().WriteOptionsReference() ? TString() : gTools().Color("reset");
   TString brk    = gConfig().WriteOptionsReference() ? "<br>" : "";

   Log() << Endl;
   Log() << col << "--- Short description:" << colres << Endl;
   Log() << Endl;
   Log() << "This method uses a collection of so called rules to create a" << Endl;
   Log() << "discriminating scoring function. Each rule consists of a series" << Endl;
   Log() << "of cuts in parameter space. The ensemble of rules are created" << Endl;
   Log() << "from a forest of decision trees, trained using the training data." << Endl;
   Log() << "Each node (apart from the root) corresponds to one rule." << Endl;
   Log() << "The scoring function is then obtained by linearly combining" << Endl;
   Log() << "the rules. A fitting procedure is applied to find the optimum" << Endl;
   Log() << "set of coefficients. The goal is to find a model with few rules" << Endl;
   Log() << "but with a strong discriminating power." << Endl;
   Log() << Endl;
   Log() << col << "--- Performance optimisation:" << colres << Endl;
   Log() << Endl;
   Log() << "There are two important considerations to make when optimising:" << Endl;
   Log() << Endl;
   Log() << "  1. Topology of the decision tree forest" << brk << Endl;
   Log() << "  2. Fitting of the coefficients" << Endl;
   Log() << Endl;
   Log() << "The maximum complexity of the rules is defined by the size of" << Endl;
   Log() << "the trees. Large trees will yield many complex rules and capture" << Endl;
   Log() << "higher order correlations. On the other hand, small trees will" << Endl;
   Log() << "lead to a smaller ensemble with simple rules, only capable of" << Endl;
   Log() << "modeling simple structures." << Endl;
   Log() << "Several parameters exists for controlling the complexity of the" << Endl;
   Log() << "rule ensemble." << Endl;
   Log() << Endl;
   Log() << "The fitting procedure searches for a minimum using a gradient" << Endl;
   Log() << "directed path. Apart from step size and number of steps, the" << Endl;
   Log() << "evolution of the path is defined by a cut-off parameter, tau." << Endl;
   Log() << "This parameter is unknown and depends on the training data." << Endl;
   Log() << "A large value will tend to give large weights to a few rules." << Endl;
   Log() << "Similarly, a small value will lead to a large set of rules" << Endl;
   Log() << "with similar weights." << Endl;
   Log() << Endl;
   Log() << "A final point is the model used; rules and/or linear terms." << Endl;
   Log() << "For a given training sample, the result may improve by adding" << Endl;
   Log() << "linear terms. If best performance is obtained using only linear" << Endl;
   Log() << "terms, it is very likely that the Fisher discriminant would be" << Endl;
   Log() << "a better choice. Ideally the fitting procedure should be able to" << Endl;
   Log() << "make this choice by giving appropriate weights for either terms." << Endl;
   Log() << Endl;
   Log() << col << "--- Performance tuning via configuration options:" << colres << Endl;
   Log() << Endl;
   Log() << "I.  TUNING OF RULE ENSEMBLE:" << Endl;
   Log() << Endl;
   Log() << "   " << col << "ForestType  " << colres
         << ": Recommended is to use the default \"AdaBoost\"." << brk << Endl;
   Log() << "   " << col << "nTrees      " << colres
         << ": More trees leads to more rules but also slow" << Endl;
   Log() << "                 performance. With too few trees the risk is" << Endl;
   Log() << "                 that the rule ensemble becomes too simple." << brk << Endl;
   Log() << "   " << col << "fEventsMin  " << colres << brk << Endl;
   Log() << "   " << col << "fEventsMax  " << colres
         << ": With a lower min, more large trees will be generated" << Endl;
   Log() << "                 leading to more complex rules." << Endl;
   Log() << "                 With a higher max, more small trees will be" << Endl;
   Log() << "                 generated leading to more simple rules." << Endl;
   Log() << "                 By changing this range, the average complexity" << Endl;
   Log() << "                 of the rule ensemble can be controlled." << brk << Endl;
   Log() << "   " << col << "RuleMinDist " << colres
         << ": By increasing the minimum distance between" << Endl;
   Log() << "                 rules, fewer and more diverse rules will remain." << Endl;
   Log() << "                 Initially it is a good idea to keep this small" << Endl;
   Log() << "                 or zero and let the fitting do the selection of" << Endl;
   Log() << "                 rules. In order to reduce the ensemble size," << Endl;
   Log() << "                 the value can then be increased." << Endl;
   Log() << Endl;
   //         "|--------------------------------------------------------------|"
   Log() << "II. TUNING OF THE FITTING:" << Endl;
   Log() << Endl;
   Log() << "   " << col << "GDPathEveFrac " << colres
         << ": fraction of events in path evaluation" << Endl;
   Log() << "                 Increasing this fraction will improve the path" << Endl;
   Log() << "                 finding. However, a too high value will give few" << Endl;
   Log() << "                 unique events available for error estimation." << Endl;
   Log() << "                 It is recommended to use the default = 0.5." << brk << Endl;
   Log() << "   " << col << "GDTau         " << colres
         << ": cutoff parameter tau" << Endl;
   Log() << "                 By default this value is set to -1.0." << Endl;
   //         "|----------------|---------------------------------------------|"
   Log() << "                 This means that the cut off parameter is" << Endl;
   Log() << "                 automatically estimated. In most cases" << Endl;
   Log() << "                 this should be fine. However, you may want" << Endl;
   Log() << "                 to fix this value if you already know it" << Endl;
   Log() << "                 and want to reduce on training time." << brk << Endl;
   Log() << "   " << col << "GDTauPrec     " << colres
         << ": precision of estimated tau" << Endl;
   Log() << "                 Increase this precision to find a more" << Endl;
   Log() << "                 optimum cut-off parameter." << brk << Endl;
   Log() << "   " << col << "GDNStep       " << colres
         << ": number of steps in path search" << Endl;
   Log() << "                 If the number of steps is too small, then" << Endl;
   Log() << "                 the program will give a warning message." << Endl;
   Log() << Endl;
   Log() << "III. WARNING MESSAGES" << Endl;
   Log() << Endl;
   Log() << col << "Risk(i+1)>=Risk(i) in path" << colres << brk << Endl;
   Log() << col << "Chaotic behaviour of risk evolution." << colres << Endl;
   //         "|----------------|---------------------------------------------|"
   Log() << "                 The error rate was still decreasing at the end" << Endl;
   Log() << "                 By construction the Risk should always decrease." << Endl;
   Log() << "                 However, if the training sample is too small or" << Endl;
   Log() << "                 the model is overtrained, such warnings can" << Endl;
   Log() << "                 occur." << Endl;
   Log() << "                 The warnings can safely be ignored if only a" << Endl;
   Log() << "                 few (<3) occur. If more warnings are generated," << Endl;
   Log() << "                 the fitting fails." << Endl;
   Log() << "                 A remedy may be to increase the value" << brk << Endl;
   Log() << "                 "
         << col << "GDValidEveFrac" << colres
         << " to 1.0 (or a larger value)." << brk << Endl;
   Log() << "                 In addition, if "
         << col << "GDPathEveFrac" << colres
         << " is too high" << Endl;
   Log() << "                 the same warnings may occur since the events" << Endl;
   Log() << "                 used for error estimation are also used for" << Endl;
   Log() << "                 path estimation." << Endl;
   Log() << "                 Another possibility is to modify the model - " << Endl;
   Log() << "                 See above on tuning the rule ensemble." << Endl;
   Log() << Endl;
   Log() << col << "The error rate was still decreasing at the end of the path"
         << colres << Endl;
   Log() << "                 Too few steps in path! Increase "
         << col << "GDNSteps" <<  colres << "." << Endl;
   Log() << Endl;
   Log() << col << "Reached minimum early in the search" << colres << Endl;

   Log() << "                 Minimum was found early in the fitting. This" << Endl;
   Log() << "                 may indicate that the used step size "
         << col << "GDStep" <<  colres << "." << Endl;
   Log() << "                 was too large. Reduce it and rerun." << Endl;
   Log() << "                 If the results still are not OK, modify the" << Endl;
   Log() << "                 model either by modifying the rule ensemble" << Endl;
   Log() << "                 or add/remove linear terms" << Endl;
}
