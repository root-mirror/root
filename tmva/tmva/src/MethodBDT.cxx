// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard v. Toerne, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodBDT (BDT = Boosted Decision Trees)                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Analysis of Boosted Decision Trees                                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Doug Schouten   <dschoute@sfu.ca>        - Simon Fraser U., Canada        *
 *      Jan Therhaag    <jan.therhaag@cern.ch>   - U. of Bonn, Germany            *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
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

/*! \class TMVA::MethodBDT
\ingroup TMVA

Analysis of Boosted Decision Trees

Boosted decision trees have been successfully used in High Energy
Physics analysis for example by the MiniBooNE experiment
(Yang-Roe-Zhu, physics/0508045). In Boosted Decision Trees, the
selection is done on a majority vote on the result of several decision
trees, which are all derived from the same training sample by
supplying different event weights during the training.

### Decision trees:

Successive decision nodes are used to categorize the
events out of the sample as either signal or background. Each node
uses only a single discriminating variable to decide if the event is
signal-like ("goes right") or background-like ("goes left"). This
forms a tree like structure with "baskets" at the end (leave nodes),
and an event is classified as either signal or background according to
whether the basket where it ends up has been classified signal or
background during the training. Training of a decision tree is the
process to define the "cut criteria" for each node. The training
starts with the root node. Here one takes the full training event
sample and selects the variable and corresponding cut value that gives
the best separation between signal and background at this stage. Using
this cut criterion, the sample is then divided into two subsamples, a
signal-like (right) and a background-like (left) sample. Two new nodes
are then created for each of the two sub-samples and they are
constructed using the same mechanism as described for the root
node. The devision is stopped once a certain node has reached either a
minimum number of events, or a minimum or maximum signal purity. These
leave nodes are then called "signal" or "background" if they contain
more signal respective background events from the training sample.

### Boosting:

The idea behind adaptive boosting (AdaBoost) is, that signal events
from the training sample, that end up in a background node
(and vice versa) are given a larger weight than events that are in
the correct leave node. This results in a re-weighed training event
sample, with which then a new decision tree can be developed.
The boosting can be applied several times (typically 100-500 times)
and one ends up with a set of decision trees (a forest).
Gradient boosting works more like a function expansion approach, where
each tree corresponds to a summand. The parameters for each summand (tree)
are determined by the minimization of a error function (binomial log-
likelihood for classification and Huber loss for regression).
A greedy algorithm is used, which means, that only one tree is modified
at a time, while the other trees stay fixed.

### Bagging:

In this particular variant of the Boosted Decision Trees the boosting
is not done on the basis of previous training results, but by a simple
stochastic re-sampling of the initial training event sample.

### Random Trees:

Similar to the "Random Forests" from Leo Breiman and Adele Cutler, it
uses the bagging algorithm together and bases the determination of the
best node-split during the training on a random subset of variables only
which is individually chosen for each split.

### Analysis:

Applying an individual decision tree to a test event results in a
classification of the event as either signal or background. For the
boosted decision tree selection, an event is successively subjected to
the whole set of decision trees and depending on how often it is
classified as signal, a "likelihood" estimator is constructed for the
event being signal or background. The value of this estimator is the
one which is then used to select the events from an event sample, and
the cut value on this estimator defines the efficiency and purity of
the selection.

*/


#include "TMVA/MethodBDT.h"
#include "TMVA/Config.h"

#include "TMVA/BDTEventWrapper.h"
#include "TMVA/BinarySearchTree.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/CrossEntropy.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/DataSet.h"
#include "TMVA/GiniIndex.h"
#include "TMVA/GiniIndexWithLaplace.h"
#include "TMVA/Interval.h"
#include "TMVA/IMethod.h"
#include "TMVA/LogInterval.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MisClassificationError.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/OptimizeConfigParameters.h"
#include "TMVA/PDF.h"
#include "TMVA/Ranking.h"
#include "TMVA/Results.h"
#include "TMVA/ResultsMulticlass.h"
#include "TMVA/SdivSqrtSplusB.h"
#include "TMVA/SeparationBase.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TRandom3.h"
#include "TMath.h"
#include "TMatrixTSym.h"
#include "TGraph.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

using std::vector;
using std::make_pair;

REGISTER_METHOD(BDT)

ClassImp(TMVA::MethodBDT);

   const Int_t TMVA::MethodBDT::fgDebugLevel = 0;

////////////////////////////////////////////////////////////////////////////////
/// The standard constructor for the "boosted decision trees".

TMVA::MethodBDT::MethodBDT( const TString& jobName,
                            const TString& methodTitle,
                            DataSetInfo& theData,
                            const TString& theOption ) :
   TMVA::MethodBase( jobName, Types::kBDT, methodTitle, theData, theOption)
   , fTrainSample(0)
   , fNTrees(0)
   , fSigToBkgFraction(0)
   , fAdaBoostBeta(0)
//   , fTransitionPoint(0)
   , fShrinkage(0)
   , fBaggedBoost(kFALSE)
   , fBaggedGradBoost(kFALSE)
//   , fSumOfWeights(0)
   , fMinNodeEvents(0)
   , fMinNodeSize(5)
   , fMinNodeSizeS("5%")
   , fNCuts(0)
   , fUseFisherCuts(0)        // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fMinLinCorrForFisher(.8) // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fUseExclusiveVars(0)     // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fUseYesNoLeaf(kFALSE)
   , fNodePurityLimit(0)
   , fNNodesMax(0)
   , fMaxDepth(0)
   , fPruneMethod(DecisionTree::kNoPruning)
   , fPruneStrength(0)
   , fFValidationEvents(0)
   , fAutomatic(kFALSE)
   , fRandomisedTrees(kFALSE)
   , fUseNvars(0)
   , fUsePoissonNvars(0)  // don't use this initialisation, only here to make  Coverity happy. Is set in Init()
   , fUseNTrainEvents(0)
   , fBaggedSampleFraction(0)
   , fNoNegWeightsInTraining(kFALSE)
   , fInverseBoostNegWeights(kFALSE)
   , fPairNegWeightsGlobal(kFALSE)
   , fTrainWithNegWeights(kFALSE)
   , fDoBoostMonitor(kFALSE)
   , fITree(0)
   , fBoostWeight(0)
   , fErrorFraction(0)
   , fCss(0)
   , fCts_sb(0)
   , fCtb_ss(0)
   , fCbb(0)
   , fDoPreselection(kFALSE)
   , fSkipNormalization(kFALSE)
   , fHistoricBool(kFALSE)
{
   fMonitorNtuple = NULL;
   fSepType = NULL;
   fRegressionLossFunctionBDTG = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::MethodBDT::MethodBDT( DataSetInfo& theData,
                            const TString& theWeightFile)
   : TMVA::MethodBase( Types::kBDT, theData, theWeightFile)
   , fTrainSample(0)
   , fNTrees(0)
   , fSigToBkgFraction(0)
   , fAdaBoostBeta(0)
//   , fTransitionPoint(0)
   , fShrinkage(0)
   , fBaggedBoost(kFALSE)
   , fBaggedGradBoost(kFALSE)
//   , fSumOfWeights(0)
   , fMinNodeEvents(0)
   , fMinNodeSize(5)
   , fMinNodeSizeS("5%")
   , fNCuts(0)
   , fUseFisherCuts(0)        // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fMinLinCorrForFisher(.8) // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fUseExclusiveVars(0)     // don't use this initialisation, only here to make  Coverity happy. Is set in DeclarOptions()
   , fUseYesNoLeaf(kFALSE)
   , fNodePurityLimit(0)
   , fNNodesMax(0)
   , fMaxDepth(0)
   , fPruneMethod(DecisionTree::kNoPruning)
   , fPruneStrength(0)
   , fFValidationEvents(0)
   , fAutomatic(kFALSE)
   , fRandomisedTrees(kFALSE)
   , fUseNvars(0)
   , fUsePoissonNvars(0)  // don't use this initialisation, only here to make  Coverity happy. Is set in Init()
   , fUseNTrainEvents(0)
   , fBaggedSampleFraction(0)
   , fNoNegWeightsInTraining(kFALSE)
   , fInverseBoostNegWeights(kFALSE)
   , fPairNegWeightsGlobal(kFALSE)
   , fTrainWithNegWeights(kFALSE)
   , fDoBoostMonitor(kFALSE)
   , fITree(0)
   , fBoostWeight(0)
   , fErrorFraction(0)
   , fCss(0)
   , fCts_sb(0)
   , fCtb_ss(0)
   , fCbb(0)
   , fDoPreselection(kFALSE)
   , fSkipNormalization(kFALSE)
   , fHistoricBool(kFALSE)
{
   fMonitorNtuple = NULL;
   fSepType = NULL;
   fRegressionLossFunctionBDTG = nullptr;
   // constructor for calculating BDT-MVA using previously generated decision trees
   // the result of the previous training (the decision trees) are read in via the
   // weight file. Make sure the the variables correspond to the ones used in
   // creating the "weight"-file
}

////////////////////////////////////////////////////////////////////////////////
/// BDT can handle classification with multiple classes and regression with one regression-target.

Bool_t TMVA::MethodBDT::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass ) return kTRUE;
   if( type == Types::kRegression && numberTargets == 1 ) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Define the options (their key words). That can be set in the option string.
///
/// know options:
///
///  - nTrees        number of trees in the forest to be created
///  - BoostType     the boosting type for the trees in the forest (AdaBoost e.t.c..).
///                  Known:
///                        - AdaBoost
///                        - AdaBoostR2 (Adaboost for regression)
///                        - Bagging
///                        - GradBoost
///  - AdaBoostBeta     the boosting parameter, beta, for AdaBoost
///  - UseRandomisedTrees  choose at each node splitting a random set of variables
///  - UseNvars         use UseNvars variables in randomised trees
///  - UsePoisson Nvars use UseNvars not as fixed number but as mean of a poisson distribution
///  - SeparationType   the separation criterion applied in the node splitting.
///                  Known:
///                        - GiniIndex
///                        - MisClassificationError
///                        - CrossEntropy
///                        - SDivSqrtSPlusB
///  - MinNodeSize:     minimum percentage of training events in a leaf node (leaf criteria, stop splitting)
///  - nCuts:           the number of steps in the optimisation of the cut for a node (if < 0, then
///                  step size is determined by the events)
///  - UseFisherCuts:   use multivariate splits using the Fisher criterion
///  - UseYesNoLeaf     decide if the classification is done simply by the node type, or the S/B
///                  (from the training) in the leaf node
///  - NodePurityLimit  the minimum purity to classify a node as a signal node (used in pruning and boosting to determine
///                  misclassification error rate)
///  - PruneMethod      The Pruning method.
///                  Known:
///                        - NoPruning  // switch off pruning completely
///                        - ExpectedError
///                        - CostComplexity
///  - PruneStrength    a parameter to adjust the amount of pruning. Should be large enough such that overtraining is avoided.
///  - PruningValFraction   number of events to use for optimizing pruning (only if PruneStrength < 0, i.e. automatic pruning)
///  - NegWeightTreatment
///                       - IgnoreNegWeightsInTraining  Ignore negative weight events in the training.
///                       -  DecreaseBoostWeight     Boost ev. with neg. weight with 1/boostweight instead of boostweight
///                       -  PairNegWeightsGlobal    Pair ev. with neg. and pos. weights in training sample and "annihilate" them
///  - MaxDepth         maximum depth of the decision tree allowed before further splitting is stopped
///  - SkipNormalization       Skip normalization at initialization, to keep expectation value of BDT output
///             according to the fraction of events

void TMVA::MethodBDT::DeclareOptions()
{
   DeclareOptionRef(fNTrees, "NTrees", "Number of trees in the forest");
   if (DoRegression()) {
      DeclareOptionRef(fMaxDepth=50,"MaxDepth","Max depth of the decision tree allowed");
   }else{
      DeclareOptionRef(fMaxDepth=3,"MaxDepth","Max depth of the decision tree allowed");
   }

   TString tmp="5%"; if (DoRegression()) tmp="0.2%";
   DeclareOptionRef(fMinNodeSizeS=tmp, "MinNodeSize", "Minimum percentage of training events required in a leaf node (default: Classification: 5%, Regression: 0.2%)");
   // MinNodeSize:     minimum percentage of training events in a leaf node (leaf criteria, stop splitting)
   DeclareOptionRef(fNCuts, "nCuts", "Number of grid points in variable range used in finding optimal cut in node splitting");

   DeclareOptionRef(fBoostType, "BoostType", "Boosting type for the trees in the forest (note: AdaCost is still experimental)");

   AddPreDefVal(TString("AdaBoost"));
   AddPreDefVal(TString("RealAdaBoost"));
   AddPreDefVal(TString("AdaCost"));
   AddPreDefVal(TString("Bagging"));
   //   AddPreDefVal(TString("RegBoost"));
   AddPreDefVal(TString("AdaBoostR2"));
   AddPreDefVal(TString("Grad"));
   if (DoRegression()) {
      fBoostType = "AdaBoostR2";
   }else{
      fBoostType = "AdaBoost";
   }
   DeclareOptionRef(fAdaBoostR2Loss="Quadratic", "AdaBoostR2Loss", "Type of Loss function in AdaBoostR2");
   AddPreDefVal(TString("Linear"));
   AddPreDefVal(TString("Quadratic"));
   AddPreDefVal(TString("Exponential"));

   DeclareOptionRef(fBaggedBoost=kFALSE, "UseBaggedBoost","Use only a random subsample of all events for growing the trees in each boost iteration.");
   DeclareOptionRef(fShrinkage = 1.0, "Shrinkage", "Learning rate for BoostType=Grad algorithm");
   DeclareOptionRef(fAdaBoostBeta=.5, "AdaBoostBeta", "Learning rate  for AdaBoost algorithm");
   DeclareOptionRef(fRandomisedTrees,"UseRandomisedTrees","Determine at each node splitting the cut variable only as the best out of a random subset of variables (like in RandomForests)");
   DeclareOptionRef(fUseNvars,"UseNvars","Size of the subset of variables used with RandomisedTree option");
   DeclareOptionRef(fUsePoissonNvars,"UsePoissonNvars", "Interpret \"UseNvars\" not as fixed number but as mean of a Poisson distribution in each split with RandomisedTree option");
   DeclareOptionRef(fBaggedSampleFraction=.6,"BaggedSampleFraction","Relative size of bagged event sample to original size of the data sample (used whenever bagging is used (i.e. UseBaggedBoost, Bagging,)" );

   DeclareOptionRef(fUseYesNoLeaf=kTRUE, "UseYesNoLeaf",
                    "Use Sig or Bkg categories, or the purity=S/(S+B) as classification of the leaf node -> Real-AdaBoost");
   if (DoRegression()) {
      fUseYesNoLeaf = kFALSE;
   }

   DeclareOptionRef(fNegWeightTreatment="InverseBoostNegWeights","NegWeightTreatment","How to treat events with negative weights in the BDT training (particular the boosting) : IgnoreInTraining;  Boost With inverse boostweight; Pair events with negative and positive weights in training sample and *annihilate* them (experimental!)");
   AddPreDefVal(TString("InverseBoostNegWeights"));
   AddPreDefVal(TString("IgnoreNegWeightsInTraining"));
   AddPreDefVal(TString("NoNegWeightsInTraining"));    // well, let's be nice to users and keep at least this old name anyway ..
   AddPreDefVal(TString("PairNegWeightsGlobal"));
   AddPreDefVal(TString("Pray"));



   DeclareOptionRef(fCss=1.,   "Css",   "AdaCost: cost of true signal selected signal");
   DeclareOptionRef(fCts_sb=1.,"Cts_sb","AdaCost: cost of true signal selected bkg");
   DeclareOptionRef(fCtb_ss=1.,"Ctb_ss","AdaCost: cost of true bkg    selected signal");
   DeclareOptionRef(fCbb=1.,   "Cbb",   "AdaCost: cost of true bkg    selected bkg ");

   DeclareOptionRef(fNodePurityLimit=0.5, "NodePurityLimit", "In boosting/pruning, nodes with purity > NodePurityLimit are signal; background otherwise.");


   DeclareOptionRef(fSepTypeS, "SeparationType", "Separation criterion for node splitting");
   AddPreDefVal(TString("CrossEntropy"));
   AddPreDefVal(TString("GiniIndex"));
   AddPreDefVal(TString("GiniIndexWithLaplace"));
   AddPreDefVal(TString("MisClassificationError"));
   AddPreDefVal(TString("SDivSqrtSPlusB"));
   AddPreDefVal(TString("RegressionVariance"));
   if (DoRegression()) {
      fSepTypeS = "RegressionVariance";
   }else{
      fSepTypeS = "GiniIndex";
   }

   DeclareOptionRef(fRegressionLossFunctionBDTGS = "Huber", "RegressionLossFunctionBDTG", "Loss function for BDTG regression.");
   AddPreDefVal(TString("Huber"));
   AddPreDefVal(TString("AbsoluteDeviation"));
   AddPreDefVal(TString("LeastSquares"));

   DeclareOptionRef(fHuberQuantile = 0.7, "HuberQuantile", "In the Huber loss function this is the quantile that separates the core from the tails in the residuals distribution.");

   DeclareOptionRef(fDoBoostMonitor=kFALSE,"DoBoostMonitor","Create control plot with ROC integral vs tree number");

   DeclareOptionRef(fUseFisherCuts=kFALSE, "UseFisherCuts", "Use multivariate splits using the Fisher criterion");
   DeclareOptionRef(fMinLinCorrForFisher=.8,"MinLinCorrForFisher", "The minimum linear correlation between two variables demanded for use in Fisher criterion in node splitting");
   DeclareOptionRef(fUseExclusiveVars=kFALSE,"UseExclusiveVars","Variables already used in fisher criterion are not anymore analysed individually for node splitting");


   DeclareOptionRef(fDoPreselection=kFALSE,"DoPreselection","and and apply automatic pre-selection for 100% efficient signal (bkg) cuts prior to training");


   DeclareOptionRef(fSigToBkgFraction=1,"SigToBkgFraction","Sig to Bkg ratio used in Training (similar to NodePurityLimit, which cannot be used in real adaboost");

   DeclareOptionRef(fPruneMethodS, "PruneMethod", "Note: for BDTs use small trees (e.g.MaxDepth=3) and NoPruning:  Pruning: Method used for pruning (removal) of statistically insignificant branches ");
   AddPreDefVal(TString("NoPruning"));
   AddPreDefVal(TString("ExpectedError"));
   AddPreDefVal(TString("CostComplexity"));

   DeclareOptionRef(fPruneStrength, "PruneStrength", "Pruning strength");

   DeclareOptionRef(fFValidationEvents=0.5, "PruningValFraction", "Fraction of events to use for optimizing automatic pruning.");

   DeclareOptionRef(fSkipNormalization=kFALSE, "SkipNormalization", "Skip normalization at initialization, to keep expectation value of BDT output according to the fraction of events");

    // deprecated options, still kept for the moment:
   DeclareOptionRef(fMinNodeEvents=0, "nEventsMin", "deprecated: Use MinNodeSize (in % of training events) instead");

   DeclareOptionRef(fBaggedGradBoost=kFALSE, "UseBaggedGrad","deprecated: Use *UseBaggedBoost* instead:  Use only a random subsample of all events for growing the trees in each iteration.");
   DeclareOptionRef(fBaggedSampleFraction, "GradBaggingFraction","deprecated: Use *BaggedSampleFraction* instead: Defines the fraction of events to be used in each iteration, e.g. when UseBaggedGrad=kTRUE. ");
   DeclareOptionRef(fUseNTrainEvents,"UseNTrainEvents","deprecated: Use *BaggedSampleFraction* instead: Number of randomly picked training events used in randomised (and bagged) trees");
   DeclareOptionRef(fNNodesMax,"NNodesMax","deprecated: Use MaxDepth instead to limit the tree size" );


}

////////////////////////////////////////////////////////////////////////////////
/// Options that are used ONLY for the READER to ensure backward compatibility.

void TMVA::MethodBDT::DeclareCompatibilityOptions() {
   MethodBase::DeclareCompatibilityOptions();


   DeclareOptionRef(fHistoricBool=kTRUE, "UseWeightedTrees",
                    "Use weighted trees or simple average in classification from the forest");
   DeclareOptionRef(fHistoricBool=kFALSE, "PruneBeforeBoost", "Flag to prune the tree before applying boosting algorithm");
   DeclareOptionRef(fHistoricBool=kFALSE,"RenormByClass","Individually re-normalize each event class to the original size after boosting");

   AddPreDefVal(TString("NegWeightTreatment"),TString("IgnoreNegWeights"));

}

////////////////////////////////////////////////////////////////////////////////
/// The option string is decoded, for available options see "DeclareOptions".

void TMVA::MethodBDT::ProcessOptions()
{
   fSepTypeS.ToLower();
   if      (fSepTypeS == "misclassificationerror") fSepType = new MisClassificationError();
   else if (fSepTypeS == "giniindex")              fSepType = new GiniIndex();
   else if (fSepTypeS == "giniindexwithlaplace")   fSepType = new GiniIndexWithLaplace();
   else if (fSepTypeS == "crossentropy")           fSepType = new CrossEntropy();
   else if (fSepTypeS == "sdivsqrtsplusb")         fSepType = new SdivSqrtSplusB();
   else if (fSepTypeS == "regressionvariance")     fSepType = NULL;
   else {
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<ProcessOptions> unknown Separation Index option " << fSepTypeS << " called" << Endl;
   }

   if(!(fHuberQuantile >= 0.0 && fHuberQuantile <= 1.0)){
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<ProcessOptions> Huber Quantile must be in range [0,1]. Value given, " << fHuberQuantile << ", does not match this criteria" << Endl;
   }


   fRegressionLossFunctionBDTGS.ToLower();
   if      (fRegressionLossFunctionBDTGS == "huber")                  fRegressionLossFunctionBDTG = new HuberLossFunctionBDT(fHuberQuantile);
   else if (fRegressionLossFunctionBDTGS == "leastsquares")           fRegressionLossFunctionBDTG = new LeastSquaresLossFunctionBDT();
   else if (fRegressionLossFunctionBDTGS == "absolutedeviation")      fRegressionLossFunctionBDTG = new AbsoluteDeviationLossFunctionBDT();
   else {
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<ProcessOptions> unknown Regression Loss Function BDT option " << fRegressionLossFunctionBDTGS << " called" << Endl;
   }

   fPruneMethodS.ToLower();
   if      (fPruneMethodS == "expectederror")  fPruneMethod = DecisionTree::kExpectedErrorPruning;
   else if (fPruneMethodS == "costcomplexity") fPruneMethod = DecisionTree::kCostComplexityPruning;
   else if (fPruneMethodS == "nopruning")      fPruneMethod = DecisionTree::kNoPruning;
   else {
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<ProcessOptions> unknown PruneMethod " << fPruneMethodS << " option called" << Endl;
   }
   if (fPruneStrength < 0 && (fPruneMethod != DecisionTree::kNoPruning) && fBoostType!="Grad") fAutomatic = kTRUE;
   else fAutomatic = kFALSE;
   if (fAutomatic && fPruneMethod==DecisionTree::kExpectedErrorPruning){
      Log() << kFATAL
            <<  "Sorry automatic pruning strength determination is not implemented yet for ExpectedErrorPruning" << Endl;
   }


   if (fMinNodeEvents > 0){
      fMinNodeSize = Double_t(fMinNodeEvents*100.) / Data()->GetNTrainingEvents();
      Log() << kWARNING << "You have explicitly set ** nEventsMin = " << fMinNodeEvents<<" ** the min absolute number \n"
            << "of events in a leaf node. This is DEPRECATED, please use the option \n"
            << "*MinNodeSize* giving the relative number as percentage of training \n"
            << "events instead. \n"
            << "nEventsMin="<<fMinNodeEvents<< "--> MinNodeSize="<<fMinNodeSize<<"%"
            << Endl;
      Log() << kWARNING << "Note also that explicitly setting *nEventsMin* so far OVERWRITES the option recommended \n"
            << " *MinNodeSize* = " << fMinNodeSizeS << " option !!" << Endl ;
      fMinNodeSizeS = Form("%F3.2",fMinNodeSize);

   }else{
      SetMinNodeSize(fMinNodeSizeS);
   }


   fAdaBoostR2Loss.ToLower();

   if (fBoostType=="Grad") {
      fPruneMethod = DecisionTree::kNoPruning;
      if (fNegWeightTreatment=="InverseBoostNegWeights"){
         Log() << kINFO << "the option NegWeightTreatment=InverseBoostNegWeights does"
               << " not exist for BoostType=Grad" << Endl;
         Log() << kINFO << "--> change to new default NegWeightTreatment=Pray" << Endl;
         Log() << kDEBUG << "i.e. simply keep them as if which should work fine for Grad Boost" << Endl;
         fNegWeightTreatment="Pray";
         fNoNegWeightsInTraining=kFALSE;
      }
   } else if (fBoostType=="RealAdaBoost"){
      fBoostType    = "AdaBoost";
      fUseYesNoLeaf = kFALSE;
   } else if (fBoostType=="AdaCost"){
      fUseYesNoLeaf = kFALSE;
   }

   if (fFValidationEvents < 0.0) fFValidationEvents = 0.0;
   if (fAutomatic && fFValidationEvents > 0.5) {
      Log() << kWARNING << "You have chosen to use more than half of your training sample "
            << "to optimize the automatic pruning algorithm. This is probably wasteful "
            << "and your overall results will be degraded. Are you sure you want this?"
            << Endl;
   }


   if (this->Data()->HasNegativeEventWeights()){
      Log() << kINFO << " You are using a Monte Carlo that has also negative weights. "
            << "That should in principle be fine as long as on average you end up with "
            << "something positive. For this you have to make sure that the minimal number "
            << "of (un-weighted) events demanded for a tree node (currently you use: MinNodeSize="
            << fMinNodeSizeS << "  ("<< fMinNodeSize << "%)"
            <<", (or the deprecated equivalent nEventsMin) you can set this via the "
            <<"BDT option string when booking the "
            << "classifier) is large enough to allow for reasonable averaging!!! "
            << " If this does not help.. maybe you want to try the option: IgnoreNegWeightsInTraining  "
            << "which ignores events with negative weight in the training. " << Endl
            << Endl << "Note: You'll get a WARNING message during the training if that should ever happen" << Endl;
   }

   if (DoRegression()) {
      if (fUseYesNoLeaf && !IsConstructedFromWeightFile()){
         Log() << kWARNING << "Regression Trees do not work with fUseYesNoLeaf=TRUE --> I will set it to FALSE" << Endl;
         fUseYesNoLeaf = kFALSE;
      }

      if (fSepType != NULL){
         Log() << kWARNING << "Regression Trees do not work with Separation type other than <RegressionVariance> --> I will use it instead" << Endl;
         fSepType = NULL;
      }
      if (fUseFisherCuts){
         Log() << kWARNING << "Sorry, UseFisherCuts is not available for regression analysis, I will ignore it!" << Endl;
         fUseFisherCuts = kFALSE;
      }
      if (fNCuts < 0) {
         Log() << kWARNING << "Sorry, the option of nCuts<0 using a more elaborate node splitting algorithm " << Endl;
         Log() << kWARNING << "is not implemented for regression analysis ! " << Endl;
         Log() << kWARNING << "--> I switch do default nCuts = 20 and use standard node splitting"<<Endl;
         fNCuts=20;
      }
   }
   if (fRandomisedTrees){
      Log() << kINFO << " Randomised trees use no pruning" << Endl;
      fPruneMethod = DecisionTree::kNoPruning;
      //      fBoostType   = "Bagging";
   }

   if (fUseFisherCuts) {
      Log() << kWARNING << "When using the option UseFisherCuts, the other option nCuts<0 (i.e. using" << Endl;
      Log() << " a more elaborate node splitting algorithm) is not implemented. " << Endl;
      //I will switch o " << Endl;
      //Log() << "--> I switch do default nCuts = 20 and use standard node splitting WITH possible Fisher criteria"<<Endl;
      fNCuts=20;
   }

   if (fNTrees==0){
      Log() << kERROR << " Zero Decision Trees demanded... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      fNTrees = 1;
   }

   fNegWeightTreatment.ToLower();
   if      (fNegWeightTreatment == "ignorenegweightsintraining")   fNoNegWeightsInTraining = kTRUE;
   else if (fNegWeightTreatment == "nonegweightsintraining")   fNoNegWeightsInTraining = kTRUE;
   else if (fNegWeightTreatment == "inverseboostnegweights") fInverseBoostNegWeights = kTRUE;
   else if (fNegWeightTreatment == "pairnegweightsglobal")   fPairNegWeightsGlobal   = kTRUE;
   else if (fNegWeightTreatment == "pray")   Log() << kDEBUG << "Yes, good luck with praying " << Endl;
   else {
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<ProcessOptions> unknown option for treating negative event weights during training " << fNegWeightTreatment << " requested" << Endl;
   }

   if (fNegWeightTreatment == "pairnegweightsglobal")
      Log() << kWARNING << " you specified the option NegWeightTreatment=PairNegWeightsGlobal : This option is still considered EXPERIMENTAL !! " << Endl;


   // dealing with deprecated options !
   if (fNNodesMax>0) {
      UInt_t tmp=1; // depth=0  == 1 node
      fMaxDepth=0;
      while (tmp < fNNodesMax){
         tmp+=2*tmp;
         fMaxDepth++;
      }
      Log() << kWARNING << "You have specified a deprecated option *NNodesMax="<<fNNodesMax
            << "* \n this has been translated to MaxDepth="<<fMaxDepth<<Endl;
   }


   if (fUseNTrainEvents>0){
      fBaggedSampleFraction  = (Double_t) fUseNTrainEvents/Data()->GetNTrainingEvents();
      Log() << kWARNING << "You have specified a deprecated option *UseNTrainEvents="<<fUseNTrainEvents
            << "* \n this has been translated to BaggedSampleFraction="<<fBaggedSampleFraction<<"(%)"<<Endl;
   }

   if (fBoostType=="Bagging") fBaggedBoost = kTRUE;
   if (fBaggedGradBoost){
      fBaggedBoost = kTRUE;
      Log() << kWARNING << "You have specified a deprecated option *UseBaggedGrad* --> please use  *UseBaggedBoost* instead" << Endl;
   }

}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBDT::SetMinNodeSize(Double_t sizeInPercent){
   if (sizeInPercent > 0 && sizeInPercent < 50){
      fMinNodeSize=sizeInPercent;

   } else {
      Log() << kFATAL << "you have demanded a minimal node size of "
            << sizeInPercent << "% of the training events.. \n"
            << " that somehow does not make sense "<<Endl;
   }

}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBDT::SetMinNodeSize(TString sizeInPercent){
   sizeInPercent.ReplaceAll("%","");
   sizeInPercent.ReplaceAll(" ","");
   if (sizeInPercent.IsFloat()) SetMinNodeSize(sizeInPercent.Atof());
   else {
      Log() << kFATAL << "I had problems reading the option MinNodeEvents, which "
            << "after removing a possible % sign now reads " << sizeInPercent << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialisation with defaults for the BDT-Method.

void TMVA::MethodBDT::Init( void )
{
   fNTrees         = 800;
   if (fAnalysisType == Types::kClassification || fAnalysisType == Types::kMulticlass ) {
      fMaxDepth        = 3;
      fBoostType      = "AdaBoost";
      if(DataInfo().GetNClasses()!=0) //workaround for multiclass application
         fMinNodeSize = 5.;
   }else {
      fMaxDepth = 50;
      fBoostType      = "AdaBoostR2";
      fAdaBoostR2Loss = "Quadratic";
      if(DataInfo().GetNClasses()!=0) //workaround for multiclass application
         fMinNodeSize  = .2;
   }


   fNCuts          = 20;
   fPruneMethodS   = "NoPruning";
   fPruneMethod    = DecisionTree::kNoPruning;
   fPruneStrength  = 0;
   fAutomatic      = kFALSE;
   fFValidationEvents = 0.5;
   fRandomisedTrees = kFALSE;
   //   fUseNvars        =  (GetNvar()>12) ? UInt_t(GetNvar()/8) : TMath::Max(UInt_t(2),UInt_t(GetNvar()/3));
   fUseNvars        =  UInt_t(TMath::Sqrt(GetNvar())+0.6);
   fUsePoissonNvars = kTRUE;
   fShrinkage       = 1.0;
//   fSumOfWeights    = 0.0;

   // reference cut value to distinguish signal-like from background-like events
   SetSignalReferenceCut( 0 );
}


////////////////////////////////////////////////////////////////////////////////
/// Reset the method, as if it had just been instantiated (forget all training etc.).

void TMVA::MethodBDT::Reset( void )
{
   // I keep the BDT EventSample and its Validation sample (eventually they should all
   // disappear and just use the DataSet samples ..

   // remove all the trees
   for (UInt_t i=0; i<fForest.size();           i++) delete fForest[i];
   fForest.clear();

   fBoostWeights.clear();
   if (fMonitorNtuple) { fMonitorNtuple->Delete(); fMonitorNtuple=NULL; }
   fVariableImportance.clear();
   fResiduals.clear();
   fLossFunctionEventInfo.clear();
   // now done in "InitEventSample" which is called in "Train"
   // reset all previously stored/accumulated BOOST weights in the event sample
   //for (UInt_t iev=0; iev<fEventSample.size(); iev++) fEventSample[iev]->SetBoostWeight(1.);
   if (Data()) Data()->DeleteResults(GetMethodName(), Types::kTraining, GetAnalysisType());
   Log() << kDEBUG << " successfully(?) reset the method " << Endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.
///
///  - Note: fEventSample and ValidationSample are already deleted at the end of TRAIN
///         When they are not used anymore

TMVA::MethodBDT::~MethodBDT( void )
{
   for (UInt_t i=0; i<fForest.size();           i++) delete fForest[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the event sample (i.e. reset the boost-weights... etc).

void TMVA::MethodBDT::InitEventSample( void )
{
   if (!HasTrainingTree()) Log() << kFATAL << "<Init> Data().TrainingTree() is zero pointer" << Endl;

   if (fEventSample.size() > 0) { // do not re-initialise the event sample, just set all boostweights to 1. as if it were untouched
      // reset all previously stored/accumulated BOOST weights in the event sample
      for (UInt_t iev=0; iev<fEventSample.size(); iev++) fEventSample[iev]->SetBoostWeight(1.);
   } else {
      Data()->SetCurrentType(Types::kTraining);
      UInt_t nevents = Data()->GetNTrainingEvents();

      std::vector<const TMVA::Event*> tmpEventSample;
      for (Long64_t ievt=0; ievt<nevents; ievt++) {
         //  const Event *event = new Event(*(GetEvent(ievt)));
         Event* event = new Event( *GetTrainingEvent(ievt) );
         tmpEventSample.push_back(event);
      }

      if (!DoRegression()) DeterminePreselectionCuts(tmpEventSample);
      else fDoPreselection = kFALSE; // just to make sure...

      for (UInt_t i=0; i<tmpEventSample.size(); i++) delete tmpEventSample[i];


      Bool_t firstNegWeight=kTRUE;
      Bool_t firstZeroWeight=kTRUE;
      for (Long64_t ievt=0; ievt<nevents; ievt++) {
         //         const Event *event = new Event(*(GetEvent(ievt)));
         // const Event* event = new Event( *GetTrainingEvent(ievt) );
         Event* event = new Event( *GetTrainingEvent(ievt) );
         if (fDoPreselection){
            if (TMath::Abs(ApplyPreselectionCuts(event)) > 0.05) {
               delete event;
               continue;
            }
         }

         if (event->GetWeight() < 0 && (IgnoreEventsWithNegWeightsInTraining() || fNoNegWeightsInTraining)){
            if (firstNegWeight) {
               Log() << kWARNING << " Note, you have events with negative event weight in the sample, but you've chosen to ignore them" << Endl;
               firstNegWeight=kFALSE;
            }
            delete event;
         }else if (event->GetWeight()==0){
            if (firstZeroWeight) {
               firstZeroWeight = kFALSE;
               Log() << "Events with weight == 0 are going to be simply ignored " << Endl;
            }
            delete event;
         }else{
            if (event->GetWeight() < 0) {
               fTrainWithNegWeights=kTRUE;
               if (firstNegWeight){
                  firstNegWeight = kFALSE;
                  if (fPairNegWeightsGlobal){
                     Log() << kWARNING << "Events with negative event weights are found and "
                           << " will be removed prior to the actual BDT training by global "
                           << " paring (and subsequent annihilation) with positiv weight events"
                           << Endl;
                  }else{
                     Log() << kWARNING << "Events with negative event weights are USED during "
                           << "the BDT training. This might cause problems with small node sizes "
                           << "or with the boosting. Please remove negative events from training "
                           << "using the option *IgnoreEventsWithNegWeightsInTraining* in case you "
                           << "observe problems with the boosting"
                           << Endl;
                  }
               }
            }
            // if fAutomatic == true you need a validation sample to optimize pruning
            if (fAutomatic) {
               Double_t modulo = 1.0/(fFValidationEvents);
               Int_t   imodulo = static_cast<Int_t>( fmod(modulo,1.0) > 0.5 ? ceil(modulo) : floor(modulo) );
               if (ievt % imodulo == 0) fValidationSample.push_back( event );
               else                     fEventSample.push_back( event );
            }
            else {
               fEventSample.push_back(event);
            }
         }
      }

      if (fAutomatic) {
         Log() << kINFO << "<InitEventSample> Internally I use " << fEventSample.size()
               << " for Training  and " << fValidationSample.size()
               << " for Pruning Validation (" << ((Float_t)fValidationSample.size())/((Float_t)fEventSample.size()+fValidationSample.size())*100.0
               << "% of training used for validation)" << Endl;
      }

      // some pre-processing for events with negative weights
      if (fPairNegWeightsGlobal) PreProcessNegativeEventWeights();
   }

   if (DoRegression()) {
      // Regression, no reweighting to do
   } else if (DoMulticlass()) {
      // Multiclass, only gradboost is supported. No reweighting.
   } else if (!fSkipNormalization) {
      // Binary classification.
      Log() << kDEBUG << "\t<InitEventSample> For classification trees, "<< Endl;
      Log() << kDEBUG << " \tthe effective number of backgrounds is scaled to match "<<Endl;
      Log() << kDEBUG << " \tthe signal. Otherwise the first boosting step would do 'just that'!"<<Endl;
      // it does not make sense in decision trees to start with unequal number of signal/background
      // events (weights) .. hence normalize them now (happens otherwise in first 'boosting step'
      // anyway..
      // Also make sure, that the sum_of_weights == sample.size() .. as this is assumed in
      // the DecisionTree to derive a sensible number for "fMinSize" (min.#events in node)
      // that currently is an OR between "weighted" and "unweighted number"
      // I want:
      //    nS + nB = n
      //   a*SW + b*BW = n
      //   (a*SW)/(b*BW) = fSigToBkgFraction
      //
      // ==> b = n/((1+f)BW)  and a = (nf/(1+f))/SW

      Double_t nevents = fEventSample.size();
      Double_t sumSigW=0, sumBkgW=0;
      Int_t    sumSig=0, sumBkg=0;
      for (UInt_t ievt=0; ievt<fEventSample.size(); ievt++) {
         if ((DataInfo().IsSignal(fEventSample[ievt])) ) {
            sumSigW += fEventSample[ievt]->GetWeight();
            sumSig++;
         } else {
            sumBkgW += fEventSample[ievt]->GetWeight();
            sumBkg++;
         }
      }
      if (sumSigW && sumBkgW){
         Double_t normSig = nevents/((1+fSigToBkgFraction)*sumSigW)*fSigToBkgFraction;
         Double_t normBkg = nevents/((1+fSigToBkgFraction)*sumBkgW); ;
         Log() << kDEBUG << "\tre-normalise events such that Sig and Bkg have respective sum of weights = "
               << fSigToBkgFraction << Endl;
         Log() << kDEBUG << "  \tsig->sig*"<<normSig << "ev. bkg->bkg*"<<normBkg << "ev." <<Endl;
         Log() << kHEADER << "#events: (reweighted) sig: "<< sumSigW*normSig << " bkg: " << sumBkgW*normBkg << Endl;
         Log() << kINFO << "#events: (unweighted) sig: "<< sumSig << " bkg: " << sumBkg << Endl;
         for (Long64_t ievt=0; ievt<nevents; ievt++) {
            if ((DataInfo().IsSignal(fEventSample[ievt])) ) fEventSample[ievt]->SetBoostWeight(normSig);
            else                                            fEventSample[ievt]->SetBoostWeight(normBkg);
         }
      }else{
         Log() << kINFO << "--> could not determine scaling factors as either there are " << Endl;
         Log() << kINFO << " no signal events (sumSigW="<<sumSigW<<") or no bkg ev. (sumBkgW="<<sumBkgW<<")"<<Endl;
      }

   }

   fTrainSample = &fEventSample;
   if (fBaggedBoost){
      GetBaggedSubSample(fEventSample);
      fTrainSample = &fSubSample;
   }

   //just for debug purposes..
   /*
     sumSigW=0;
     sumBkgW=0;
     for (UInt_t ievt=0; ievt<fEventSample.size(); ievt++) {
     if ((DataInfo().IsSignal(fEventSample[ievt])) ) sumSigW += fEventSample[ievt]->GetWeight();
     else sumBkgW += fEventSample[ievt]->GetWeight();
     }
     Log() << kWARNING << "sigSumW="<<sumSigW<<"bkgSumW="<<sumBkgW<< Endl;
   */
}

////////////////////////////////////////////////////////////////////////////////
/// O.k. you know there are events with negative event weights. This routine will remove
/// them by pairing them with the closest event(s) of the same event class with positive
/// weights
/// A first attempt is "brute force", I dont' try to be clever using search trees etc,
/// just quick and dirty to see if the result is any good

void TMVA::MethodBDT::PreProcessNegativeEventWeights(){
   Double_t totalNegWeights = 0;
   Double_t totalPosWeights = 0;
   Double_t totalWeights    = 0;
   std::vector<const Event*> negEvents;
   for (UInt_t iev = 0; iev < fEventSample.size(); iev++){
      if (fEventSample[iev]->GetWeight() < 0) {
         totalNegWeights += fEventSample[iev]->GetWeight();
         negEvents.push_back(fEventSample[iev]);
      } else {
         totalPosWeights += fEventSample[iev]->GetWeight();
      }
      totalWeights += fEventSample[iev]->GetWeight();
   }
   if (totalNegWeights == 0 ) {
      Log() << kINFO << "no negative event weights found .. no preprocessing necessary" << Endl;
      return;
   } else {
      Log() << kINFO << "found a total of " << totalNegWeights << " of negative event weights which I am going to try to pair with positive events to annihilate them" << Endl;
      Log() << kINFO << "found a total of " << totalPosWeights << " of events with positive weights" << Endl;
      Log() << kINFO << "--> total sum of weights = " << totalWeights << " = " << totalNegWeights+totalPosWeights << Endl;
   }

   std::vector<TMatrixDSym*>* cov = gTools().CalcCovarianceMatrices( fEventSample, 2);

   TMatrixDSym *invCov;

   for (Int_t i=0; i<2; i++){
      invCov = ((*cov)[i]);
      if ( TMath::Abs(invCov->Determinant()) < 10E-24 ) {
         std::cout << "<MethodBDT::PreProcessNeg...> matrix is almost singular with determinant="
                   << TMath::Abs(invCov->Determinant())
                   << " did you use the variables that are linear combinations or highly correlated?"
                   << std::endl;
      }
      if ( TMath::Abs(invCov->Determinant()) < 10E-120 ) {
         std::cout << "<MethodBDT::PreProcessNeg...> matrix is singular with determinant="
                   << TMath::Abs(invCov->Determinant())
                   << " did you use the variables that are linear combinations?"
                   << std::endl;
      }

      invCov->Invert();
   }



   Log() << kINFO << "Found a total of " << totalNegWeights << " in negative weights out of " << fEventSample.size() << " training events "  << Endl;
   Timer timer(negEvents.size(),"Negative Event paired");
   for (UInt_t nev = 0; nev < negEvents.size(); nev++){
      timer.DrawProgressBar( nev );
      Double_t weight = negEvents[nev]->GetWeight();
      UInt_t  iClassID = negEvents[nev]->GetClass();
      invCov = ((*cov)[iClassID]);
      while (weight < 0){
         // find closest event with positive event weight and "pair" it with the negative event
         // (add their weight) until there is no negative weight anymore
         Int_t iMin=-1;
         Double_t dist, minDist=10E270;
         for (UInt_t iev = 0; iev < fEventSample.size(); iev++){
            if (iClassID==fEventSample[iev]->GetClass() && fEventSample[iev]->GetWeight() > 0){
               dist=0;
               for (UInt_t ivar=0; ivar < GetNvar(); ivar++){
                  for (UInt_t jvar=0; jvar<GetNvar(); jvar++){
                     dist += (negEvents[nev]->GetValue(ivar)-fEventSample[iev]->GetValue(ivar))*
                        (*invCov)[ivar][jvar]*
                        (negEvents[nev]->GetValue(jvar)-fEventSample[iev]->GetValue(jvar));
                  }
               }
               if (dist < minDist) { iMin=iev; minDist=dist;}
            }
         }

         if (iMin > -1) {
            //     std::cout << "Happily pairing .. weight before : " << negEvents[nev]->GetWeight() << " and " << fEventSample[iMin]->GetWeight();
            Double_t newWeight = (negEvents[nev]->GetWeight() + fEventSample[iMin]->GetWeight());
            if (newWeight > 0){
               negEvents[nev]->SetBoostWeight( 0 );
               fEventSample[iMin]->SetBoostWeight( newWeight/fEventSample[iMin]->GetOriginalWeight() );  // note the weight*boostweight should be "newWeight"
            } else {
               negEvents[nev]->SetBoostWeight( newWeight/negEvents[nev]->GetOriginalWeight() ); // note the weight*boostweight should be "newWeight"
               fEventSample[iMin]->SetBoostWeight( 0 );
            }
            //     std::cout << " and afterwards " <<  negEvents[nev]->GetWeight() <<  " and the paired " << fEventSample[iMin]->GetWeight() << " dist="<<minDist<< std::endl;
         } else Log() << kFATAL << "preprocessing didn't find event to pair with the negative weight ... probably a bug" << Endl;
         weight = negEvents[nev]->GetWeight();
      }
   }
   Log() << kINFO << "<Negative Event Pairing> took: " << timer.GetElapsedTime()
         << "                              " << Endl;

   // just check.. now there should be no negative event weight left anymore
   totalNegWeights = 0;
   totalPosWeights = 0;
   totalWeights    = 0;
   Double_t sigWeight=0;
   Double_t bkgWeight=0;
   Int_t    nSig=0;
   Int_t    nBkg=0;

   std::vector<const Event*> newEventSample;

   for (UInt_t iev = 0; iev < fEventSample.size(); iev++){
      if (fEventSample[iev]->GetWeight() < 0) {
         totalNegWeights += fEventSample[iev]->GetWeight();
         totalWeights    += fEventSample[iev]->GetWeight();
      } else {
         totalPosWeights += fEventSample[iev]->GetWeight();
         totalWeights    += fEventSample[iev]->GetWeight();
      }
      if (fEventSample[iev]->GetWeight() > 0) {
         newEventSample.push_back(new Event(*fEventSample[iev]));
         if (fEventSample[iev]->GetClass() == fSignalClass){
            sigWeight += fEventSample[iev]->GetWeight();
            nSig+=1;
         }else{
            bkgWeight += fEventSample[iev]->GetWeight();
            nBkg+=1;
         }
      }
   }
   if (totalNegWeights < 0) Log() << kFATAL << " compensation of negative event weights with positive ones did not work " << totalNegWeights << Endl;

   for (UInt_t i=0; i<fEventSample.size();      i++) delete fEventSample[i];
   fEventSample = newEventSample;

   Log() << kINFO  << " after PreProcessing, the Event sample is left with " << fEventSample.size() << " events (unweighted), all with positive weights, adding up to " << totalWeights << Endl;
   Log() << kINFO  << " nSig="<<nSig << " sigWeight="<<sigWeight <<  " nBkg="<<nBkg << " bkgWeight="<<bkgWeight << Endl;


}

////////////////////////////////////////////////////////////////////////////////
/// Call the Optimizer with the set of parameters and ranges that
/// are meant to be tuned.

std::map<TString,Double_t>  TMVA::MethodBDT::OptimizeTuningParameters(TString fomType, TString fitType)
{
   // fill all the tuning parameters that should be optimized into a map:
   std::map<TString,TMVA::Interval*> tuneParameters;
   std::map<TString,Double_t> tunedParameters;

   // note: the 3rd parameter in the interval is the "number of bins", NOT the stepsize !!
   //       the actual VALUES at (at least for the scan, guess also in GA) are always
   //       read from the middle of the bins. Hence.. the choice of Intervals e.g. for the
   //       MaxDepth, in order to make nice integer values!!!

   // find some reasonable ranges for the optimisation of MinNodeEvents:

   tuneParameters.insert(std::pair<TString,Interval*>("NTrees",         new Interval(10,1000,5))); //  stepsize 50
   tuneParameters.insert(std::pair<TString,Interval*>("MaxDepth",       new Interval(2,4,3)));    // stepsize 1
   tuneParameters.insert(std::pair<TString,Interval*>("MinNodeSize",    new LogInterval(1,30,30)));    //
   //tuneParameters.insert(std::pair<TString,Interval*>("NodePurityLimit",new Interval(.4,.6,3)));   // stepsize .1
   //tuneParameters.insert(std::pair<TString,Interval*>("BaggedSampleFraction",new Interval(.4,.9,6)));   // stepsize .1

   // method-specific parameters
   if        (fBoostType=="AdaBoost"){
      tuneParameters.insert(std::pair<TString,Interval*>("AdaBoostBeta",   new Interval(.2,1.,5)));

   }else if (fBoostType=="Grad"){
      tuneParameters.insert(std::pair<TString,Interval*>("Shrinkage",      new Interval(0.05,0.50,5)));

   }else if (fBoostType=="Bagging" && fRandomisedTrees){
      Int_t min_var  = TMath::FloorNint( GetNvar() * .25 );
      Int_t max_var  = TMath::CeilNint(  GetNvar() * .75 );
      tuneParameters.insert(std::pair<TString,Interval*>("UseNvars",       new Interval(min_var,max_var,4)));

   }

   Log()<<kINFO << " the following BDT parameters will be tuned on the respective *grid*\n"<<Endl;
   std::map<TString,TMVA::Interval*>::iterator it;
   for(it=tuneParameters.begin(); it!= tuneParameters.end(); ++it){
      Log() << kWARNING << it->first << Endl;
      std::ostringstream oss;
      (it->second)->Print(oss);
      Log()<<oss.str();
      Log()<<Endl;
   }

   OptimizeConfigParameters optimize(this, tuneParameters, fomType, fitType);
   tunedParameters=optimize.optimize();

   return tunedParameters;

}

////////////////////////////////////////////////////////////////////////////////
/// Set the tuning parameters according to the argument.

void TMVA::MethodBDT::SetTuneParameters(std::map<TString,Double_t> tuneParameters)
{
   std::map<TString,Double_t>::iterator it;
   for(it=tuneParameters.begin(); it!= tuneParameters.end(); ++it){
      Log() << kWARNING << it->first << " = " << it->second << Endl;
      if (it->first ==  "MaxDepth"       ) SetMaxDepth        ((Int_t)it->second);
      else if (it->first ==  "MinNodeSize"    ) SetMinNodeSize     (it->second);
      else if (it->first ==  "NTrees"         ) SetNTrees          ((Int_t)it->second);
      else if (it->first ==  "NodePurityLimit") SetNodePurityLimit (it->second);
      else if (it->first ==  "AdaBoostBeta"   ) SetAdaBoostBeta    (it->second);
      else if (it->first ==  "Shrinkage"      ) SetShrinkage       (it->second);
      else if (it->first ==  "UseNvars"       ) SetUseNvars        ((Int_t)it->second);
      else if (it->first ==  "BaggedSampleFraction" ) SetBaggedSampleFraction (it->second);
      else Log() << kFATAL << " SetParameter for " << it->first << " not yet implemented " <<Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// BDT training.


void TMVA::MethodBDT::Train()
{
   TMVA::DecisionTreeNode::SetIsTraining(true);

   // fill the STL Vector with the event sample
   // (needs to be done here and cannot be done in "init" as the options need to be
   // known).
   InitEventSample();

   if (fNTrees==0){
      Log() << kERROR << " Zero Decision Trees demanded... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      fNTrees = 1;
   }

   if (fInteractive && fInteractive->NotInitialized()){
     std::vector<TString> titles = {"Boost weight", "Error Fraction"};
     fInteractive->Init(titles);
   }
   fIPyMaxIter = fNTrees;
   fExitFromTraining = false;

   // HHV (it's been here since looong but I really don't know why we cannot handle
   // normalized variables in BDTs...  todo
   if (IsNormalised()) Log() << kFATAL << "\"Normalise\" option cannot be used with BDT; "
                             << "please remove the option from the configuration string, or "
                             << "use \"!Normalise\""
                             << Endl;

   if(DoRegression())
      Log() << kINFO << "Regression Loss Function: "<< fRegressionLossFunctionBDTG->Name() << Endl;

   Log() << kINFO << "Training "<< fNTrees << " Decision Trees ... patience please" << Endl;

   Log() << kDEBUG << "Training with maximal depth = " <<fMaxDepth
         << ", MinNodeEvents=" << fMinNodeEvents
         << ", NTrees="<<fNTrees
         << ", NodePurityLimit="<<fNodePurityLimit
         << ", AdaBoostBeta="<<fAdaBoostBeta
         << Endl;

   // weights applied in boosting
   Int_t nBins;
   Double_t xMin,xMax;
   TString hname = "AdaBooost weight distribution";

   nBins= 100;
   xMin = 0;
   xMax = 30;

   if (DoRegression()) {
      nBins= 100;
      xMin = 0;
      xMax = 1;
      hname="Boost event weights distribution";
   }

   // book monitoring histograms (for AdaBost only)

   TH1* h = new TH1F(Form("%s_BoostWeight",DataInfo().GetName()),hname,nBins,xMin,xMax);
   TH1* nodesBeforePruningVsTree = new TH1I(Form("%s_NodesBeforePruning",DataInfo().GetName()),"nodes before pruning",fNTrees,0,fNTrees);
   TH1* nodesAfterPruningVsTree = new TH1I(Form("%s_NodesAfterPruning",DataInfo().GetName()),"nodes after pruning",fNTrees,0,fNTrees);



   if(!DoMulticlass()){
      Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());

      h->SetXTitle("boost weight");
      results->Store(h, "BoostWeights");


      // Monitor the performance (on TEST sample) versus number of trees
      if (fDoBoostMonitor){
         TH2* boostMonitor = new TH2F("BoostMonitor","ROC Integral Vs iTree",2,0,fNTrees,2,0,1.05);
         boostMonitor->SetXTitle("#tree");
         boostMonitor->SetYTitle("ROC Integral");
         results->Store(boostMonitor, "BoostMonitor");
         TGraph *boostMonitorGraph = new TGraph();
         boostMonitorGraph->SetName("BoostMonitorGraph");
         boostMonitorGraph->SetTitle("ROCIntegralVsNTrees");
         results->Store(boostMonitorGraph, "BoostMonitorGraph");
      }

      // weights applied in boosting vs tree number
      h = new TH1F("BoostWeightVsTree","Boost weights vs tree",fNTrees,0,fNTrees);
      h->SetXTitle("#tree");
      h->SetYTitle("boost weight");
      results->Store(h, "BoostWeightsVsTree");

      // error fraction vs tree number
      h = new TH1F("ErrFractHist","error fraction vs tree number",fNTrees,0,fNTrees);
      h->SetXTitle("#tree");
      h->SetYTitle("error fraction");
      results->Store(h, "ErrorFrac");

      // nNodesBeforePruning vs tree number
      nodesBeforePruningVsTree->SetXTitle("#tree");
      nodesBeforePruningVsTree->SetYTitle("#tree nodes");
      results->Store(nodesBeforePruningVsTree);

      // nNodesAfterPruning vs tree number
      nodesAfterPruningVsTree->SetXTitle("#tree");
      nodesAfterPruningVsTree->SetYTitle("#tree nodes");
      results->Store(nodesAfterPruningVsTree);

   }

   fMonitorNtuple= new TTree("MonitorNtuple","BDT variables");
   fMonitorNtuple->Branch("iTree",&fITree,"iTree/I");
   fMonitorNtuple->Branch("boostWeight",&fBoostWeight,"boostWeight/D");
   fMonitorNtuple->Branch("errorFraction",&fErrorFraction,"errorFraction/D");

   Timer timer( fNTrees, GetName() );
   Int_t nNodesBeforePruningCount = 0;
   Int_t nNodesAfterPruningCount = 0;

   Int_t nNodesBeforePruning = 0;
   Int_t nNodesAfterPruning = 0;

   if(fBoostType=="Grad"){
      InitGradBoost(fEventSample);
   }

   Int_t itree=0;
   Bool_t continueBoost=kTRUE;
   //for (int itree=0; itree<fNTrees; itree++) {

   while (itree < fNTrees && continueBoost){
     if (fExitFromTraining) break;
     fIPyCurrentIter = itree;
      timer.DrawProgressBar( itree );
      // Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());
      // TH1 *hxx = new TH1F(Form("swdist%d",itree),Form("swdist%d",itree),10000,0,15);
      // results->Store(hxx,Form("swdist%d",itree));
      // TH1 *hxy = new TH1F(Form("bwdist%d",itree),Form("bwdist%d",itree),10000,0,15);
      // results->Store(hxy,Form("bwdist%d",itree));
      // for (Int_t iev=0; iev<fEventSample.size(); iev++) {
      //    if (fEventSample[iev]->GetClass()!=0) hxy->Fill((fEventSample[iev])->GetWeight());
      //    else          hxx->Fill((fEventSample[iev])->GetWeight());
      // }

      if(DoMulticlass()){
         if (fBoostType!="Grad"){
            Log() << kFATAL << "Multiclass is currently only supported by gradient boost. "
                  << "Please change boost option accordingly (BoostType=Grad)." << Endl;
         }

         UInt_t nClasses = DataInfo().GetNClasses();
         for (UInt_t i=0;i<nClasses;i++){
            // Careful: If fSepType is nullptr, the tree will be considered a regression tree and
            // use the correct output for gradboost (response rather than yesnoleaf) in checkEvent.
            // See TMVA::MethodBDT::InitGradBoost.
            fForest.push_back( new DecisionTree( fSepType, fMinNodeSize, fNCuts, &(DataInfo()), i,
                                                 fRandomisedTrees, fUseNvars, fUsePoissonNvars, fMaxDepth,
                                                 itree*nClasses+i, fNodePurityLimit, itree*nClasses+1));
            fForest.back()->SetNVars(GetNvar());
            if (fUseFisherCuts) {
               fForest.back()->SetUseFisherCuts();
               fForest.back()->SetMinLinCorrForFisher(fMinLinCorrForFisher);
               fForest.back()->SetUseExclusiveVars(fUseExclusiveVars);
            }
            // the minimum linear correlation between two variables demanded for use in fisher criterion in node splitting

            nNodesBeforePruning = fForest.back()->BuildTree(*fTrainSample);
            Double_t bw = this->Boost(*fTrainSample, fForest.back(),i);
            if (bw > 0) {
               fBoostWeights.push_back(bw);
            }else{
               fBoostWeights.push_back(0);
               Log() << kWARNING << "stopped boosting at itree="<<itree << Endl;
               //               fNTrees = itree+1; // that should stop the boosting
               continueBoost=kFALSE;
            }
         }
      }
      else{

         DecisionTree* dt = new DecisionTree( fSepType, fMinNodeSize, fNCuts, &(DataInfo()), fSignalClass,
                                              fRandomisedTrees, fUseNvars, fUsePoissonNvars, fMaxDepth,
                                              itree, fNodePurityLimit, itree);

         fForest.push_back(dt);
         fForest.back()->SetNVars(GetNvar());
         if (fUseFisherCuts) {
            fForest.back()->SetUseFisherCuts();
            fForest.back()->SetMinLinCorrForFisher(fMinLinCorrForFisher);
            fForest.back()->SetUseExclusiveVars(fUseExclusiveVars);
         }

         nNodesBeforePruning = fForest.back()->BuildTree(*fTrainSample);

         if (fUseYesNoLeaf && !DoRegression() && fBoostType!="Grad") { // remove leaf nodes where both daughter nodes are of same type
            nNodesBeforePruning = fForest.back()->CleanTree();
         }

         nNodesBeforePruningCount += nNodesBeforePruning;
         nodesBeforePruningVsTree->SetBinContent(itree+1,nNodesBeforePruning);

         fForest.back()->SetPruneMethod(fPruneMethod); // set the pruning method for the tree
         fForest.back()->SetPruneStrength(fPruneStrength); // set the strength parameter

         std::vector<const Event*> * validationSample = NULL;
         if(fAutomatic) validationSample = &fValidationSample;
         Double_t bw = this->Boost(*fTrainSample, fForest.back());
         if (bw > 0) {
            fBoostWeights.push_back(bw);
         }else{
            fBoostWeights.push_back(0);
            Log() << kWARNING << "stopped boosting at itree="<<itree << Endl;
            continueBoost=kFALSE;
         }

         // if fAutomatic == true, pruneStrength will be the optimal pruning strength
         // determined by the pruning algorithm; otherwise, it is simply the strength parameter
         // set by the user
         if  (fPruneMethod != DecisionTree::kNoPruning) fForest.back()->PruneTree(validationSample);

         if (fUseYesNoLeaf && !DoRegression() && fBoostType!="Grad"){ // remove leaf nodes where both daughter nodes are of same type
            fForest.back()->CleanTree();
         }
         nNodesAfterPruning = fForest.back()->GetNNodes();
         nNodesAfterPruningCount += nNodesAfterPruning;
         nodesAfterPruningVsTree->SetBinContent(itree+1,nNodesAfterPruning);

         if (fInteractive){
           fInteractive->AddPoint(itree, fBoostWeight, fErrorFraction);
         }
         fITree = itree;
         fMonitorNtuple->Fill();
         if (fDoBoostMonitor){
            if (! DoRegression() ){
               if (  itree==fNTrees-1 ||  (!(itree%500)) ||
                     (!(itree%250) && itree <1000)||
                     (!(itree%100) && itree < 500)||
                     (!(itree%50)  && itree < 250)||
                     (!(itree%25)  && itree < 150)||
                     (!(itree%10)  && itree <  50)||
                     (!(itree%5)   && itree <  20)
                     ) BoostMonitor(itree);
            }
         }
      }
      itree++;
   }

   // get elapsed time
   Log() << kDEBUG << "\t<Train> elapsed time: " << timer.GetElapsedTime()
         << "                              " << Endl;
   if (fPruneMethod == DecisionTree::kNoPruning) {
      Log() << kDEBUG << "\t<Train> average number of nodes (w/o pruning) : "
            << nNodesBeforePruningCount/GetNTrees() << Endl;
   }
   else {
      Log() << kDEBUG << "\t<Train> average number of nodes before/after pruning : "
            << nNodesBeforePruningCount/GetNTrees() << " / "
            << nNodesAfterPruningCount/GetNTrees()
            << Endl;
   }
   TMVA::DecisionTreeNode::SetIsTraining(false);


   // reset all previously stored/accumulated BOOST weights in the event sample
   //   for (UInt_t iev=0; iev<fEventSample.size(); iev++) fEventSample[iev]->SetBoostWeight(1.);
   Log() << kDEBUG << "Now I delete the privat data sample"<< Endl;
   for (UInt_t i=0; i<fEventSample.size();      i++) delete fEventSample[i];
   for (UInt_t i=0; i<fValidationSample.size(); i++) delete fValidationSample[i];
   fEventSample.clear();
   fValidationSample.clear();

   if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
   ExitFromTraining();
}


////////////////////////////////////////////////////////////////////////////////
/// Returns MVA value: -1 for background, 1 for signal.

Double_t TMVA::MethodBDT::GetGradBoostMVA(const TMVA::Event* e, UInt_t nTrees)
{
   Double_t sum=0;
   for (UInt_t itree=0; itree<nTrees; itree++) {
      //loop over all trees in forest
      sum += fForest[itree]->CheckEvent(e,kFALSE);

   }
   return 2.0/(1.0+exp(-2.0*sum))-1; //MVA output between -1 and 1
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate residual for all events.

void TMVA::MethodBDT::UpdateTargets(std::vector<const TMVA::Event*>& eventSample, UInt_t cls)
{
   if (DoMulticlass()) {
      UInt_t nClasses = DataInfo().GetNClasses();
      Bool_t isLastClass = (cls == nClasses - 1);

      #ifdef R__USE_IMT
      //
      // This is the multi-threaded multiclass version
      //
      // Note: we only need to update the predicted probabilities every
      // `nClasses` tree. Let's call a set of `nClasses` trees a "round". Thus
      // the algortihm is split in two parts `update_residuals` and
      // `update_residuals_last` where the latter is inteded to be run instead
      // of the former for the last tree in a "round".
      //
      std::map<const TMVA::Event *, std::vector<double>> & residuals = this->fResiduals;
      DecisionTree & lastTree = *(this->fForest.back());

      auto update_residuals = [&residuals, &lastTree, cls](const TMVA::Event * e) {
         residuals[e].at(cls) += lastTree.CheckEvent(e, kFALSE);
      };

      auto update_residuals_last = [&residuals, &lastTree, cls, nClasses](const TMVA::Event * e) {
         residuals[e].at(cls) += lastTree.CheckEvent(e, kFALSE);

         auto &residualsThisEvent = residuals[e];

         std::vector<Double_t> expCache(nClasses, 0.0);
         std::transform(residualsThisEvent.begin(),
                        residualsThisEvent.begin() + nClasses,
                        expCache.begin(), [](Double_t d) { return exp(d); });

         Double_t exp_sum = std::accumulate(expCache.begin(),
                                            expCache.begin() + nClasses,
                                            0.0);

         for (UInt_t i = 0; i < nClasses; i++) {
            Double_t p_cls = expCache[i] / exp_sum;

            Double_t res = (e->GetClass() == i) ? (1.0 - p_cls) : (-p_cls);
            const_cast<TMVA::Event *>(e)->SetTarget(i, res);
         }
      };

      if (isLastClass) {
         TMVA::Config::Instance().GetThreadExecutor()
                                 .Foreach(update_residuals_last, eventSample);
      } else {
         TMVA::Config::Instance().GetThreadExecutor()
                                 .Foreach(update_residuals, eventSample);
      }
      #else
      //
      // Single-threaded multiclass version
      //
      std::vector<Double_t> expCache;
      if (isLastClass) {
         expCache.resize(nClasses);
      }

      for (auto e : eventSample) {
         fResiduals[e].at(cls) += fForest.back()->CheckEvent(e, kFALSE);
         if (isLastClass) {
            auto &residualsThisEvent = fResiduals[e];
            std::transform(residualsThisEvent.begin(),
                           residualsThisEvent.begin() + nClasses,
                           expCache.begin(), [](Double_t d) { return exp(d); });

            Double_t exp_sum = std::accumulate(expCache.begin(),
                                               expCache.begin() + nClasses,
                                               0.0);

            for (UInt_t i = 0; i < nClasses; i++) {
               Double_t p_cls = expCache[i] / exp_sum;

               Double_t res = (e->GetClass() == i) ? (1.0 - p_cls) : (-p_cls);
               const_cast<TMVA::Event *>(e)->SetTarget(i, res);
            }
         }
      }
      #endif
   } else {
      std::map<const TMVA::Event *, std::vector<double>> & residuals = this->fResiduals;
      DecisionTree & lastTree = *(this->fForest.back());

      UInt_t signalClass = DataInfo().GetSignalClassIndex();

      #ifdef R__USE_IMT
      auto update_residuals = [&residuals, &lastTree, signalClass](const TMVA::Event * e) {
         double & residualAt0 = residuals[e].at(0);
         residualAt0 += lastTree.CheckEvent(e, kFALSE);

         Double_t p_sig = 1.0 / (1.0 + exp(-2.0 * residualAt0));
         Double_t res = ((e->GetClass() == signalClass) ? (1.0 - p_sig) : (-p_sig));

         const_cast<TMVA::Event *>(e)->SetTarget(0, res);
      };

      TMVA::Config::Instance().GetThreadExecutor()
                              .Foreach(update_residuals, eventSample);
      #else
      for (auto e : eventSample) {
         double & residualAt0 = residuals[e].at(0);
         residualAt0 += lastTree.CheckEvent(e, kFALSE);

         Double_t p_sig = 1.0 / (1.0 + exp(-2.0 * residualAt0));
         Double_t res = ((e->GetClass() == signalClass) ? (1.0 - p_sig) : (-p_sig));

         const_cast<TMVA::Event *>(e)->SetTarget(0, res);
      }
      #endif
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Calculate residuals for all events and update targets for next iter.
///
/// \param[in] eventSample The collection of events currently under training.
/// \param[in] first Should be true when called before the first boosting
///                  iteration has been run
///
void TMVA::MethodBDT::UpdateTargetsRegression(std::vector<const TMVA::Event*>& eventSample, Bool_t first)
{
   if (!first) {
#ifdef R__USE_IMT
      UInt_t nPartitions = TMVA::Config::Instance().GetThreadExecutor().GetPoolSize();
      auto seeds = ROOT::TSeqU(nPartitions);

      // need a lambda function to pass to TThreadExecutor::MapReduce
      auto f = [this, &nPartitions](UInt_t partition = 0) -> Int_t {
         Int_t start = 1.0 * partition / nPartitions * this->fEventSample.size();
         Int_t end = (partition + 1.0) / nPartitions * this->fEventSample.size();

         for (Int_t i = start; i < end; ++i) {
            const TMVA::Event *e = fEventSample[i];
            LossFunctionEventInfo & lossInfo = fLossFunctionEventInfo.at(e);
            lossInfo.predictedValue += fForest.back()->CheckEvent(e, kFALSE);
         }

         return 0;
      };

      TMVA::Config::Instance().GetThreadExecutor().Map(f, seeds);
#else
      for (const TMVA::Event *e : fEventSample) {
         LossFunctionEventInfo & lossInfo = fLossFunctionEventInfo.at(e);
         lossInfo.predictedValue += fForest.back()->CheckEvent(e, kFALSE);
      }
#endif
   }

   // NOTE: Set targets are also parallelised internally
   fRegressionLossFunctionBDTG->SetTargets(eventSample, fLossFunctionEventInfo);

}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the desired response value for each region.

Double_t TMVA::MethodBDT::GradBoost(std::vector<const TMVA::Event*>& eventSample, DecisionTree *dt, UInt_t cls)
{
   struct LeafInfo {
      Double_t sumWeightTarget = 0;
      Double_t sum2 = 0;
   };

   std::unordered_map<TMVA::DecisionTreeNode*, LeafInfo> leaves;
   for (auto e : eventSample) {
      Double_t weight = e->GetWeight();
      TMVA::DecisionTreeNode* node = dt->GetEventNode(*e);
      auto &v = leaves[node];
      auto target = e->GetTarget(cls);
      v.sumWeightTarget += target * weight;
      v.sum2 += fabs(target) * (1.0 - fabs(target)) * weight;
   }
   for (auto &iLeave : leaves) {
      constexpr auto minValue = 1e-30;
      if (iLeave.second.sum2 < minValue) {
         iLeave.second.sum2 = minValue;
      }
      const Double_t K = DataInfo().GetNClasses();
      iLeave.first->SetResponse(fShrinkage * (K - 1) / K * iLeave.second.sumWeightTarget / iLeave.second.sum2);
   }

   //call UpdateTargets before next tree is grown

   DoMulticlass() ? UpdateTargets(fEventSample, cls) : UpdateTargets(fEventSample);
   return 1; //trees all have the same weight
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation of M_TreeBoost using any loss function as described by Friedman 1999.

Double_t TMVA::MethodBDT::GradBoostRegression(std::vector<const TMVA::Event*>& eventSample, DecisionTree *dt )
{
   // get the vector of events for each terminal so that we can calculate the constant fit value in each
   // terminal node
   // #### Not sure how many events are in each node in advance, so I can't parallelize this easily
   std::map<TMVA::DecisionTreeNode*,vector< TMVA::LossFunctionEventInfo > > leaves;
   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      TMVA::DecisionTreeNode* node = dt->GetEventNode(*(*e));
      (leaves[node]).push_back(fLossFunctionEventInfo[*e]);
   }

   // calculate the constant fit for each terminal node based upon the events in the node
   // node (iLeave->first), vector of event information (iLeave->second)
   // #### could parallelize this and do the leaves at the same time, but this doesn't take very long compared
   // #### to the other processes
   for (std::map<TMVA::DecisionTreeNode*,vector< TMVA::LossFunctionEventInfo > >::iterator iLeave=leaves.begin();
        iLeave!=leaves.end();++iLeave){
      Double_t fit = fRegressionLossFunctionBDTG->Fit(iLeave->second);
      (iLeave->first)->SetResponse(fShrinkage*fit);
   }

   UpdateTargetsRegression(*fTrainSample);

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize targets for first tree.

void TMVA::MethodBDT::InitGradBoost( std::vector<const TMVA::Event*>& eventSample)
{
   // Should get rid of this line. It's just for debugging.
   //std::sort(eventSample.begin(), eventSample.end(), [](const TMVA::Event* a, const TMVA::Event* b){
   //                                     return (a->GetTarget(0) < b->GetTarget(0)); });
   fSepType=NULL; //set fSepType to NULL (regression trees are used for both classification an regression)
   if(DoRegression()){
      for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
         fLossFunctionEventInfo[*e]= TMVA::LossFunctionEventInfo((*e)->GetTarget(0), 0, (*e)->GetWeight());
      }

      fRegressionLossFunctionBDTG->Init(fLossFunctionEventInfo, fBoostWeights);
      UpdateTargetsRegression(*fTrainSample,kTRUE);

      return;
   }
   else if(DoMulticlass()){
      UInt_t nClasses = DataInfo().GetNClasses();
      for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
         for (UInt_t i=0;i<nClasses;i++){
            //Calculate initial residua, assuming equal probability for all classes
            Double_t r = (*e)->GetClass()==i?(1-1.0/nClasses):(-1.0/nClasses);
            const_cast<TMVA::Event*>(*e)->SetTarget(i,r);
            fResiduals[*e].push_back(0);
         }
      }
   }
   else{
      for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
         Double_t r = (DataInfo().IsSignal(*e)?1:0)-0.5; //Calculate initial residua
         const_cast<TMVA::Event*>(*e)->SetTarget(0,r);
         fResiduals[*e].push_back(0);
      }
   }

}
////////////////////////////////////////////////////////////////////////////////
/// Test the tree quality.. in terms of Misclassification.

Double_t TMVA::MethodBDT::TestTreeQuality( DecisionTree *dt )
{
   Double_t ncorrect=0, nfalse=0;
   for (UInt_t ievt=0; ievt<fValidationSample.size(); ievt++) {
      Bool_t isSignalType= (dt->CheckEvent(fValidationSample[ievt]) > fNodePurityLimit ) ? 1 : 0;

      if (isSignalType == (DataInfo().IsSignal(fValidationSample[ievt])) ) {
         ncorrect += fValidationSample[ievt]->GetWeight();
      }
      else{
         nfalse += fValidationSample[ievt]->GetWeight();
      }
   }

   return  ncorrect / (ncorrect + nfalse);
}

////////////////////////////////////////////////////////////////////////////////
/// Apply the boosting algorithm (the algorithm is selecte via the the "option" given
/// in the constructor. The return value is the boosting weight.

Double_t TMVA::MethodBDT::Boost( std::vector<const TMVA::Event*>& eventSample, DecisionTree *dt, UInt_t cls )
{
   Double_t returnVal=-1;

   if      (fBoostType=="AdaBoost")    returnVal = this->AdaBoost  (eventSample, dt);
   else if (fBoostType=="AdaCost")     returnVal = this->AdaCost   (eventSample, dt);
   else if (fBoostType=="Bagging")     returnVal = this->Bagging   ( );
   else if (fBoostType=="RegBoost")    returnVal = this->RegBoost  (eventSample, dt);
   else if (fBoostType=="AdaBoostR2")  returnVal = this->AdaBoostR2(eventSample, dt);
   else if (fBoostType=="Grad"){
      if(DoRegression())
         returnVal = this->GradBoostRegression(eventSample, dt);
      else if(DoMulticlass())
         returnVal = this->GradBoost (eventSample, dt, cls);
      else
         returnVal = this->GradBoost (eventSample, dt);
   }
   else {
      Log() << kINFO << GetOptions() << Endl;
      Log() << kFATAL << "<Boost> unknown boost option " << fBoostType<< " called" << Endl;
   }

   if (fBaggedBoost){
      GetBaggedSubSample(fEventSample);
   }


   return returnVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the ROCIntegral vs Itree from the testSample for the monitoring plots
/// during the training .. but using the testing events

void TMVA::MethodBDT::BoostMonitor(Int_t iTree)
{
   Results* results = Data()->GetResults(GetMethodName(),Types::kTraining, Types::kMaxAnalysisType);

   TH1F *tmpS = new TH1F( "tmpS", "",     100 , -1., 1.00001 );
   TH1F *tmpB = new TH1F( "tmpB", "",     100 , -1., 1.00001 );
   TH1F *tmp;


   UInt_t signalClassNr = DataInfo().GetClassInfo("Signal")->GetNumber();

   // const std::vector<Event*> events=Data()->GetEventCollection(Types::kTesting);
   // //   fMethod->GetTransformationHandler().CalcTransformations(fMethod->Data()->GetEventCollection(Types::kTesting));
   // for (UInt_t iev=0; iev < events.size() ; iev++){
   //    if (events[iev]->GetClass() == signalClassNr) tmp=tmpS;
   //    else                                          tmp=tmpB;
   //    tmp->Fill(PrivateGetMvaValue(*(events[iev])),events[iev]->GetWeight());
   // }

   UInt_t nevents = Data()->GetNTestEvents();
   for (UInt_t iev=0; iev < nevents; iev++){
      const Event* event = GetTestingEvent(iev);

      if (event->GetClass() == signalClassNr) {tmp=tmpS;}
      else                                    {tmp=tmpB;}
      tmp->Fill(PrivateGetMvaValue(event),event->GetWeight());
   }
   Double_t max=1;

   std::vector<TH1F*> hS;
   std::vector<TH1F*> hB;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++){
      hS.push_back(new TH1F(Form("SigVar%dAtTree%d",ivar,iTree),Form("SigVar%dAtTree%d",ivar,iTree),100,DataInfo().GetVariableInfo(ivar).GetMin(),DataInfo().GetVariableInfo(ivar).GetMax()));
      hB.push_back(new TH1F(Form("BkgVar%dAtTree%d",ivar,iTree),Form("BkgVar%dAtTree%d",ivar,iTree),100,DataInfo().GetVariableInfo(ivar).GetMin(),DataInfo().GetVariableInfo(ivar).GetMax()));
      results->Store(hS.back(),hS.back()->GetTitle());
      results->Store(hB.back(),hB.back()->GetTitle());
   }


   for (UInt_t iev=0; iev < fEventSample.size(); iev++){
      if (fEventSample[iev]->GetBoostWeight() > max) max = 1.01*fEventSample[iev]->GetBoostWeight();
   }
   TH1F *tmpBoostWeightsS = new TH1F(Form("BoostWeightsInTreeS%d",iTree),Form("BoostWeightsInTreeS%d",iTree),100,0.,max);
   TH1F *tmpBoostWeightsB = new TH1F(Form("BoostWeightsInTreeB%d",iTree),Form("BoostWeightsInTreeB%d",iTree),100,0.,max);
   results->Store(tmpBoostWeightsS,tmpBoostWeightsS->GetTitle());
   results->Store(tmpBoostWeightsB,tmpBoostWeightsB->GetTitle());

   TH1F *tmpBoostWeights;
   std::vector<TH1F*> *h;

   for (UInt_t iev=0; iev < fEventSample.size(); iev++){
      if (fEventSample[iev]->GetClass() == signalClassNr) {
         tmpBoostWeights=tmpBoostWeightsS;
         h=&hS;
      }else{
         tmpBoostWeights=tmpBoostWeightsB;
         h=&hB;
      }
      tmpBoostWeights->Fill(fEventSample[iev]->GetBoostWeight());
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++){
         (*h)[ivar]->Fill(fEventSample[iev]->GetValue(ivar),fEventSample[iev]->GetWeight());
      }
   }


   TMVA::PDF *sig = new TMVA::PDF( " PDF Sig", tmpS, TMVA::PDF::kSpline3 );
   TMVA::PDF *bkg = new TMVA::PDF( " PDF Bkg", tmpB, TMVA::PDF::kSpline3 );


   TGraph*  gr=results->GetGraph("BoostMonitorGraph");
   Int_t nPoints = gr->GetN();
   gr->Set(nPoints+1);
   gr->SetPoint(nPoints,(Double_t)iTree+1,GetROCIntegral(sig,bkg));

   tmpS->Delete();
   tmpB->Delete();

   delete sig;
   delete bkg;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// The AdaBoost implementation.
/// a new training sample is generated by weighting
/// events that are misclassified by the decision tree. The weight
/// applied is \f$ w = \frac{(1-err)}{err} \f$ or more general:
///            \f$ w = (\frac{(1-err)}{err})^\beta \f$
/// where \f$err\f$ is the fraction of misclassified events in the tree ( <0.5 assuming
/// demanding the that previous selection was better than random guessing)
/// and "beta" being a free parameter (standard: beta = 1) that modifies the
/// boosting.

Double_t TMVA::MethodBDT::AdaBoost( std::vector<const TMVA::Event*>& eventSample, DecisionTree *dt )
{
   Double_t err=0, sumGlobalw=0, sumGlobalwfalse=0, sumGlobalwfalse2=0;

   std::vector<Double_t> sumw(DataInfo().GetNClasses(),0); //for individually re-scaling  each class

   Double_t maxDev=0;
   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      Double_t w = (*e)->GetWeight();
      sumGlobalw += w;
      UInt_t iclass=(*e)->GetClass();
      sumw[iclass] += w;

      if ( DoRegression() ) {
         Double_t tmpDev = TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) );
         sumGlobalwfalse += w * tmpDev;
         sumGlobalwfalse2 += w * tmpDev*tmpDev;
         if (tmpDev > maxDev) maxDev = tmpDev;
      }else{

         if (fUseYesNoLeaf){
            Bool_t isSignalType = (dt->CheckEvent(*e,fUseYesNoLeaf) > fNodePurityLimit );
            if (!(isSignalType == DataInfo().IsSignal(*e))) {
               sumGlobalwfalse+= w;
            }
         }else{
            Double_t dtoutput = (dt->CheckEvent(*e,fUseYesNoLeaf) - 0.5)*2.;
            Int_t    trueType;
            if (DataInfo().IsSignal(*e)) trueType = 1;
            else trueType = -1;
            sumGlobalwfalse+= w*trueType*dtoutput;
         }
      }
   }

   err = sumGlobalwfalse/sumGlobalw ;
   if ( DoRegression() ) {
      //if quadratic loss:
      if (fAdaBoostR2Loss=="linear"){
         err = sumGlobalwfalse/maxDev/sumGlobalw ;
      }
      else if (fAdaBoostR2Loss=="quadratic"){
         err = sumGlobalwfalse2/maxDev/maxDev/sumGlobalw ;
      }
      else if (fAdaBoostR2Loss=="exponential"){
         err = 0;
         for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
            Double_t w = (*e)->GetWeight();
            Double_t  tmpDev = TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) );
            err += w * (1 - exp (-tmpDev/maxDev)) / sumGlobalw;
         }

      }
      else {
         Log() << kFATAL << " you've chosen a Loss type for Adaboost other than linear, quadratic or exponential "
               << " namely " << fAdaBoostR2Loss << "\n"
               << "and this is not implemented... a typo in the options ??" <<Endl;
      }
   }

   Log() << kDEBUG << "BDT AdaBoos  wrong/all: " << sumGlobalwfalse << "/" << sumGlobalw << Endl;


   Double_t newSumGlobalw=0;
   std::vector<Double_t> newSumw(sumw.size(),0);

   Double_t boostWeight=1.;
   if (err >= 0.5 && fUseYesNoLeaf) { // sanity check ... should never happen as otherwise there is apparently
      // something odd with the assignment of the leaf nodes (rem: you use the training
      // events for this determination of the error rate)
      if (dt->GetNNodes() == 1){
         Log() << kERROR << " YOUR tree has only 1 Node... kind of a funny *tree*. I cannot "
               << "boost such a thing... if after 1 step the error rate is == 0.5"
               << Endl
               << "please check why this happens, maybe too many events per node requested ?"
               << Endl;

      }else{
         Log() << kERROR << " The error rate in the BDT boosting is > 0.5. ("<< err
               << ") That should not happen, please check your code (i.e... the BDT code), I "
               << " stop boosting here" <<  Endl;
         return -1;
      }
      err = 0.5;
   } else if (err < 0) {
      Log() << kERROR << " The error rate in the BDT boosting is < 0. That can happen"
            << " due to improper treatment of negative weights in a Monte Carlo.. (if you have"
            << " an idea on how to do it in a better way, please let me know (Helge.Voss@cern.ch)"
            << " for the time being I set it to its absolute value.. just to continue.." <<  Endl;
      err = TMath::Abs(err);
   }
   if (fUseYesNoLeaf)
      boostWeight = TMath::Log((1.-err)/err)*fAdaBoostBeta;
   else
      boostWeight = TMath::Log((1.+err)/(1-err))*fAdaBoostBeta;


   Log() << kDEBUG << "BDT AdaBoos  wrong/all: " << sumGlobalwfalse << "/" << sumGlobalw << " 1-err/err="<<boostWeight<< " log.."<<TMath::Log(boostWeight)<<Endl;

   Results* results = Data()->GetResults(GetMethodName(),Types::kTraining, Types::kMaxAnalysisType);


   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {

      if (fUseYesNoLeaf||DoRegression()){
         if ((!( (dt->CheckEvent(*e,fUseYesNoLeaf) > fNodePurityLimit ) == DataInfo().IsSignal(*e))) || DoRegression()) {
            Double_t boostfactor = TMath::Exp(boostWeight);

            if (DoRegression()) boostfactor = TMath::Power(1/boostWeight,(1.-TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) )/maxDev ) );
            if ( (*e)->GetWeight() > 0 ){
               (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);
               // Helge change back            (*e)->ScaleBoostWeight(boostfactor);
               if (DoRegression()) results->GetHist("BoostWeights")->Fill(boostfactor);
            } else {
               if ( fInverseBoostNegWeights )(*e)->ScaleBoostWeight( 1. / boostfactor); // if the original event weight is negative, and you want to "increase" the events "positive" influence, you'd rather make the event weight "smaller" in terms of it's absolute value while still keeping it something "negative"
               else (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);

            }
         }

      }else{
         Double_t dtoutput = (dt->CheckEvent(*e,fUseYesNoLeaf) - 0.5)*2.;
         Int_t    trueType;
         if (DataInfo().IsSignal(*e)) trueType = 1;
         else trueType = -1;
         Double_t boostfactor = TMath::Exp(-1*boostWeight*trueType*dtoutput);

         if ( (*e)->GetWeight() > 0 ){
            (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);
            // Helge change back            (*e)->ScaleBoostWeight(boostfactor);
            if (DoRegression()) results->GetHist("BoostWeights")->Fill(boostfactor);
         } else {
            if ( fInverseBoostNegWeights )(*e)->ScaleBoostWeight( 1. / boostfactor); // if the original event weight is negative, and you want to "increase" the events "positive" influence, you'd rather make the event weight "smaller" in terms of it's absolute value while still keeping it something "negative"
            else (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);
         }
      }
      newSumGlobalw+=(*e)->GetWeight();
      newSumw[(*e)->GetClass()] += (*e)->GetWeight();
   }


   //   Double_t globalNormWeight=sumGlobalw/newSumGlobalw;
   Double_t globalNormWeight=( (Double_t) eventSample.size())/newSumGlobalw;
   Log() << kDEBUG << "new Nsig="<<newSumw[0]*globalNormWeight << " new Nbkg="<<newSumw[1]*globalNormWeight << Endl;


   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      //      if (fRenormByClass) (*e)->ScaleBoostWeight( normWeightByClass[(*e)->GetClass()] );
      //      else                (*e)->ScaleBoostWeight( globalNormWeight );
      //      else                (*e)->ScaleBoostWeight( globalNormWeight );
      if (DataInfo().IsSignal(*e))(*e)->ScaleBoostWeight( globalNormWeight * fSigToBkgFraction );
      else                (*e)->ScaleBoostWeight( globalNormWeight );
   }

   if (!(DoRegression()))results->GetHist("BoostWeights")->Fill(boostWeight);
   results->GetHist("BoostWeightsVsTree")->SetBinContent(fForest.size(),boostWeight);
   results->GetHist("ErrorFrac")->SetBinContent(fForest.size(),err);

   fBoostWeight = boostWeight;
   fErrorFraction = err;

   return boostWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// The AdaCost boosting algorithm takes a simple cost Matrix  (currently fixed for
/// all events... later could be modified to use individual cost matrices for each
/// events as in the original paper...
///
///                   true_signal true_bkg
///     ----------------------------------
///     sel_signal |   Css         Ctb_ss    Cxx.. in the range [0,1]
///     sel_bkg    |   Cts_sb      Cbb
///
/// and takes this into account when calculating the mis class. cost (former: error fraction):
///
///     err = sum_events ( weight* y_true*y_sel * beta(event)

Double_t TMVA::MethodBDT::AdaCost( vector<const TMVA::Event*>& eventSample, DecisionTree *dt )
{
   Double_t Css = fCss;
   Double_t Cbb = fCbb;
   Double_t Cts_sb = fCts_sb;
   Double_t Ctb_ss = fCtb_ss;

   Double_t err=0, sumGlobalWeights=0, sumGlobalCost=0;

   std::vector<Double_t> sumw(DataInfo().GetNClasses(),0);      //for individually re-scaling  each class

   for (vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      Double_t w = (*e)->GetWeight();
      sumGlobalWeights += w;
      UInt_t iclass=(*e)->GetClass();

      sumw[iclass] += w;

      if ( DoRegression() ) {
         Log() << kFATAL << " AdaCost not implemented for regression"<<Endl;
      }else{

         Double_t dtoutput = (dt->CheckEvent(*e,false) - 0.5)*2.;
         Int_t    trueType;
         Bool_t   isTrueSignal = DataInfo().IsSignal(*e);
         Bool_t   isSelectedSignal = (dtoutput>0);
         if (isTrueSignal) trueType = 1;
         else trueType = -1;

         Double_t cost=0;
         if       (isTrueSignal  && isSelectedSignal)  cost=Css;
         else if  (isTrueSignal  && !isSelectedSignal) cost=Cts_sb;
         else if  (!isTrueSignal  && isSelectedSignal) cost=Ctb_ss;
         else if  (!isTrueSignal && !isSelectedSignal) cost=Cbb;
         else Log() << kERROR << "something went wrong in AdaCost" << Endl;

         sumGlobalCost+= w*trueType*dtoutput*cost;

      }
   }

   if ( DoRegression() ) {
      Log() << kFATAL << " AdaCost not implemented for regression"<<Endl;
   }

   //   Log() << kDEBUG << "BDT AdaBoos  wrong/all: " << sumGlobalCost << "/" << sumGlobalWeights << Endl;
   //      Log() << kWARNING << "BDT AdaBoos  wrong/all: " << sumGlobalCost << "/" << sumGlobalWeights << Endl;
   sumGlobalCost /= sumGlobalWeights;
   //   Log() << kWARNING << "BDT AdaBoos  wrong/all: " << sumGlobalCost << "/" << sumGlobalWeights << Endl;


   Double_t newSumGlobalWeights=0;
   vector<Double_t> newSumClassWeights(sumw.size(),0);

   Double_t boostWeight = TMath::Log((1+sumGlobalCost)/(1-sumGlobalCost)) * fAdaBoostBeta;

   Results* results = Data()->GetResults(GetMethodName(),Types::kTraining, Types::kMaxAnalysisType);

   for (vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      Double_t dtoutput = (dt->CheckEvent(*e,false) - 0.5)*2.;
      Int_t    trueType;
      Bool_t   isTrueSignal = DataInfo().IsSignal(*e);
      Bool_t   isSelectedSignal = (dtoutput>0);
      if (isTrueSignal) trueType = 1;
      else trueType = -1;

      Double_t cost=0;
      if       (isTrueSignal  && isSelectedSignal)  cost=Css;
      else if  (isTrueSignal  && !isSelectedSignal) cost=Cts_sb;
      else if  (!isTrueSignal  && isSelectedSignal) cost=Ctb_ss;
      else if  (!isTrueSignal && !isSelectedSignal) cost=Cbb;
      else Log() << kERROR << "something went wrong in AdaCost" << Endl;

      Double_t boostfactor = TMath::Exp(-1*boostWeight*trueType*dtoutput*cost);
      if (DoRegression())Log() << kFATAL << " AdaCost not implemented for regression"<<Endl;
      if ( (*e)->GetWeight() > 0 ){
         (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);
         // Helge change back            (*e)->ScaleBoostWeight(boostfactor);
         if (DoRegression())Log() << kFATAL << " AdaCost not implemented for regression"<<Endl;
      } else {
         if ( fInverseBoostNegWeights )(*e)->ScaleBoostWeight( 1. / boostfactor); // if the original event weight is negative, and you want to "increase" the events "positive" influence, you'd rather make the event weight "smaller" in terms of it's absolute value while still keeping it something "negative"
      }

      newSumGlobalWeights+=(*e)->GetWeight();
      newSumClassWeights[(*e)->GetClass()] += (*e)->GetWeight();
   }


   //  Double_t globalNormWeight=sumGlobalWeights/newSumGlobalWeights;
   Double_t globalNormWeight=Double_t(eventSample.size())/newSumGlobalWeights;
   Log() << kDEBUG << "new Nsig="<<newSumClassWeights[0]*globalNormWeight << " new Nbkg="<<newSumClassWeights[1]*globalNormWeight << Endl;


   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      // if (fRenormByClass) (*e)->ScaleBoostWeight( normWeightByClass[(*e)->GetClass()] );
      // else                (*e)->ScaleBoostWeight( globalNormWeight );
      if (DataInfo().IsSignal(*e))(*e)->ScaleBoostWeight( globalNormWeight * fSigToBkgFraction );
      else                (*e)->ScaleBoostWeight( globalNormWeight );
   }


   if (!(DoRegression()))results->GetHist("BoostWeights")->Fill(boostWeight);
   results->GetHist("BoostWeightsVsTree")->SetBinContent(fForest.size(),boostWeight);
   results->GetHist("ErrorFrac")->SetBinContent(fForest.size(),err);

   fBoostWeight = boostWeight;
   fErrorFraction = err;


   return boostWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Call it boot-strapping, re-sampling or whatever you like, in the end it is nothing
/// else but applying "random" poisson weights to each event.

Double_t TMVA::MethodBDT::Bagging( )
{
   // this is now done in "MethodBDT::Boost  as it might be used by other boost methods, too
   // GetBaggedSample(eventSample);

   return 1.;  //here as there are random weights for each event, just return a constant==1;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills fEventSample with fBaggedSampleFraction*NEvents random training events.

void TMVA::MethodBDT::GetBaggedSubSample(std::vector<const TMVA::Event*>& eventSample)
{

   Double_t n;
   TRandom3 *trandom   = new TRandom3(100*fForest.size()+1234);

   if (!fSubSample.empty()) fSubSample.clear();

   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      n = trandom->PoissonD(fBaggedSampleFraction);
      for (Int_t i=0;i<n;i++) fSubSample.push_back(*e);
   }

   delete trandom;
   return;

   /*
     UInt_t nevents = fEventSample.size();

     if (!fSubSample.empty()) fSubSample.clear();
     TRandom3 *trandom   = new TRandom3(fForest.size()+1);

     for (UInt_t ievt=0; ievt<nevents; ievt++) { // recreate new random subsample
     if(trandom->Rndm()<fBaggedSampleFraction)
     fSubSample.push_back(fEventSample[ievt]);
     }
     delete trandom;
   */

}

////////////////////////////////////////////////////////////////////////////////
/// A special boosting only for Regression (not implemented).

Double_t TMVA::MethodBDT::RegBoost( std::vector<const TMVA::Event*>& /* eventSample */, DecisionTree* /* dt */ )
{
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Adaption of the AdaBoost to regression problems (see H.Drucker 1997).

Double_t TMVA::MethodBDT::AdaBoostR2( std::vector<const TMVA::Event*>& eventSample, DecisionTree *dt )
{
   if ( !DoRegression() ) Log() << kFATAL << "Somehow you chose a regression boost method for a classification job" << Endl;

   Double_t err=0, sumw=0, sumwfalse=0, sumwfalse2=0;
   Double_t maxDev=0;
   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      Double_t w = (*e)->GetWeight();
      sumw += w;

      Double_t  tmpDev = TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) );
      sumwfalse  += w * tmpDev;
      sumwfalse2 += w * tmpDev*tmpDev;
      if (tmpDev > maxDev) maxDev = tmpDev;
   }

   //if quadratic loss:
   if (fAdaBoostR2Loss=="linear"){
      err = sumwfalse/maxDev/sumw ;
   }
   else if (fAdaBoostR2Loss=="quadratic"){
      err = sumwfalse2/maxDev/maxDev/sumw ;
   }
   else if (fAdaBoostR2Loss=="exponential"){
      err = 0;
      for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
         Double_t w = (*e)->GetWeight();
         Double_t  tmpDev = TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) );
         err += w * (1 - exp (-tmpDev/maxDev)) / sumw;
      }

   }
   else {
      Log() << kFATAL << " you've chosen a Loss type for Adaboost other than linear, quadratic or exponential "
            << " namely " << fAdaBoostR2Loss << "\n"
            << "and this is not implemented... a typo in the options ??" <<Endl;
   }


   if (err >= 0.5) { // sanity check ... should never happen as otherwise there is apparently
      // something odd with the assignment of the leaf nodes (rem: you use the training
      // events for this determination of the error rate)
      if (dt->GetNNodes() == 1){
         Log() << kERROR << " YOUR tree has only 1 Node... kind of a funny *tree*. I cannot "
               << "boost such a thing... if after 1 step the error rate is == 0.5"
               << Endl
               << "please check why this happens, maybe too many events per node requested ?"
               << Endl;

      }else{
         Log() << kERROR << " The error rate in the BDT boosting is > 0.5. ("<< err
               << ") That should not happen, but is possible for regression trees, and"
               << " should trigger a stop for the boosting. please check your code (i.e... the BDT code), I "
               << " stop boosting " <<  Endl;
         return -1;
      }
      err = 0.5;
   } else if (err < 0) {
      Log() << kERROR << " The error rate in the BDT boosting is < 0. That can happen"
            << " due to improper treatment of negative weights in a Monte Carlo.. (if you have"
            << " an idea on how to do it in a better way, please let me know (Helge.Voss@cern.ch)"
            << " for the time being I set it to its absolute value.. just to continue.." <<  Endl;
      err = TMath::Abs(err);
   }

   Double_t boostWeight = err / (1.-err);
   Double_t newSumw=0;

   Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, Types::kMaxAnalysisType);

   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      Double_t boostfactor =  TMath::Power(boostWeight,(1.-TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) )/maxDev ) );
      results->GetHist("BoostWeights")->Fill(boostfactor);
      //      std::cout << "R2  " << boostfactor << "   " << boostWeight << "   " << (1.-TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) )/maxDev)  << std::endl;
      if ( (*e)->GetWeight() > 0 ){
         Float_t newBoostWeight = (*e)->GetBoostWeight() * boostfactor;
         Float_t newWeight = (*e)->GetWeight() * (*e)->GetBoostWeight() * boostfactor;
         if (newWeight == 0) {
            Log() << kINFO << "Weight=    "   <<   (*e)->GetWeight() << Endl;
            Log() << kINFO  << "BoostWeight= " <<   (*e)->GetBoostWeight() << Endl;
            Log() << kINFO  << "boostweight="<<boostWeight << "  err= " <<err << Endl;
            Log() << kINFO  << "NewBoostWeight= " <<   newBoostWeight << Endl;
            Log() << kINFO  << "boostfactor= " <<  boostfactor << Endl;
            Log() << kINFO  << "maxDev     = " <<  maxDev << Endl;
            Log() << kINFO  << "tmpDev     = " <<  TMath::Abs(dt->CheckEvent(*e,kFALSE) - (*e)->GetTarget(0) ) << Endl;
            Log() << kINFO  << "target     = " <<  (*e)->GetTarget(0)  << Endl;
            Log() << kINFO  << "estimate   = " <<  dt->CheckEvent(*e,kFALSE)  << Endl;
         }
         (*e)->SetBoostWeight( newBoostWeight );
         //         (*e)->SetBoostWeight( (*e)->GetBoostWeight() * boostfactor);
      } else {
         (*e)->SetBoostWeight( (*e)->GetBoostWeight() / boostfactor);
      }
      newSumw+=(*e)->GetWeight();
   }

   // re-normalise the weights
   Double_t normWeight =  sumw / newSumw;
   for (std::vector<const TMVA::Event*>::const_iterator e=eventSample.begin(); e!=eventSample.end();++e) {
      //Helge    (*e)->ScaleBoostWeight( sumw/newSumw);
      // (*e)->ScaleBoostWeight( normWeight);
      (*e)->SetBoostWeight( (*e)->GetBoostWeight() * normWeight );
   }


   results->GetHist("BoostWeightsVsTree")->SetBinContent(fForest.size(),1./boostWeight);
   results->GetHist("ErrorFrac")->SetBinContent(fForest.size(),err);

   fBoostWeight = boostWeight;
   fErrorFraction = err;

   return TMath::Log(1./boostWeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Write weights to XML.

void TMVA::MethodBDT::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");

   if (fDoPreselection){
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++){
         gTools().AddAttr( wght, Form("PreselectionLowBkgVar%d",ivar),      fIsLowBkgCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionLowBkgVar%dValue",ivar), fLowBkgCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionLowSigVar%d",ivar),      fIsLowSigCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionLowSigVar%dValue",ivar), fLowSigCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionHighBkgVar%d",ivar),     fIsHighBkgCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionHighBkgVar%dValue",ivar),fHighBkgCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionHighSigVar%d",ivar),     fIsHighSigCut[ivar]);
         gTools().AddAttr( wght, Form("PreselectionHighSigVar%dValue",ivar),fHighSigCut[ivar]);
      }
   }


   gTools().AddAttr( wght, "NTrees", fForest.size() );
   gTools().AddAttr( wght, "AnalysisType", fForest.back()->GetAnalysisType() );

   for (UInt_t i=0; i< fForest.size(); i++) {
      void* trxml = fForest[i]->AddXMLTo(wght);
      gTools().AddAttr( trxml, "boostWeight", fBoostWeights[i] );
      gTools().AddAttr( trxml, "itree", i );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reads the BDT from the xml file.

void TMVA::MethodBDT::ReadWeightsFromXML(void* parent) {
   UInt_t i;
   for (i=0; i<fForest.size(); i++) delete fForest[i];
   fForest.clear();
   fBoostWeights.clear();

   UInt_t ntrees;
   UInt_t analysisType;
   Float_t boostWeight;


   if (gTools().HasAttr( parent, Form("PreselectionLowBkgVar%d",0))) {
      fIsLowBkgCut.resize(GetNvar());
      fLowBkgCut.resize(GetNvar());
      fIsLowSigCut.resize(GetNvar());
      fLowSigCut.resize(GetNvar());
      fIsHighBkgCut.resize(GetNvar());
      fHighBkgCut.resize(GetNvar());
      fIsHighSigCut.resize(GetNvar());
      fHighSigCut.resize(GetNvar());

      Bool_t tmpBool;
      Double_t tmpDouble;
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++){
         gTools().ReadAttr( parent, Form("PreselectionLowBkgVar%d",ivar), tmpBool);
         fIsLowBkgCut[ivar]=tmpBool;
         gTools().ReadAttr( parent, Form("PreselectionLowBkgVar%dValue",ivar), tmpDouble);
         fLowBkgCut[ivar]=tmpDouble;
         gTools().ReadAttr( parent, Form("PreselectionLowSigVar%d",ivar), tmpBool);
         fIsLowSigCut[ivar]=tmpBool;
         gTools().ReadAttr( parent, Form("PreselectionLowSigVar%dValue",ivar), tmpDouble);
         fLowSigCut[ivar]=tmpDouble;
         gTools().ReadAttr( parent, Form("PreselectionHighBkgVar%d",ivar), tmpBool);
         fIsHighBkgCut[ivar]=tmpBool;
         gTools().ReadAttr( parent, Form("PreselectionHighBkgVar%dValue",ivar), tmpDouble);
         fHighBkgCut[ivar]=tmpDouble;
         gTools().ReadAttr( parent, Form("PreselectionHighSigVar%d",ivar),tmpBool);
         fIsHighSigCut[ivar]=tmpBool;
         gTools().ReadAttr( parent, Form("PreselectionHighSigVar%dValue",ivar), tmpDouble);
         fHighSigCut[ivar]=tmpDouble;
      }
   }


   gTools().ReadAttr( parent, "NTrees", ntrees );

   if(gTools().HasAttr(parent, "TreeType")) { // pre 4.1.0 version
      gTools().ReadAttr( parent, "TreeType", analysisType );
   } else {                                 // from 4.1.0 onwards
      gTools().ReadAttr( parent, "AnalysisType", analysisType );
   }

   void* ch = gTools().GetChild(parent);
   i=0;
   while(ch) {
      fForest.push_back( dynamic_cast<DecisionTree*>( DecisionTree::CreateFromXML(ch, GetTrainingTMVAVersionCode()) ) );
      fForest.back()->SetAnalysisType(Types::EAnalysisType(analysisType));
      fForest.back()->SetTreeID(i++);
      gTools().ReadAttr(ch,"boostWeight",boostWeight);
      fBoostWeights.push_back(boostWeight);
      ch = gTools().GetNextChild(ch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read the weights (BDT coefficients).

void  TMVA::MethodBDT::ReadWeightsFromStream( std::istream& istr )
{
   TString dummy;
   //   Types::EAnalysisType analysisType;
   Int_t analysisType(0);

   // coverity[tainted_data_argument]
   istr >> dummy >> fNTrees;
   Log() << kINFO << "Read " << fNTrees << " Decision trees" << Endl;

   for (UInt_t i=0;i<fForest.size();i++) delete fForest[i];
   fForest.clear();
   fBoostWeights.clear();
   Int_t iTree;
   Double_t boostWeight;
   for (int i=0;i<fNTrees;i++) {
      istr >> dummy >> iTree >> dummy >> boostWeight;
      if (iTree != i) {
         fForest.back()->Print( std::cout );
         Log() << kFATAL << "Error while reading weight file; mismatch iTree="
               << iTree << " i=" << i
               << " dummy " << dummy
               << " boostweight " << boostWeight
               << Endl;
      }
      fForest.push_back( new DecisionTree() );
      fForest.back()->SetAnalysisType(Types::EAnalysisType(analysisType));
      fForest.back()->SetTreeID(i);
      fForest.back()->Read(istr, GetTrainingTMVAVersionCode());
      fBoostWeights.push_back(boostWeight);
   }
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::MethodBDT::GetMvaValue( Double_t* err, Double_t* errUpper ){
   return this->GetMvaValue( err, errUpper, 0 );
}

////////////////////////////////////////////////////////////////////////////////
/// Return the MVA value (range [-1;1]) that classifies the
/// event according to the majority vote from the total number of
/// decision trees.

Double_t TMVA::MethodBDT::GetMvaValue( Double_t* err, Double_t* errUpper, UInt_t useNTrees )
{
   const Event* ev = GetEvent();
   if (fDoPreselection) {
      Double_t val = ApplyPreselectionCuts(ev);
      if (TMath::Abs(val)>0.05) return val;
   }
   return PrivateGetMvaValue(ev, err, errUpper, useNTrees);

}

////////////////////////////////////////////////////////////////////////////////
/// Return the MVA value (range [-1;1]) that classifies the
/// event according to the majority vote from the total number of
/// decision trees.

Double_t TMVA::MethodBDT::PrivateGetMvaValue(const TMVA::Event* ev, Double_t* err, Double_t* errUpper, UInt_t useNTrees )
{
   // cannot determine error
   NoErrorCalc(err, errUpper);

   // allow for the possibility to use less trees in the actual MVA calculation
   // than have been originally trained.
   UInt_t nTrees = fForest.size();

   if (useNTrees > 0 ) nTrees = useNTrees;

   if (fBoostType=="Grad") return GetGradBoostMVA(ev,nTrees);

   Double_t myMVA = 0;
   Double_t norm  = 0;
   for (UInt_t itree=0; itree<nTrees; itree++) {
      //
      myMVA += fBoostWeights[itree] * fForest[itree]->CheckEvent(ev,fUseYesNoLeaf);
      norm  += fBoostWeights[itree];
   }
   return ( norm > std::numeric_limits<double>::epsilon() ) ? myMVA /= norm : 0 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the multiclass MVA response for the BDT classifier.

const std::vector<Float_t>& TMVA::MethodBDT::GetMulticlassValues()
{
   const TMVA::Event *e = GetEvent();
   if (fMulticlassReturnVal == NULL) fMulticlassReturnVal = new std::vector<Float_t>();
   fMulticlassReturnVal->clear();

   UInt_t nClasses = DataInfo().GetNClasses();
   std::vector<Double_t> temp(nClasses);
   auto forestSize = fForest.size();

   #ifdef R__USE_IMT
   std::vector<TMVA::DecisionTree *> forest = fForest;
   auto get_output = [&e, &forest, &temp, forestSize, nClasses](UInt_t iClass) {
      for (UInt_t itree = iClass; itree < forestSize; itree += nClasses) {
         temp[iClass] += forest[itree]->CheckEvent(e, kFALSE);
      }
   };

   TMVA::Config::Instance().GetThreadExecutor()
                           .Foreach(get_output, ROOT::TSeqU(nClasses));
   #else
   // trees 0, nClasses, 2*nClasses, ... belong to class 0
   // trees 1, nClasses+1, 2*nClasses+1, ... belong to class 1 and so forth
   UInt_t classOfTree = 0;
   for (UInt_t itree = 0; itree < forestSize; ++itree) {
      temp[classOfTree] += fForest[itree]->CheckEvent(e, kFALSE);
      if (++classOfTree == nClasses) classOfTree = 0; // cheap modulo
   }
   #endif

   // we want to calculate sum of exp(temp[j] - temp[i]) for all i,j (i!=j)
   // first calculate exp(), then replace minus with division.
   std::transform(temp.begin(), temp.end(), temp.begin(), [](Double_t d){return exp(d);});

   Double_t exp_sum = std::accumulate(temp.begin(), temp.end(), 0.0);

   for (UInt_t i = 0; i < nClasses; i++) {
      Double_t p_cls = temp[i] / exp_sum;
      (*fMulticlassReturnVal).push_back(p_cls);
   }

   return *fMulticlassReturnVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the regression value generated by the BDTs.

const std::vector<Float_t> & TMVA::MethodBDT::GetRegressionValues()
{

   if (fRegressionReturnVal == NULL) fRegressionReturnVal = new std::vector<Float_t>();
   fRegressionReturnVal->clear();

   const Event * ev = GetEvent();
   Event * evT = new Event(*ev);

   Double_t myMVA = 0;
   Double_t norm  = 0;
   if (fBoostType=="AdaBoostR2") {
      // rather than using the weighted average of the tree respones in the forest
      // H.Decker(1997) proposed to use the "weighted median"

      // sort all individual tree responses according to the prediction value
      //   (keep the association to their tree weight)
      // the sum up all the associated weights (starting from the one whose tree
      //   yielded the smalles response) up to the tree "t" at which you've
      //   added enough tree weights to have more than half of the sum of all tree weights.
      // choose as response of the forest that one which belongs to this "t"

      vector< Double_t > response(fForest.size());
      vector< Double_t > weight(fForest.size());
      Double_t           totalSumOfWeights = 0;

      for (UInt_t itree=0; itree<fForest.size(); itree++) {
         response[itree]    = fForest[itree]->CheckEvent(ev,kFALSE);
         weight[itree]      = fBoostWeights[itree];
         totalSumOfWeights += fBoostWeights[itree];
      }

      std::vector< std::vector<Double_t> > vtemp;
      vtemp.push_back( response ); // this is the vector that will get sorted
      vtemp.push_back( weight );
      gTools().UsefulSortAscending( vtemp );

      Int_t t=0;
      Double_t sumOfWeights = 0;
      while (sumOfWeights <= totalSumOfWeights/2.) {
         sumOfWeights += vtemp[1][t];
         t++;
      }

      Double_t rVal=0;
      Int_t    count=0;
      for (UInt_t i= TMath::Max(UInt_t(0),UInt_t(t-(fForest.size()/6)-0.5));
           i< TMath::Min(UInt_t(fForest.size()),UInt_t(t+(fForest.size()/6)+0.5)); i++) {
         count++;
         rVal+=vtemp[0][i];
      }
      //      fRegressionReturnVal->push_back( rVal/Double_t(count));
      evT->SetTarget(0, rVal/Double_t(count) );
   }
   else if(fBoostType=="Grad"){
      for (UInt_t itree=0; itree<fForest.size(); itree++) {
         myMVA += fForest[itree]->CheckEvent(ev,kFALSE);
      }
      //      fRegressionReturnVal->push_back( myMVA+fBoostWeights[0]);
      evT->SetTarget(0, myMVA+fBoostWeights[0] );
   }
   else{
      for (UInt_t itree=0; itree<fForest.size(); itree++) {
         //
         myMVA += fBoostWeights[itree] * fForest[itree]->CheckEvent(ev,kFALSE);
         norm  += fBoostWeights[itree];
      }
      //      fRegressionReturnVal->push_back( ( norm > std::numeric_limits<double>::epsilon() ) ? myMVA /= norm : 0 );
      evT->SetTarget(0, ( norm > std::numeric_limits<double>::epsilon() ) ? myMVA /= norm : 0 );
   }



   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );
   fRegressionReturnVal->push_back( evT2->GetTarget(0) );

   delete evT;


   return *fRegressionReturnVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Here we could write some histograms created during the processing
/// to the output file.

void  TMVA::MethodBDT::WriteMonitoringHistosToFile( void ) const
{
   Log() << kDEBUG << "\tWrite monitoring histograms to file: " << BaseDir()->GetPath() << Endl;

   //Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, Types::kMaxAnalysisType);
   //results->GetStorage()->Write();
   fMonitorNtuple->Write();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the relative variable importance, normalized to all
/// variables together having the importance 1. The importance in
/// evaluated as the total separation-gain that this variable had in
/// the decision trees (weighted by the number of events)

vector< Double_t > TMVA::MethodBDT::GetVariableImportance()
{
   fVariableImportance.resize(GetNvar());
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) {
      fVariableImportance[ivar]=0;
   }
   Double_t  sum=0;
   for (UInt_t itree = 0; itree < GetNTrees(); itree++) {
      std::vector<Double_t> relativeImportance(fForest[itree]->GetVariableImportance());
      for (UInt_t i=0; i< relativeImportance.size(); i++) {
         fVariableImportance[i] +=  fBoostWeights[itree] * relativeImportance[i];
      }
   }

   for (UInt_t ivar=0; ivar< fVariableImportance.size(); ivar++){
      fVariableImportance[ivar] = TMath::Sqrt(fVariableImportance[ivar]);
      sum += fVariableImportance[ivar];
   }
   for (UInt_t ivar=0; ivar< fVariableImportance.size(); ivar++) fVariableImportance[ivar] /= sum;

   return fVariableImportance;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the measure for the variable importance of variable "ivar"
/// which is later used in GetVariableImportance() to calculate the
/// relative variable importances.

Double_t TMVA::MethodBDT::GetVariableImportance( UInt_t ivar )
{
   std::vector<Double_t> relativeImportance = this->GetVariableImportance();
   if (ivar < (UInt_t)relativeImportance.size()) return relativeImportance[ivar];
   else Log() << kFATAL << "<GetVariableImportance> ivar = " << ivar << " is out of range " << Endl;

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute ranking of input variables

const TMVA::Ranking* TMVA::MethodBDT::CreateRanking()
{
   // create the ranking object
   fRanking = new Ranking( GetName(), "Variable Importance" );
   vector< Double_t> importance(this->GetVariableImportance());

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {

      fRanking->AddRank( Rank( GetInputLabel(ivar), importance[ivar] ) );
   }

   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////
/// Get help message text.

void TMVA::MethodBDT::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Boosted Decision Trees are a collection of individual decision" << Endl;
   Log() << "trees which form a multivariate classifier by (weighted) majority " << Endl;
   Log() << "vote of the individual trees. Consecutive decision trees are  " << Endl;
   Log() << "trained using the original training data set with re-weighted " << Endl;
   Log() << "events. By default, the AdaBoost method is employed, which gives " << Endl;
   Log() << "events that were misclassified in the previous tree a larger " << Endl;
   Log() << "weight in the training of the following tree." << Endl;
   Log() << Endl;
   Log() << "Decision trees are a sequence of binary splits of the data sample" << Endl;
   Log() << "using a single discriminant variable at a time. A test event " << Endl;
   Log() << "ending up after the sequence of left-right splits in a final " << Endl;
   Log() << "(\"leaf\") node is classified as either signal or background" << Endl;
   Log() << "depending on the majority type of training events in that node." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "By the nature of the binary splits performed on the individual" << Endl;
   Log() << "variables, decision trees do not deal well with linear correlations" << Endl;
   Log() << "between variables (they need to approximate the linear split in" << Endl;
   Log() << "the two dimensional space by a sequence of splits on the two " << Endl;
   Log() << "variables individually). Hence decorrelation could be useful " << Endl;
   Log() << "to optimise the BDT performance." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The two most important parameters in the configuration are the  " << Endl;
   Log() << "minimal number of events requested by a leaf node as percentage of the " <<Endl;
   Log() << "   number of training events (option \"MinNodeSize\"  replacing the actual number " << Endl;
   Log() << " of events \"nEventsMin\" as given in earlier versions" << Endl;
   Log() << "If this number is too large, detailed features " << Endl;
   Log() << "in the parameter space are hard to be modelled. If it is too small, " << Endl;
   Log() << "the risk to overtrain rises and boosting seems to be less effective" << Endl;
   Log() << "  typical values from our current experience for best performance  " << Endl;
   Log() << "  are between 0.5(%) and 10(%) " << Endl;
   Log() << Endl;
   Log() << "The default minimal number is currently set to " << Endl;
   Log() << "   max(20, (N_training_events / N_variables^2 / 10)) " << Endl;
   Log() << "and can be changed by the user." << Endl;
   Log() << Endl;
   Log() << "The other crucial parameter, the pruning strength (\"PruneStrength\")," << Endl;
   Log() << "is also related to overtraining. It is a regularisation parameter " << Endl;
   Log() << "that is used when determining after the training which splits " << Endl;
   Log() << "are considered statistically insignificant and are removed. The" << Endl;
   Log() << "user is advised to carefully watch the BDT screen output for" << Endl;
   Log() << "the comparison between efficiencies obtained on the training and" << Endl;
   Log() << "the independent test sample. They should be equal within statistical" << Endl;
   Log() << "errors, in order to minimize statistical fluctuations in different samples." << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Make ROOT-independent C++ class for classifier response (classifier-specific implementation).

void TMVA::MethodBDT::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   TString nodeName = className;
   nodeName.ReplaceAll("Read","");
   nodeName.Append("Node");
   // write BDT-specific classifier response
   fout << "   std::vector<"<<nodeName<<"*> fForest;       // i.e. root nodes of decision trees" << std::endl;
   fout << "   std::vector<double>                fBoostWeights; // the weights applied in the individual boosts" << std::endl;
   fout << "};" << std::endl << std::endl;

   if(GetAnalysisType() == Types::kMulticlass) {
      fout << "std::vector<double> ReadBDTG::GetMulticlassValues__( const std::vector<double>& inputValues ) const" << std::endl;
      fout << "{" << std::endl;
      fout << "   uint nClasses = " << DataInfo().GetNClasses() << ";" << std::endl;
      fout << "   std::vector<double> fMulticlassReturnVal;" << std::endl;
      fout << "   fMulticlassReturnVal.reserve(nClasses);" << std::endl;
      fout << std::endl;
      fout << "   std::vector<double> temp(nClasses);" << std::endl;
      fout << "   auto forestSize = fForest.size();" << std::endl;
      fout << "   // trees 0, nClasses, 2*nClasses, ... belong to class 0" << std::endl;
      fout << "   // trees 1, nClasses+1, 2*nClasses+1, ... belong to class 1 and so forth" << std::endl;
      fout << "   uint classOfTree = 0;" << std::endl;
      fout << "   for (uint itree = 0; itree < forestSize; ++itree) {" << std::endl;
      fout << "      BDTGNode *current = fForest[itree];" << std::endl;
      fout << "      while (current->GetNodeType() == 0) { //intermediate node" << std::endl;
      fout << "         if (current->GoesRight(inputValues)) current=(BDTGNode*)current->GetRight();" << std::endl;
      fout << "         else current=(BDTGNode*)current->GetLeft();" << std::endl;
      fout << "      }" << std::endl;
      fout << "      temp[classOfTree] += current->GetResponse();" << std::endl;
      fout << "      if (++classOfTree == nClasses) classOfTree = 0; // cheap modulo" << std::endl;
      fout << "   }" << std::endl;
      fout << std::endl;
      fout << "   // we want to calculate sum of exp(temp[j] - temp[i]) for all i,j (i!=j)" << std::endl;
      fout << "   // first calculate exp(), then replace minus with division." << std::endl;
      fout << "   std::transform(temp.begin(), temp.end(), temp.begin(), [](double d){return exp(d);});" << std::endl;
      fout << std::endl;
      fout << "   for(uint iClass=0; iClass<nClasses; iClass++){" << std::endl;
      fout << "      double norm = 0.0;" << std::endl;
      fout << "      for(uint j=0;j<nClasses;j++){" << std::endl;
      fout << "         if(iClass!=j)" << std::endl;
      fout << "            norm += temp[j] / temp[iClass];" << std::endl;
      fout << "      }" << std::endl;
      fout << "      fMulticlassReturnVal.push_back(1.0/(1.0+norm));" << std::endl;
      fout << "   }" << std::endl;
      fout << std::endl;
      fout << "   return fMulticlassReturnVal;" << std::endl;
      fout << "}" << std::endl;
   } else {
      fout << "double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
      fout << "{" << std::endl;
      fout << "   double myMVA = 0;" << std::endl;
      if (fDoPreselection){
         for (UInt_t ivar = 0; ivar< fIsLowBkgCut.size(); ivar++){
            if (fIsLowBkgCut[ivar]){
               fout << "   if (inputValues["<<ivar<<"] < " << fLowBkgCut[ivar] << ") return -1;  // is background preselection cut" << std::endl;
            }
            if (fIsLowSigCut[ivar]){
               fout << "   if (inputValues["<<ivar<<"] < "<< fLowSigCut[ivar] << ") return  1;  // is signal preselection cut" << std::endl;
            }
            if (fIsHighBkgCut[ivar]){
               fout << "   if (inputValues["<<ivar<<"] > "<<fHighBkgCut[ivar] <<")  return -1;  // is background preselection cut" << std::endl;
            }
            if (fIsHighSigCut[ivar]){
               fout << "   if (inputValues["<<ivar<<"] > "<<fHighSigCut[ivar]<<")  return  1;  // is signal preselection cut" << std::endl;
            }
         }
      }

      if (fBoostType!="Grad"){
         fout << "   double norm  = 0;" << std::endl;
      }
      fout << "   for (unsigned int itree=0; itree<fForest.size(); itree++){" << std::endl;
      fout << "      "<<nodeName<<" *current = fForest[itree];" << std::endl;
      fout << "      while (current->GetNodeType() == 0) { //intermediate node" << std::endl;
      fout << "         if (current->GoesRight(inputValues)) current=("<<nodeName<<"*)current->GetRight();" << std::endl;
      fout << "         else current=("<<nodeName<<"*)current->GetLeft();" << std::endl;
      fout << "      }" << std::endl;
      if (fBoostType=="Grad"){
         fout << "      myMVA += current->GetResponse();" << std::endl;
      }else{
         if (fUseYesNoLeaf) fout << "      myMVA += fBoostWeights[itree] *  current->GetNodeType();" << std::endl;
         else               fout << "      myMVA += fBoostWeights[itree] *  current->GetPurity();" << std::endl;
         fout << "      norm  += fBoostWeights[itree];" << std::endl;
      }
      fout << "   }" << std::endl;
      if (fBoostType=="Grad"){
         fout << "   return 2.0/(1.0+exp(-2.0*myMVA))-1.0;" << std::endl;
      }
      else fout << "   return myMVA /= norm;" << std::endl;
      fout << "}" << std::endl << std::endl;
   }

   fout << "void " << className << "::Initialize()" << std::endl;
   fout << "{" << std::endl;
   fout << "  double inf = std::numeric_limits<double>::infinity();" << std::endl;
   fout << "  double nan = std::numeric_limits<double>::quiet_NaN();" << std::endl;
   //Now for each decision tree, write directly the constructors of the nodes in the tree structure
   for (UInt_t itree=0; itree<GetNTrees(); itree++) {
      fout << "  // itree = " << itree << std::endl;
      fout << "  fBoostWeights.push_back(" << fBoostWeights[itree] << ");" << std::endl;
      fout << "  fForest.push_back( " << std::endl;
      this->MakeClassInstantiateNode((DecisionTreeNode*)fForest[itree]->GetRoot(), fout, className);
      fout <<"   );" << std::endl;
   }
   fout << "   return;" << std::endl;
   fout << "};" << std::endl;
   fout << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   for (unsigned int itree=0; itree<fForest.size(); itree++) { " << std::endl;
   fout << "      delete fForest[itree]; " << std::endl;
   fout << "   }" << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Specific class header.

void TMVA::MethodBDT::MakeClassSpecificHeader(  std::ostream& fout, const TString& className) const
{
   TString nodeName = className;
   nodeName.ReplaceAll("Read","");
   nodeName.Append("Node");
   fout << "#include <algorithm>" << std::endl;
   fout << "#include <limits>" << std::endl;
   fout << std::endl;
   //fout << "#ifndef NN" << std::endl; commented out on purpose see next line
   fout << "#define NN new "<<nodeName << std::endl; // NN definition depends on individual methods. Important to have NO #ifndef if several BDT methods compile together
   //fout << "#endif" << std::endl; commented out on purpose see previous line
   fout << std::endl;
   fout << "#ifndef "<<nodeName<<"__def" << std::endl;
   fout << "#define "<<nodeName<<"__def" << std::endl;
   fout << std::endl;
   fout << "class "<<nodeName<<" {" << std::endl;
   fout << std::endl;
   fout << "public:" << std::endl;
   fout << std::endl;
   fout << "   // constructor of an essentially \"empty\" node floating in space" << std::endl;
   fout << "   "<<nodeName<<" ( "<<nodeName<<"* left,"<<nodeName<<"* right," << std::endl;
   if (fUseFisherCuts){
      fout << "                          int nFisherCoeff," << std::endl;
      for (UInt_t i=0;i<GetNVariables()+1;i++){
         fout << "                          double fisherCoeff"<<i<<"," << std::endl;
      }
   }
   fout << "                          int selector, double cutValue, bool cutType, " << std::endl;
   fout << "                          int nodeType, double purity, double response ) :" << std::endl;
   fout << "   fLeft         ( left         )," << std::endl;
   fout << "   fRight        ( right        )," << std::endl;
   if (fUseFisherCuts) fout << "   fNFisherCoeff ( nFisherCoeff )," << std::endl;
   fout << "   fSelector     ( selector     )," << std::endl;
   fout << "   fCutValue     ( cutValue     )," << std::endl;
   fout << "   fCutType      ( cutType      )," << std::endl;
   fout << "   fNodeType     ( nodeType     )," << std::endl;
   fout << "   fPurity       ( purity       )," << std::endl;
   fout << "   fResponse     ( response     ){" << std::endl;
   if (fUseFisherCuts){
      for (UInt_t i=0;i<GetNVariables()+1;i++){
         fout << "     fFisherCoeff.push_back(fisherCoeff"<<i<<");" << std::endl;
      }
   }
   fout << "   }" << std::endl << std::endl;
   fout << "   virtual ~"<<nodeName<<"();" << std::endl << std::endl;
   fout << "   // test event if it descends the tree at this node to the right" << std::endl;
   fout << "   virtual bool GoesRight( const std::vector<double>& inputValues ) const;" << std::endl;
   fout << "   "<<nodeName<<"* GetRight( void )  {return fRight; };" << std::endl << std::endl;
   fout << "   // test event if it descends the tree at this node to the left " << std::endl;
   fout << "   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;" << std::endl;
   fout << "   "<<nodeName<<"* GetLeft( void ) { return fLeft; };   " << std::endl << std::endl;
   fout << "   // return  S/(S+B) (purity) at this node (from  training)" << std::endl << std::endl;
   fout << "   double GetPurity( void ) const { return fPurity; } " << std::endl;
   fout << "   // return the node type" << std::endl;
   fout << "   int    GetNodeType( void ) const { return fNodeType; }" << std::endl;
   fout << "   double GetResponse(void) const {return fResponse;}" << std::endl << std::endl;
   fout << "private:" << std::endl << std::endl;
   fout << "   "<<nodeName<<"*   fLeft;     // pointer to the left daughter node" << std::endl;
   fout << "   "<<nodeName<<"*   fRight;    // pointer to the right daughter node" << std::endl;
   if (fUseFisherCuts){
      fout << "   int                     fNFisherCoeff; // =0 if this node doesn't use fisher, else =nvar+1 " << std::endl;
      fout << "   std::vector<double>     fFisherCoeff;  // the fisher coeff (offset at the last element)" << std::endl;
   }
   fout << "   int                     fSelector; // index of variable used in node selection (decision tree)   " << std::endl;
   fout << "   double                  fCutValue; // cut value applied on this node to discriminate bkg against sig" << std::endl;
   fout << "   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise" << std::endl;
   fout << "   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal " << std::endl;
   fout << "   double                  fPurity;   // Purity of node from training"<< std::endl;
   fout << "   double                  fResponse; // Regression response value of node" << std::endl;
   fout << "}; " << std::endl;
   fout << std::endl;
   fout << "//_______________________________________________________________________" << std::endl;
   fout << "   "<<nodeName<<"::~"<<nodeName<<"()" << std::endl;
   fout << "{" << std::endl;
   fout << "   if (fLeft  != NULL) delete fLeft;" << std::endl;
   fout << "   if (fRight != NULL) delete fRight;" << std::endl;
   fout << "}; " << std::endl;
   fout << std::endl;
   fout << "//_______________________________________________________________________" << std::endl;
   fout << "bool "<<nodeName<<"::GoesRight( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   // test event if it descends the tree at this node to the right" << std::endl;
   fout << "   bool result;" << std::endl;
   if (fUseFisherCuts){
      fout << "   if (fNFisherCoeff == 0){" << std::endl;
      fout << "     result = (inputValues[fSelector] >= fCutValue );" << std::endl;
      fout << "   }else{" << std::endl;
      fout << "     double fisher = fFisherCoeff.at(fFisherCoeff.size()-1);" << std::endl;
      fout << "     for (unsigned int ivar=0; ivar<fFisherCoeff.size()-1; ivar++)" << std::endl;
      fout << "       fisher += fFisherCoeff.at(ivar)*inputValues.at(ivar);" << std::endl;
      fout << "     result = fisher > fCutValue;" << std::endl;
      fout << "   }" << std::endl;
   }else{
      fout << "     result = (inputValues[fSelector] >= fCutValue );" << std::endl;
   }
   fout << "   if (fCutType == true) return result; //the cuts are selecting Signal ;" << std::endl;
   fout << "   else return !result;" << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "//_______________________________________________________________________" << std::endl;
   fout << "bool "<<nodeName<<"::GoesLeft( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   // test event if it descends the tree at this node to the left" << std::endl;
   fout << "   if (!this->GoesRight(inputValues)) return true;" << std::endl;
   fout << "   else return false;" << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "#endif" << std::endl;
   fout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively descends a tree and writes the node instance to the output stream.

void TMVA::MethodBDT::MakeClassInstantiateNode( DecisionTreeNode *n, std::ostream& fout, const TString& className ) const
{
   if (n == NULL) {
      Log() << kFATAL << "MakeClassInstantiateNode: started with undefined node" <<Endl;
      return ;
   }
   fout << "NN("<<std::endl;
   if (n->GetLeft() != NULL){
      this->MakeClassInstantiateNode( (DecisionTreeNode*)n->GetLeft() , fout, className);
   }
   else {
      fout << "0";
   }
   fout << ", " <<std::endl;
   if (n->GetRight() != NULL){
      this->MakeClassInstantiateNode( (DecisionTreeNode*)n->GetRight(), fout, className );
   }
   else {
      fout << "0";
   }
   fout << ", " <<  std::endl
        << std::setprecision(6);
   if (fUseFisherCuts){
      fout << n->GetNFisherCoeff() << ", ";
      for (UInt_t i=0; i< GetNVariables()+1; i++) {
         if (n->GetNFisherCoeff() == 0 ){
            fout <<  "0, ";
         }else{
            fout << n->GetFisherCoeff(i) << ", ";
         }
      }
   }
   fout << n->GetSelector() << ", "
        << n->GetCutValue() << ", "
        << n->GetCutType() << ", "
        << n->GetNodeType() << ", "
        << n->GetPurity() << ","
        << n->GetResponse() << ") ";
}

////////////////////////////////////////////////////////////////////////////////
/// Find useful preselection cuts that will be applied before
/// and Decision Tree training.. (and of course also applied
/// in the GetMVA .. --> -1 for background +1 for Signal)

void TMVA::MethodBDT::DeterminePreselectionCuts(const std::vector<const TMVA::Event*>& eventSample)
{
   Double_t nTotS = 0.0, nTotB = 0.0;
   Int_t nTotS_unWeighted = 0, nTotB_unWeighted = 0;

   std::vector<TMVA::BDTEventWrapper> bdtEventSample;

   fIsLowSigCut.assign(GetNvar(),kFALSE);
   fIsLowBkgCut.assign(GetNvar(),kFALSE);
   fIsHighSigCut.assign(GetNvar(),kFALSE);
   fIsHighBkgCut.assign(GetNvar(),kFALSE);

   fLowSigCut.assign(GetNvar(),0.);   //  ---------------| -->  in var is signal (accept all above lower cut)
   fLowBkgCut.assign(GetNvar(),0.);   //  ---------------| -->  in var is bkg    (accept all above lower cut)
   fHighSigCut.assign(GetNvar(),0.);  //  <-- | --------------  in var is signal (accept all blow cut)
   fHighBkgCut.assign(GetNvar(),0.);  //  <-- | --------------  in var is blg    (accept all blow cut)


   // Initialize (un)weighted counters for signal & background
   // Construct a list of event wrappers that point to the original data
   for( std::vector<const TMVA::Event*>::const_iterator it = eventSample.begin(); it != eventSample.end(); ++it ) {
      if (DataInfo().IsSignal(*it)){
         nTotS += (*it)->GetWeight();
         ++nTotS_unWeighted;
      }
      else {
         nTotB += (*it)->GetWeight();
         ++nTotB_unWeighted;
      }
      bdtEventSample.push_back(TMVA::BDTEventWrapper(*it));
   }

   for( UInt_t ivar = 0; ivar < GetNvar(); ivar++ ) { // loop over all discriminating variables
      TMVA::BDTEventWrapper::SetVarIndex(ivar); // select the variable to sort by
      std::sort( bdtEventSample.begin(),bdtEventSample.end() ); // sort the event data

      Double_t bkgWeightCtr = 0.0, sigWeightCtr = 0.0;
      std::vector<TMVA::BDTEventWrapper>::iterator it = bdtEventSample.begin(), it_end = bdtEventSample.end();
      for( ; it != it_end; ++it ) {
         if (DataInfo().IsSignal(**it))
            sigWeightCtr += (**it)->GetWeight();
         else
            bkgWeightCtr += (**it)->GetWeight();
         // Store the accumulated signal (background) weights
         it->SetCumulativeWeight(false,bkgWeightCtr);
         it->SetCumulativeWeight(true,sigWeightCtr);
      }

      //variable that determines how "exact" you cut on the preselection found in the training data. Here I chose
      //1% of the variable range...
      Double_t dVal = (DataInfo().GetVariableInfo(ivar).GetMax() - DataInfo().GetVariableInfo(ivar).GetMin())/100. ;
      Double_t nSelS, nSelB, effS=0.05, effB=0.05, rejS=0.05, rejB=0.05;
      Double_t tmpEffS, tmpEffB, tmpRejS, tmpRejB;
      // Locate the optimal cut for this (ivar-th) variable



      for(UInt_t iev = 1; iev < bdtEventSample.size(); iev++) {
         //dVal = bdtEventSample[iev].GetVal() - bdtEventSample[iev-1].GetVal();

         nSelS = bdtEventSample[iev].GetCumulativeWeight(true);
         nSelB = bdtEventSample[iev].GetCumulativeWeight(false);
         // you look for some 100% efficient pre-selection cut to remove background.. i.e. nSelS=0 && nSelB>5%nTotB or ( nSelB=0 nSelS>5%nTotS)
         tmpEffS=nSelS/nTotS;
         tmpEffB=nSelB/nTotB;
         tmpRejS=1-tmpEffS;
         tmpRejB=1-tmpEffB;
         if      (nSelS==0     && tmpEffB>effB)  {effB=tmpEffB; fLowBkgCut[ivar]  = bdtEventSample[iev].GetVal() - dVal; fIsLowBkgCut[ivar]=kTRUE;}
         else if (nSelB==0     && tmpEffS>effS)  {effS=tmpEffS; fLowSigCut[ivar]  = bdtEventSample[iev].GetVal() - dVal; fIsLowSigCut[ivar]=kTRUE;}
         else if (nSelB==nTotB && tmpRejS>rejS)  {rejS=tmpRejS; fHighSigCut[ivar] = bdtEventSample[iev].GetVal() + dVal; fIsHighSigCut[ivar]=kTRUE;}
         else if (nSelS==nTotS && tmpRejB>rejB)  {rejB=tmpRejB; fHighBkgCut[ivar] = bdtEventSample[iev].GetVal() + dVal; fIsHighBkgCut[ivar]=kTRUE;}

      }
   }

   Log() << kDEBUG << " \tfound and suggest the following possible pre-selection cuts " << Endl;
   if (fDoPreselection) Log() << kDEBUG << "\tthe training will be done after these cuts... and GetMVA value returns +1, (-1) for a signal (bkg) event that passes these cuts" << Endl;
   else  Log() << kDEBUG << "\tas option DoPreselection was not used, these cuts however will not be performed, but the training will see the full sample"<<Endl;
   for (UInt_t ivar=0; ivar < GetNvar(); ivar++ ) { // loop over all discriminating variables
      if (fIsLowBkgCut[ivar]){
         Log() << kDEBUG  << " \tfound cut: Bkg if var " << ivar << " < "  << fLowBkgCut[ivar] << Endl;
      }
      if (fIsLowSigCut[ivar]){
         Log() << kDEBUG  << " \tfound cut: Sig if var " << ivar << " < "  << fLowSigCut[ivar] << Endl;
      }
      if (fIsHighBkgCut[ivar]){
         Log() << kDEBUG  << " \tfound cut: Bkg if var " << ivar << " > "  << fHighBkgCut[ivar] << Endl;
      }
      if (fIsHighSigCut[ivar]){
         Log() << kDEBUG  << " \tfound cut: Sig if var " << ivar << " > "  << fHighSigCut[ivar] << Endl;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply the  preselection cuts before even bothering about any
/// Decision Trees  in the GetMVA .. --> -1 for background +1 for Signal

Double_t TMVA::MethodBDT::ApplyPreselectionCuts(const Event* ev)
{
   Double_t result=0;

   for (UInt_t ivar=0; ivar < GetNvar(); ivar++ ) { // loop over all discriminating variables
      if (fIsLowBkgCut[ivar]){
         if (ev->GetValue(ivar) < fLowBkgCut[ivar]) result = -1;  // is background
      }
      if (fIsLowSigCut[ivar]){
         if (ev->GetValue(ivar) < fLowSigCut[ivar]) result =  1;  // is signal
      }
      if (fIsHighBkgCut[ivar]){
         if (ev->GetValue(ivar) > fHighBkgCut[ivar]) result = -1;  // is background
      }
      if (fIsHighSigCut[ivar]){
         if (ev->GetValue(ivar) > fHighSigCut[ivar]) result =  1;  // is signal
      }
   }

   return result;
}

