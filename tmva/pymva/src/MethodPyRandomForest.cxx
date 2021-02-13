// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyRandomForest                                                  *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      Random Forest Classifiear from Scikit learn                               *
 *                                                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/
#include <Python.h>    // Needs to be included first to avoid redefinition of _POSIX_C_SOURCE
#include "TMVA/MethodPyRandomForest.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "TMVA/Configurable.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include "TMVA/IMethod.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDF.h"
#include "TMVA/Ranking.h"
#include "TMVA/Results.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Timer.h"
#include "TMVA/VariableTransformBase.h"

#include "TMatrix.h"

using namespace TMVA;

namespace TMVA {
namespace Internal {
class PyGILRAII {
   PyGILState_STATE m_GILState;

public:
   PyGILRAII() : m_GILState(PyGILState_Ensure()) {}
   ~PyGILRAII() { PyGILState_Release(m_GILState); }
};
} // namespace Internal
} // namespace TMVA

REGISTER_METHOD(PyRandomForest)

ClassImp(MethodPyRandomForest);

//_______________________________________________________________________
MethodPyRandomForest::MethodPyRandomForest(const TString &jobName,
      const TString &methodTitle,
      DataSetInfo &dsi,
      const TString &theOption) :
   PyMethodBase(jobName, Types::kPyRandomForest, methodTitle, dsi, theOption),
   fNestimators(10),
   fCriterion("gini"),
   fMaxDepth("None"),
   fMinSamplesSplit(2),
   fMinSamplesLeaf(1),
   fMinWeightFractionLeaf(0),
   fMaxFeatures("'auto'"),
   fMaxLeafNodes("None"),
   fBootstrap(kTRUE),
   fOobScore(kFALSE),
   fNjobs(1),
   fRandomState("None"),
   fVerbose(0),
   fWarmStart(kFALSE),
   fClassWeight("None")
{
}

//_______________________________________________________________________
MethodPyRandomForest::MethodPyRandomForest(DataSetInfo &theData, const TString &theWeightFile)
   : PyMethodBase(Types::kPyRandomForest, theData, theWeightFile),
   fNestimators(10),
   fCriterion("gini"),
   fMaxDepth("None"),
   fMinSamplesSplit(2),
   fMinSamplesLeaf(1),
   fMinWeightFractionLeaf(0),
   fMaxFeatures("'auto'"),
   fMaxLeafNodes("None"),
   fBootstrap(kTRUE),
   fOobScore(kFALSE),
   fNjobs(1),
   fRandomState("None"),
   fVerbose(0),
   fWarmStart(kFALSE),
   fClassWeight("None")
{
}


//_______________________________________________________________________
MethodPyRandomForest::~MethodPyRandomForest(void)
{
}

//_______________________________________________________________________
Bool_t MethodPyRandomForest::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass && numberClasses >= 2) return kTRUE;
   return kFALSE;
}

//_______________________________________________________________________
void MethodPyRandomForest::DeclareOptions()
{
   MethodBase::DeclareCompatibilityOptions();

   DeclareOptionRef(fNestimators, "NEstimators", "Integer, optional (default=10). The number of trees in the forest.");
   DeclareOptionRef(fCriterion, "Criterion", "String, optional (default='gini') \
      The function to measure the quality of a split. Supported criteria are \
      'gini' for the Gini impurity and 'entropy' for the information gain. \
      Note: this parameter is tree-specific.");

   DeclareOptionRef(fMaxDepth, "MaxDepth", "integer or None, optional (default=None) \
      The maximum depth of the tree. If None, then nodes are expanded until \
      all leaves are pure or until all leaves contain less than \
      min_samples_split samples. \
      Ignored if ``max_leaf_nodes`` is not None.");

   DeclareOptionRef(fMinSamplesSplit, "MinSamplesSplit", "integer, optional (default=2)\
      The minimum number of samples required to split an internal node.");

   DeclareOptionRef(fMinSamplesLeaf, "MinSamplesLeaf", "integer, optional (default=1) \
      The minimum number of samples in newly created leaves.  A split is \
      discarded if after the split, one of the leaves would contain less then \
      ``min_samples_leaf`` samples.");
   DeclareOptionRef(fMinWeightFractionLeaf, "MinWeightFractionLeaf", "//float, optional (default=0.) \
      The minimum weighted fraction of the input samples required to be at a \
      leaf node.");
   DeclareOptionRef(fMaxFeatures, "MaxFeatures", "The number of features to consider when looking for the best split");

   DeclareOptionRef(fMaxLeafNodes, "MaxLeafNodes", "int or None, optional (default=None)\
      Grow trees with ``max_leaf_nodes`` in best-first fashion.\
      Best nodes are defined as relative reduction in impurity.\
      If None then unlimited number of leaf nodes.\
      If not None then ``max_depth`` will be ignored.");

   DeclareOptionRef(fBootstrap, "Bootstrap", "boolean, optional (default=True) \
      Whether bootstrap samples are used when building trees.");

   DeclareOptionRef(fOobScore, "OoBScore", " bool Whether to use out-of-bag samples to estimate\
      the generalization error.");

   DeclareOptionRef(fNjobs, "NJobs", " integer, optional (default=1) \
      The number of jobs to run in parallel for both `fit` and `predict`. \
      If -1, then the number of jobs is set to the number of cores.");

   DeclareOptionRef(fRandomState, "RandomState", "int, RandomState instance or None, optional (default=None)\
      If int, random_state is the seed used by the random number generator;\
      If RandomState instance, random_state is the random number generator;\
      If None, the random number generator is the RandomState instance used\
      by `np.random`.");

   DeclareOptionRef(fVerbose, "Verbose", "int, optional (default=0)\
      Controls the verbosity of the tree building process.");

   DeclareOptionRef(fWarmStart, "WarmStart", "bool, optional (default=False)\
      When set to ``True``, reuse the solution of the previous call to fit\
      and add more estimators to the ensemble, otherwise, just fit a whole\
      new forest.");

   DeclareOptionRef(fClassWeight, "ClassWeight", "dict, list of dicts, \"auto\", \"subsample\" or None, optional\
      Weights associated with classes in the form ``{class_label: weight}``.\
      If not given, all classes are supposed to have weight one. For\
      multi-output problems, a list of dicts can be provided in the same\
      order as the columns of y.\
      The \"auto\" mode uses the values of y to automatically adjust\
      weights inversely proportional to class frequencies in the input data.\
      The \"subsample\" mode is the same as \"auto\" except that weights are\
      computed based on the bootstrap sample for every tree grown.\
      For multi-output, the weights of each column of y will be multiplied.\
      Note that these weights will be multiplied with sample_weight (passed\
      through the fit method) if sample_weight is specified.");

   DeclareOptionRef(fFilenameClassifier, "FilenameClassifier",
      "Store trained classifier in this file");
}

//_______________________________________________________________________
// Check options and load them to local python namespace
void MethodPyRandomForest::ProcessOptions()
{
   if (fNestimators <= 0) {
      Log() << kFATAL << " NEstimators <=0... that does not work !! " << Endl;
   }
   pNestimators = Eval(Form("%i", fNestimators));
   PyDict_SetItemString(fLocalNS, "nEstimators", pNestimators);

   if (fCriterion != "gini" && fCriterion != "entropy") {
      Log() << kFATAL << Form(" Criterion = %s... that does not work !! ", fCriterion.Data())
            << " The options are `gini` or `entropy`." << Endl;
   }
   pCriterion = Eval(Form("'%s'", fCriterion.Data()));
   PyDict_SetItemString(fLocalNS, "criterion", pCriterion);

   pMaxDepth = Eval(fMaxDepth);
   PyDict_SetItemString(fLocalNS, "maxDepth", pMaxDepth);
   if (!pMaxDepth) {
      Log() << kFATAL << Form(" MaxDepth = %s... that does not work !! ", fMaxDepth.Data())
            << " The options are None or integer." << Endl;
   }

   if (fMinSamplesSplit < 0) {
      Log() << kFATAL << " MinSamplesSplit < 0... that does not work !! " << Endl;
   }
   pMinSamplesSplit = Eval(Form("%i", fMinSamplesSplit));
   PyDict_SetItemString(fLocalNS, "minSamplesSplit", pMinSamplesSplit);

   if (fMinSamplesLeaf < 0) {
      Log() << kFATAL << " MinSamplesLeaf < 0... that does not work !! " << Endl;
   }
   pMinSamplesLeaf = Eval(Form("%i", fMinSamplesLeaf));
   PyDict_SetItemString(fLocalNS, "minSamplesLeaf", pMinSamplesLeaf);

   if (fMinWeightFractionLeaf < 0) {
      Log() << kERROR << " MinWeightFractionLeaf < 0... that does not work !! " << Endl;
   }
   pMinWeightFractionLeaf = Eval(Form("%f", fMinWeightFractionLeaf));
   PyDict_SetItemString(fLocalNS, "minWeightFractionLeaf", pMinWeightFractionLeaf);

   if (fMaxFeatures == "auto" || fMaxFeatures == "sqrt" || fMaxFeatures == "log2"){
      fMaxFeatures = Form("'%s'", fMaxFeatures.Data());
   }
   pMaxFeatures = Eval(fMaxFeatures);
   PyDict_SetItemString(fLocalNS, "maxFeatures", pMaxFeatures);

   if (!pMaxFeatures) {
      Log() << kFATAL << Form(" MaxFeatures = %s... that does not work !! ", fMaxFeatures.Data())
            << "int, float, string or None, optional (default='auto')"
            << "The number of features to consider when looking for the best split:"
            << "If int, then consider `max_features` features at each split."
            << "If float, then `max_features` is a percentage and"
            << "`int(max_features * n_features)` features are considered at each split."
            << "If 'auto', then `max_features=sqrt(n_features)`."
            << "If 'sqrt', then `max_features=sqrt(n_features)`."
            << "If 'log2', then `max_features=log2(n_features)`."
            << "If None, then `max_features=n_features`." << Endl;
   }

   pMaxLeafNodes = Eval(fMaxLeafNodes);
   if (!pMaxLeafNodes) {
      Log() << kFATAL << Form(" MaxLeafNodes = %s... that does not work !! ", fMaxLeafNodes.Data())
            << " The options are None or integer." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "maxLeafNodes", pMaxLeafNodes);

   pRandomState = Eval(fRandomState);
   if (!pRandomState) {
      Log() << kFATAL << Form(" RandomState = %s... that does not work !! ", fRandomState.Data())
            << "If int, random_state is the seed used by the random number generator;"
            << "If RandomState instance, random_state is the random number generator;"
            << "If None, the random number generator is the RandomState instance used by `np.random`." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "randomState", pRandomState);

   pClassWeight = Eval(fClassWeight);
   if (!pClassWeight) {
      Log() << kFATAL << Form(" ClassWeight = %s... that does not work !! ", fClassWeight.Data())
            << "dict, list of dicts, 'auto', 'subsample' or None, optional" << Endl;
   }
   PyDict_SetItemString(fLocalNS, "classWeight", pClassWeight);

   if(fNjobs < 1) {
      Log() << kFATAL << Form(" NJobs = %i... that does not work !! ", fNjobs)
            << "Value has to be greater than zero." << Endl;
   }
   pNjobs = Eval(Form("%i", fNjobs));
   PyDict_SetItemString(fLocalNS, "nJobs", pNjobs);

   pBootstrap = Eval(Form("%i", UInt_t(fBootstrap)));
   PyDict_SetItemString(fLocalNS, "bootstrap", pBootstrap);
   pOobScore = Eval(Form("%i", UInt_t(fOobScore)));
   PyDict_SetItemString(fLocalNS, "oobScore", pOobScore);
   pVerbose = Eval(Form("%i", fVerbose));
   PyDict_SetItemString(fLocalNS, "verbose", pVerbose);
   pWarmStart = Eval(Form("%i", UInt_t(fWarmStart)));
   PyDict_SetItemString(fLocalNS, "warmStart", pWarmStart);

   // If no filename is given, set default
   if(fFilenameClassifier.IsNull())
   {
      fFilenameClassifier = GetWeightFileDir() + "/PyRFModel_" + GetName() + ".PyData";
   }
}

//_______________________________________________________________________
void MethodPyRandomForest::Init()
{
   TMVA::Internal::PyGILRAII raii;
   _import_array(); //require to use numpy arrays

   // Check options and load them to local python namespace
   ProcessOptions();

   // Import module for random forest classifier
   PyRunString("import sklearn.ensemble");

   // Get data properties
   fNvars = GetNVariables();
   fNoutputs = DataInfo().GetNClasses();
}

//_______________________________________________________________________
void MethodPyRandomForest::Train()
{
   // Load training data (data, classes, weights) to python arrays
   int fNrowsTraining = Data()->GetNTrainingEvents(); //every row is an event, a class type and a weight
   npy_intp dimsData[2];
   dimsData[0] = fNrowsTraining;
   dimsData[1] = fNvars;
   PyArrayObject * fTrainData = (PyArrayObject *)PyArray_SimpleNew(2, dimsData, NPY_FLOAT);
   PyDict_SetItemString(fLocalNS, "trainData", (PyObject*)fTrainData);
   float *TrainData = (float *)(PyArray_DATA(fTrainData));

   npy_intp dimsClasses = (npy_intp) fNrowsTraining;
   PyArrayObject * fTrainDataClasses = (PyArrayObject *)PyArray_SimpleNew(1, &dimsClasses, NPY_FLOAT);
   PyDict_SetItemString(fLocalNS, "trainDataClasses", (PyObject*)fTrainDataClasses);
   float *TrainDataClasses = (float *)(PyArray_DATA(fTrainDataClasses));

   PyArrayObject * fTrainDataWeights = (PyArrayObject *)PyArray_SimpleNew(1, &dimsClasses, NPY_FLOAT);
   PyDict_SetItemString(fLocalNS, "trainDataWeights", (PyObject*)fTrainDataWeights);
   float *TrainDataWeights = (float *)(PyArray_DATA(fTrainDataWeights));

   for (int i = 0; i < fNrowsTraining; i++) {
      // Fill training data matrix
      const TMVA::Event *e = Data()->GetTrainingEvent(i);
      for (UInt_t j = 0; j < fNvars; j++) {
         TrainData[j + i * fNvars] = e->GetValue(j);
      }

      // Fill target classes
      TrainDataClasses[i] = e->GetClass();

      // Get event weight
      TrainDataWeights[i] = e->GetWeight();
   }

   // Create classifier object
   PyRunString("classifier = sklearn.ensemble.RandomForestClassifier(bootstrap=bootstrap, class_weight=classWeight, criterion=criterion, max_depth=maxDepth, max_features=maxFeatures, max_leaf_nodes=maxLeafNodes, min_samples_leaf=minSamplesLeaf, min_samples_split=minSamplesSplit, min_weight_fraction_leaf=minWeightFractionLeaf, n_estimators=nEstimators, n_jobs=nJobs, oob_score=oobScore, random_state=randomState, verbose=verbose, warm_start=warmStart)",
      "Failed to setup classifier");

   // Fit classifier
   // NOTE: We dump the output to a variable so that the call does not pollute stdout
   PyRunString("dump = classifier.fit(trainData, trainDataClasses, trainDataWeights)", "Failed to train classifier");

   // Store classifier
   fClassifier = PyDict_GetItemString(fLocalNS, "classifier");
   if(fClassifier == 0) {
      Log() << kFATAL << "Can't create classifier object from RandomForestClassifier" << Endl;
      Log() << Endl;
   }

   if (IsModelPersistence()) {
      Log() << Endl;
      Log() << gTools().Color("bold") << "Saving state file: " << gTools().Color("reset") << fFilenameClassifier << Endl;
      Log() << Endl;
      Serialize(fFilenameClassifier, fClassifier);
   }
}

//_______________________________________________________________________
void MethodPyRandomForest::TestClassification()
{
   MethodBase::TestClassification();
}

//_______________________________________________________________________
std::vector<Double_t> MethodPyRandomForest::GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress)
{
   // Load model if not already done
   if (fClassifier == 0) ReadModelFromFile();

   // Determine number of events
   Long64_t nEvents = Data()->GetNEvents();
   if (firstEvt > lastEvt || lastEvt > nEvents) lastEvt = nEvents;
   if (firstEvt < 0) firstEvt = 0;
   nEvents = lastEvt-firstEvt;

     // use timer
   Timer timer( nEvents, GetName(), kTRUE );

   if (logProgress)
      Log() << kHEADER << Form("[%s] : ",DataInfo().GetName())
            << "Evaluation of " << GetMethodName() << " on "
            << (Data()->GetCurrentType() == Types::kTraining ? "training" : "testing")
            << " sample (" << nEvents << " events)" << Endl;

   // Get data
   npy_intp dims[2];
   dims[0] = nEvents;
   dims[1] = fNvars;
   PyArrayObject *pEvent= (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT);
   float *pValue = (float *)(PyArray_DATA(pEvent));

   for (Int_t ievt=0; ievt<nEvents; ievt++) {
      Data()->SetCurrentEvent(ievt);
      const TMVA::Event *e = Data()->GetEvent();
      for (UInt_t i = 0; i < fNvars; i++) {
         pValue[ievt * fNvars + i] = e->GetValue(i);
      }
   }

   // Get prediction from classifier
   PyArrayObject *result = (PyArrayObject *)PyObject_CallMethod(fClassifier, const_cast<char *>("predict_proba"), const_cast<char *>("(O)"), pEvent);
   double *proba = (double *)(PyArray_DATA(result));

   // Return signal probabilities
   if(Long64_t(mvaValues.size()) != nEvents) mvaValues.resize(nEvents);
   for (int i = 0; i < nEvents; ++i) {
      mvaValues[i] = proba[fNoutputs*i + TMVA::Types::kSignal];
   }

   Py_DECREF(pEvent);
   Py_DECREF(result);
   
   if (logProgress) {
      Log() << kINFO 
            << "Elapsed time for evaluation of " << nEvents <<  " events: "
            << timer.GetElapsedTime() << "       " << Endl;
   }

   return mvaValues;
}

//_______________________________________________________________________
Double_t MethodPyRandomForest::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
   // cannot determine error
   NoErrorCalc(errLower, errUpper);

   // Load model if not already done
   if (fClassifier == 0) ReadModelFromFile();

   // Get current event and load to python array
   const TMVA::Event *e = Data()->GetEvent();
   npy_intp dims[2];
   dims[0] = 1;
   dims[1] = fNvars;
   PyArrayObject *pEvent= (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT);
   float *pValue = (float *)(PyArray_DATA(pEvent));
   for (UInt_t i = 0; i < fNvars; i++) pValue[i] = e->GetValue(i);

   // Get prediction from classifier
   PyArrayObject *result = (PyArrayObject *)PyObject_CallMethod(fClassifier, const_cast<char *>("predict_proba"), const_cast<char *>("(O)"), pEvent);
   double *proba = (double *)(PyArray_DATA(result));

   // Return MVA value
   Double_t mvaValue;
   mvaValue = proba[TMVA::Types::kSignal]; // getting signal probability

   Py_DECREF(result);
   Py_DECREF(pEvent);

   return mvaValue;
}

//_______________________________________________________________________
std::vector<Float_t>& MethodPyRandomForest::GetMulticlassValues()
{
   // Load model if not already done
   if (fClassifier == 0) ReadModelFromFile();

   // Get current event and load to python array
   const TMVA::Event *e = Data()->GetEvent();
   npy_intp dims[2];
   dims[0] = 1;
   dims[1] = fNvars;
   PyArrayObject *pEvent= (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT);
   float *pValue = (float *)(PyArray_DATA(pEvent));
   for (UInt_t i = 0; i < fNvars; i++) pValue[i] = e->GetValue(i);

   // Get prediction from classifier
   PyArrayObject *result = (PyArrayObject *)PyObject_CallMethod(fClassifier, const_cast<char *>("predict_proba"), const_cast<char *>("(O)"), pEvent);
   double *proba = (double *)(PyArray_DATA(result));

   // Return MVA values
   if(UInt_t(classValues.size()) != fNoutputs) classValues.resize(fNoutputs);
   for(UInt_t i = 0; i < fNoutputs; i++) classValues[i] = proba[i];

   Py_DECREF(pEvent);
   Py_DECREF(result);

   return classValues;
}

//_______________________________________________________________________
void MethodPyRandomForest::ReadModelFromFile()
{
   if (!PyIsInitialized()) {
      PyInitialize();
   }

   Log() << Endl;
   Log() << gTools().Color("bold") << "Loading state file: " << gTools().Color("reset") << fFilenameClassifier << Endl;
   Log() << Endl;

   // Load classifier from file
   Int_t err = UnSerialize(fFilenameClassifier, &fClassifier);
   if(err != 0)
   {
       Log() << kFATAL << Form("Failed to load classifier from file (error code: %i): %s", err, fFilenameClassifier.Data()) << Endl;
   }

   // Book classifier object in python dict
   PyDict_SetItemString(fLocalNS, "classifier", fClassifier);

   // Load data properties
   // NOTE: This has to be repeated here for the reader application
   fNvars = GetNVariables();
   fNoutputs = DataInfo().GetNClasses();
}

//_______________________________________________________________________
const Ranking* MethodPyRandomForest::CreateRanking()
{
   // Get feature importance from classifier as an array with length equal
   // number of variables, higher value signals a higher importance
   PyArrayObject* pRanking = (PyArrayObject*) PyObject_GetAttrString(fClassifier, "feature_importances_");
   if(pRanking == 0) Log() << kFATAL << "Failed to get ranking from classifier" << Endl;

   // Fill ranking object and return it
   fRanking = new Ranking(GetName(), "Variable Importance");
   Double_t* rankingData = (Double_t*) PyArray_DATA(pRanking);
   for(UInt_t iVar=0; iVar<fNvars; iVar++){
      fRanking->AddRank(Rank(GetInputLabel(iVar), rankingData[iVar]));
   }

   Py_DECREF(pRanking);

   return fRanking;
}

//_______________________________________________________________________
void MethodPyRandomForest::GetHelpMessage() const
{
   // typical length of text line:
   //       "|--------------------------------------------------------------|"
   Log() << "A random forest is a meta estimator that fits a number of decision" << Endl;
   Log() << "tree classifiers on various sub-samples of the dataset and use" << Endl;
   Log() << "averaging to improve the predictive accuracy and control over-fitting." << Endl;
   Log() << Endl;
   Log() << "Check out the scikit-learn documentation for more information." << Endl;
}
