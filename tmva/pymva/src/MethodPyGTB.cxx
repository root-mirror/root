// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyGTB                                                           *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      GradientBoostingClassifier Classifiear from Scikit learn                  *
 *                                                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/

#include <Python.h>    // Needs to be included first to avoid redefinition of _POSIX_C_SOURCE
#include "TMVA/MethodPyGTB.h"

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
#include "TMVA/ResultsClassification.h"
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

REGISTER_METHOD(PyGTB)

ClassImp(MethodPyGTB);

//_______________________________________________________________________
MethodPyGTB::MethodPyGTB(const TString &jobName,
                         const TString &methodTitle,
                         DataSetInfo &dsi,
                         const TString &theOption) :
   PyMethodBase(jobName, Types::kPyGTB, methodTitle, dsi, theOption),
   fLoss("deviance"),
   fLearningRate(0.1),
   fNestimators(100),
   fSubsample(1.0),
   fMinSamplesSplit(2),
   fMinSamplesLeaf(1),
   fMinWeightFractionLeaf(0.0),
   fMaxDepth(3),
   fInit("None"),
   fRandomState("None"),
   fMaxFeatures("None"),
   fVerbose(0),
   fMaxLeafNodes("None"),
   fWarmStart(kFALSE)
{
}

//_______________________________________________________________________
MethodPyGTB::MethodPyGTB(DataSetInfo &theData, const TString &theWeightFile)
   : PyMethodBase(Types::kPyGTB, theData, theWeightFile),
   fLoss("deviance"),
   fLearningRate(0.1),
   fNestimators(100),
   fSubsample(1.0),
   fMinSamplesSplit(2),
   fMinSamplesLeaf(1),
   fMinWeightFractionLeaf(0.0),
   fMaxDepth(3),
   fInit("None"),
   fRandomState("None"),
   fMaxFeatures("None"),
   fVerbose(0),
   fMaxLeafNodes("None"),
   fWarmStart(kFALSE)
{
}


//_______________________________________________________________________
MethodPyGTB::~MethodPyGTB(void)
{
}

//_______________________________________________________________________
Bool_t MethodPyGTB::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass && numberClasses >= 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void MethodPyGTB::DeclareOptions()
{
   MethodBase::DeclareCompatibilityOptions();

   DeclareOptionRef(fLoss, "Loss", "{'deviance', 'exponential'}, optional (default='deviance')\
      loss function to be optimized. 'deviance' refers to\
      deviance (= logistic regression) for classification\
      with probabilistic outputs. For loss 'exponential' gradient\
      boosting recovers the AdaBoost algorithm.");

   DeclareOptionRef(fLearningRate, "LearningRate", "float, optional (default=0.1)\
      learning rate shrinks the contribution of each tree by `learning_rate`.\
      There is a trade-off between learning_rate and n_estimators.");

   DeclareOptionRef(fNestimators, "NEstimators", "int (default=100)\
      The number of boosting stages to perform. Gradient boosting\
      is fairly robust to over-fitting so a large number usually\
      results in better performance.");

   DeclareOptionRef(fSubsample, "Subsample", "float, optional (default=1.0)\
      The fraction of samples to be used for fitting the individual base\
      learners. If smaller than 1.0 this results in Stochastic Gradient\
      Boosting. `subsample` interacts with the parameter `n_estimators`.\
      Choosing `subsample < 1.0` leads to a reduction of variance\
      and an increase in bias.");

   DeclareOptionRef(fMinSamplesSplit, "MinSamplesSplit", "integer, optional (default=2)\
      The minimum number of samples required to split an internal node.");

   DeclareOptionRef(fMinSamplesLeaf, "MinSamplesLeaf", "integer, optional (default=1) \
      The minimum number of samples in newly created leaves.  A split is \
      discarded if after the split, one of the leaves would contain less then \
      ``min_samples_leaf`` samples.");

   DeclareOptionRef(fMinWeightFractionLeaf, "MinWeightFractionLeaf", "//float, optional (default=0.) \
      The minimum weighted fraction of the input samples required to be at a \
      leaf node.");

   DeclareOptionRef(fMaxDepth, "MaxDepth", "integer or None, optional (default=None) \
      The maximum depth of the tree. If None, then nodes are expanded until \
      all leaves are pure or until all leaves contain less than \
      min_samples_split samples. \
      Ignored if ``max_leaf_nodes`` is not None.");

   DeclareOptionRef(fInit, "Init", "BaseEstimator, None, optional (default=None)\
      An estimator object that is used to compute the initial\
      predictions. ``init`` has to provide ``fit`` and ``predict``.\
      If None it uses ``loss.init_estimator`");

   DeclareOptionRef(fRandomState, "RandomState", "int, RandomState instance or None, optional (default=None)\
      If int, random_state is the seed used by the random number generator;\
      If RandomState instance, random_state is the random number generator;\
      If None, the random number generator is the RandomState instance used\
      by `np.random`.");

   DeclareOptionRef(fMaxFeatures, "MaxFeatures", "The number of features to consider when looking for the best split");

   DeclareOptionRef(fVerbose, "Verbose", "int, optional (default=0)\
      Controls the verbosity of the tree building process.");

   DeclareOptionRef(fMaxLeafNodes, "MaxLeafNodes", "int or None, optional (default=None)\
      Grow trees with ``max_leaf_nodes`` in best-first fashion.\
      Best nodes are defined as relative reduction in impurity.\
      If None then unlimited number of leaf nodes.\
      If not None then ``max_depth`` will be ignored.");

   DeclareOptionRef(fWarmStart, "WarmStart", "bool, optional (default=False)\
      When set to ``True``, reuse the solution of the previous call to fit\
      and add more estimators to the ensemble, otherwise, just fit a whole\
      new forest.");

   DeclareOptionRef(fFilenameClassifier, "FilenameClassifier",
      "Store trained classifier in this file");
}

//_______________________________________________________________________
// Check options and load them to local python namespace
void MethodPyGTB::ProcessOptions()
{
   if (fLoss != "deviance" && fLoss != "exponential") {
      Log() << kFATAL << Form("Loss = %s ... that does not work!", fLoss.Data())
            << " The options are 'deviance' or 'exponential'." << Endl;
   }
   pLoss = Eval(Form("'%s'", fLoss.Data()));
   PyDict_SetItemString(fLocalNS, "loss", pLoss);

   if (fLearningRate <= 0) {
      Log() << kFATAL << "LearningRate <= 0 ... that does not work!" << Endl;
   }
   pLearningRate = Eval(Form("%f", fLearningRate));
   PyDict_SetItemString(fLocalNS, "learningRate", pLearningRate);

   if (fNestimators <= 0) {
      Log() << kFATAL << "NEstimators <= 0 ... that does not work!" << Endl;
   }
   pNestimators = Eval(Form("%i", fNestimators));
   PyDict_SetItemString(fLocalNS, "nEstimators", pNestimators);

   if (fMinSamplesSplit < 0) {
      Log() << kFATAL << "MinSamplesSplit < 0 ... that does not work!" << Endl;
   }
   pMinSamplesSplit = Eval(Form("%i", fMinSamplesSplit));
   PyDict_SetItemString(fLocalNS, "minSamplesSplit", pMinSamplesSplit);

   if (fSubsample < 0) {
      Log() << kFATAL << "Subsample < 0 ... that does not work!" << Endl;
   }
   pSubsample = Eval(Form("%f", fSubsample));
   PyDict_SetItemString(fLocalNS, "subsample", pSubsample);

   if (fMinSamplesLeaf < 0) {
      Log() << kFATAL << "MinSamplesLeaf < 0 ... that does not work!" << Endl;
   }
   pMinSamplesLeaf = Eval(Form("%i", fMinSamplesLeaf));
   PyDict_SetItemString(fLocalNS, "minSamplesLeaf", pMinSamplesLeaf);

   if (fMinSamplesSplit < 0) {
      Log() << kFATAL << "MinSamplesSplit < 0 ... that does not work!" << Endl;
   }
   pMinSamplesSplit = Eval(Form("%i", fMinSamplesSplit));
   PyDict_SetItemString(fLocalNS, "minSamplesSplit", pMinSamplesSplit);

   if (fMinWeightFractionLeaf < 0) {
      Log() << kFATAL << "MinWeightFractionLeaf < 0 ... that does not work !" << Endl;
   }
   pMinWeightFractionLeaf = Eval(Form("%f", fMinWeightFractionLeaf));
   PyDict_SetItemString(fLocalNS, "minWeightFractionLeaf", pMinWeightFractionLeaf);

   if (fMaxDepth <= 0) {
      Log() << kFATAL << " MaxDepth <= 0 ... that does not work !! " << Endl;
   }
   pMaxDepth = Eval(Form("%i", fMaxDepth));
   PyDict_SetItemString(fLocalNS, "maxDepth", pMaxDepth);

   pInit = Eval(fInit);
   if (!pInit) {
      Log() << kFATAL << Form("Init = %s ... that does not work!", fInit.Data())
            << " The options are None or BaseEstimator, which is an estimator object that"
            << "is used to compute the initial predictions. "
            << "'init' has to provide 'fit' and 'predict' methods."
            << " If None it uses 'loss.init_estimator'." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "init", pInit);

   pRandomState = Eval(fRandomState);
   if (!pRandomState) {
      Log() << kFATAL << Form(" RandomState = %s ... that does not work! ", fRandomState.Data())
            << " If int, random_state is the seed used by the random number generator;"
            << " If RandomState instance, random_state is the random number generator;"
            << " If None, the random number generator is the RandomState instance used by 'np.random'."
            << Endl;
   }
   PyDict_SetItemString(fLocalNS, "randomState", pRandomState);

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
      Log() << kFATAL << Form(" MaxLeafNodes = %s... that does not work!", fMaxLeafNodes.Data())
            << " The options are None or integer." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "maxLeafNodes", pMaxLeafNodes);

   pVerbose = Eval(Form("%i", fVerbose));
   PyDict_SetItemString(fLocalNS, "verbose", pVerbose);

   pWarmStart = Eval(Form("%i", UInt_t(fWarmStart)));
   PyDict_SetItemString(fLocalNS, "warmStart", pWarmStart);

   // If no filename is given, set default
   if(fFilenameClassifier.IsNull()) {
      fFilenameClassifier = GetWeightFileDir() + "/PyGTBModel_" + GetName() + ".PyData";
   }
}

//_______________________________________________________________________
void  MethodPyGTB::Init()
{
   TMVA::Internal::PyGILRAII raii;
   _import_array(); //require to use numpy arrays

   // Check options and load them to local python namespace
   ProcessOptions();

   // Import module for gradient tree boosting classifier
   PyRunString("import sklearn.ensemble");

   // Get data properties
   fNvars = GetNVariables();
   fNoutputs = DataInfo().GetNClasses();
}

void MethodPyGTB::Train()
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
   PyRunString("classifier = sklearn.ensemble.GradientBoostingClassifier(loss=loss, learning_rate=learningRate, n_estimators=nEstimators, max_depth=maxDepth, min_samples_split=minSamplesSplit, min_samples_leaf=minSamplesLeaf, min_weight_fraction_leaf=minWeightFractionLeaf, subsample=subsample, max_features=maxFeatures, max_leaf_nodes=maxLeafNodes, init=init, verbose=verbose, warm_start=warmStart, random_state=randomState)",
      "Failed to setup classifier");

   // Fit classifier
   // NOTE: We dump the output to a variable so that the call does not pollute stdout
   PyRunString("dump = classifier.fit(trainData, trainDataClasses, trainDataWeights)", "Failed to train classifier");

   // Store classifier
   fClassifier = PyDict_GetItemString(fLocalNS, "classifier");
   if(fClassifier == 0) {
      Log() << kFATAL << "Can't create classifier object from GradientBoostingClassifier" << Endl;
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
void MethodPyGTB::TestClassification()
{
   MethodBase::TestClassification();
}

//_______________________________________________________________________
std::vector<Double_t> MethodPyGTB::GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress)
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
Double_t MethodPyGTB::GetMvaValue(Double_t *errLower, Double_t *errUpper)
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
std::vector<Float_t>& MethodPyGTB::GetMulticlassValues()
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
void MethodPyGTB::ReadModelFromFile()
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
const Ranking* MethodPyGTB::CreateRanking()
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
void MethodPyGTB::GetHelpMessage() const
{
   // typical length of text line:
   //       "|--------------------------------------------------------------|"
   Log() << "A gradient tree boosting classifier builds a model from an ensemble" << Endl;
   Log() << "of decision trees, which are adapted each boosting step to fit better" << Endl;
   Log() << "to previously misclassified events." << Endl;
   Log() << Endl;
   Log() << "Check out the scikit-learn documentation for more information." << Endl;
}


