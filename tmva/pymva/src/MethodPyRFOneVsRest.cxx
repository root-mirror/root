// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyRFOneVsRest                                                  *
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
#include "TMVA/MethodPyRFOneVsRest.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"

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
#include "TMVA/VariableTransformBase.h"

#include "Riostream.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <iomanip>
#include <fstream>

using namespace TMVA;

REGISTER_METHOD(PyRFOneVsRest)

ClassImp(MethodPyRFOneVsRest)

//_______________________________________________________________________
MethodPyRFOneVsRest::MethodPyRFOneVsRest(const TString &jobName,
                                          const TString &methodTitle,
                                          DataSetInfo &dsi,
                                          const TString &theOption) :
  PyMethodBase(jobName, Types::kPyRFOneVsRest, methodTitle, dsi, theOption),
  n_estimators(10),
  criterion("gini"),
  max_depth("None"),
  min_samples_split(2),
  min_samples_leaf(1),
  min_weight_fraction_leaf(0),
  max_features("'auto'"),
  max_leaf_nodes("None"),
  bootstrap(kTRUE),
  oob_score(kFALSE),
  n_jobs(1),
  random_state("None"),
  verbose(0),
  warm_start(kFALSE),
  class_weight("None"),
  n_jobsOVR(1)    // for OneVsRest Classifier
{
}

//_______________________________________________________________________
MethodPyRFOneVsRest::MethodPyRFOneVsRest(DataSetInfo &theData, const TString &theWeightFile)
   :PyMethodBase(Types::kPyRFOneVsRest, theData, theWeightFile),
    n_estimators(10),
    criterion("gini"),
    max_depth("None"),
    min_samples_split(2),
    min_samples_leaf(1),
    min_weight_fraction_leaf(0),
    max_features("'auto'"),
    max_leaf_nodes("None"),
    bootstrap(kTRUE),
    oob_score(kFALSE),
    n_jobs(1),
    random_state("None"),
    verbose(0),
    warm_start(kFALSE),
    class_weight("None"),
    n_jobsOVR(1)    // for OneVsRest Classifier
{
}


//_______________________________________________________________________
MethodPyRFOneVsRest::~MethodPyRFOneVsRest(void)
{
}

//_______________________________________________________________________
Bool_t MethodPyRFOneVsRest::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void MethodPyRFOneVsRest::DeclareOptions()
{
   MethodBase::DeclareCompatibilityOptions();

   DeclareOptionRef(n_estimators, "NEstimators", "Integer, optional (default=10). The number of trees in the forest.");
   DeclareOptionRef(criterion, "Criterion", "//string, optional (default='gini') \
    The function to measure the quality of a split. Supported criteria are \
    'gini' for the Gini impurity and 'entropy' for the information gain. \
    Note: this parameter is tree-specific.");

   DeclareOptionRef(max_depth, "MaxDepth", "integer or None, optional (default=None) \
                                             The maximum depth of the tree. If None, then nodes are expanded until \
                                             all leaves are pure or until all leaves contain less than \
                                             min_samples_split samples. \
                                             Ignored if ``max_leaf_nodes`` is not None.");
   DeclareOptionRef(min_samples_split, "MinSamplesSplit", "integer, optional (default=2)\
    The minimum number of samples required to split an internal node.");

   DeclareOptionRef(min_samples_leaf, "MinSamplesLeaf", "integer, optional (default=1) \
    The minimum number of samples in newly created leaves.  A split is \
    discarded if after the split, one of the leaves would contain less then \
    ``min_samples_leaf`` samples.");
   DeclareOptionRef(min_weight_fraction_leaf, "MinWeightFractionLeaf", "//float, optional (default=0.) \
    The minimum weighted fraction of the input samples required to be at a \
    leaf node.");
   DeclareOptionRef(max_features, "MaxFeatures", "The number of features to consider when looking for the best split");
   DeclareOptionRef(max_leaf_nodes, "MaxLeafNodes", "int or None, optional (default=None)\
    Grow trees with ``max_leaf_nodes`` in best-first fashion.\
    Best nodes are defined as relative reduction in impurity.\
    If None then unlimited number of leaf nodes.\
    If not None then ``max_depth`` will be ignored.");
   DeclareOptionRef(bootstrap, "Bootstrap", "boolean, optional (default=True) \
    Whether bootstrap samples are used when building trees.");
   DeclareOptionRef(oob_score, "OoBScore", " bool Whether to use out-of-bag samples to estimate\
    the generalization error.");
   DeclareOptionRef(n_jobs, "NJobs", " integer, optional (default=1) \
    The number of jobs to run in parallel for both `fit` and `predict`. \
    If -1, then the number of jobs is set to the number of cores.");

   DeclareOptionRef(random_state, "RandomState", "int, RandomState instance or None, optional (default=None)\
    If int, random_state is the seed used by the random number generator;\
    If RandomState instance, random_state is the random number generator;\
    If None, the random number generator is the RandomState instance used\
    by `np.random`.");
   DeclareOptionRef(verbose, "Verbose", "int, optional (default=0)\
    Controls the verbosity of the tree building process.");
   DeclareOptionRef(warm_start, "WarmStart", "bool, optional (default=False)\
    When set to ``True``, reuse the solution of the previous call to fit\
    and add more estimators to the ensemble, otherwise, just fit a whole\
    new forest.");
   DeclareOptionRef(class_weight, "ClassWeight", "dict, list of dicts, \"auto\", \"subsample\" or None, optional\
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
   DeclareOptionRef(n_jobsOVR, "NJobsOVR", " integer, optional (default=1) \
    The number of jobs to use for the computation. If -1 all CPUs are used. \
    If 1 is given, no parallel computing code is used at all, which is \
    useful for debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. \
    Thus for n_jobs = -2, all CPUs but one are used.");
}

//_______________________________________________________________________
void MethodPyRFOneVsRest::ProcessOptions()
{
   if (n_estimators <= 0) {
      Log() << kERROR << " NEstimators <=0... that does not work !! "
            << " I set it to 10 .. just so that the program does not crash"
            << Endl;
      n_estimators = 10;
   }
   if (criterion != "gini" && criterion != "entropy") {
      Log() << kFATAL << Form(" Criterion = %s... that does not work !! ", criterion.Data())
            << " The options are gini of entropy."
            << Endl;
   }
   PyObject *pomax_depth = Eval(max_depth);
   if (!pomax_depth) {
      Log() << kFATAL << Form(" MaxDepth = %s... that does not work !! ", criterion.Data())
            << " The options are None or integer."
            << Endl;
   }
   Py_DECREF(pomax_depth);

   if (min_samples_split < 0) {
      Log() << kERROR << " MinSamplesSplit < 0... that does not work !! "
            << " I set it to 2 .. just so that the program does not crash"
            << Endl;
      min_samples_split = 2;
   }
   if (min_samples_leaf < 0) {
      Log() << kERROR << " MinSamplesLeaf < 0... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      min_samples_leaf = 1;
   }

   if (min_weight_fraction_leaf < 0) {
      Log() << kERROR << " MinWeightFractionLeaf < 0... that does not work !! "
            << " I set it to 0 .. just so that the program does not crash"
            << Endl;
      min_weight_fraction_leaf = 0;
   }
   if (max_features == "auto" || max_features == "sqrt" || max_features == "log2")max_features = Form("'%s'", max_features.Data());
   PyObject *pomax_features = Eval(max_features);
   if (!pomax_features) {
      Log() << kFATAL << Form(" MaxFeatures = %s... that does not work !! ", max_features.Data())
            << "int, float, string or None, optional (default='auto')"
            << "The number of features to consider when looking for the best split:"
            << "If int, then consider `max_features` features at each split."
            << "If float, then `max_features` is a percentage and"
            << "`int(max_features * n_features)` features are considered at each split."
            << "If 'auto', then `max_features=sqrt(n_features)`."
            << "If 'sqrt', then `max_features=sqrt(n_features)`."
            << "If 'log2', then `max_features=log2(n_features)`."
            << "If None, then `max_features=n_features`."
            << Endl;
   }
   Py_DECREF(pomax_features);

   PyObject *pomax_leaf_nodes = Eval(max_leaf_nodes);
   if (!pomax_leaf_nodes) {
      Log() << kFATAL << Form(" MaxLeafNodes = %s... that does not work !! ", max_leaf_nodes.Data())
            << " The options are None or integer."
            << Endl;
   }
   Py_DECREF(pomax_leaf_nodes);

//    bootstrap(kTRUE),
//    oob_score(kFALSE),
//    n_jobs(1),

   PyObject *porandom_state = Eval(random_state);
   if (!porandom_state) {
      Log() << kFATAL << Form(" RandomState = %s... that does not work !! ", random_state.Data())
            << "If int, random_state is the seed used by the random number generator;"
            << "If RandomState instance, random_state is the random number generator;"
            << "If None, the random number generator is the RandomState instance used by `np.random`."
            << Endl;
   }
   Py_DECREF(porandom_state);

//    verbose(0),
//    warm_start(kFALSE),
//    class_weight("None"),
   PyObject *poclass_weight = Eval(class_weight);
   if (!poclass_weight) {
      Log() << kFATAL << Form(" ClassWeight = %s... that does not work !! ", class_weight.Data())
            << "dict, list of dicts, 'auto', 'subsample' or None, optional"
            << Endl;
   }
   Py_DECREF(poclass_weight);

//    n_jobsOVR(1)       // For OneVsRest Classifier
}

//_______________________________________________________________________
void  MethodPyRFOneVsRest::Init()
{
  ProcessOptions();
  _import_array();//require to use numpy arrays

  //Import sklearn
  // Convert the file name to a Python string.
  PyObject *pRFName = PyUnicode_FromString("sklearn.ensemble");
  PyObject *pOVRName = PyUnicode_FromString("sklearn.multiclass");
  
  // USING SINGLE fModule VARIABLE --------------------------
  fModule = PyImport_Import(pOVRName);
  if (!fModule) {
    Log() << kFATAL << "Can't import sklearn.multiclass" << Endl;
    Log() << Endl;
  }
  Py_DECREF(pOVRName);
  
  fModule = PyImport_Import(pRFName);
  if (!fModule) {
    Log() << kFATAL << "Can't import sklearn.ensemble" << Endl;
    Log() << Endl;
  }
  Py_DECREF(pRFName);

  
  // --------------------------------------------------------

  // USING SEPERATE VARIABLES FOR RF AND OVR --------------------------
  // Import the file as a Python module.
  // fRFModule = PyImport_Import(pRFName);
  // fOVRModule = PyImport_Import(pOVRName);
  // Py_DECREF(pRFName);
  // Py_DECREF(pOVRName);

  // if (!fRFModule) {
  //   Log() << kFATAL << "Can't import sklearn.ensemble" << Endl;
  //   Log() << Endl;
  // }
  
  // if (!fOVRModule) {
  //   Log() << kFATAL << "Can't import sklearn.multiclass" << Endl;
  //   Log() << Endl;
  // }
  // ------------------------------------------------------------------


  //Training data
  UInt_t fNvars = Data()->GetNVariables();
  int fNrowsTraining = Data()->GetNTrainingEvents(); //every row is an event, a class type and a weight
  int dims[2];
  dims[0] = fNrowsTraining;
  dims[1] = fNvars;
  fTrainData = (PyArrayObject *)PyArray_FromDims(2, dims, NPY_FLOAT);
  float *TrainData = (float *)(PyArray_DATA(fTrainData));

  fTrainDataClasses = (PyArrayObject *)PyArray_FromDims(1, &fNrowsTraining, NPY_FLOAT);
  float *TrainDataClasses = (float *)(PyArray_DATA(fTrainDataClasses));

  fTrainDataWeights = (PyArrayObject *)PyArray_FromDims(1, &fNrowsTraining, NPY_FLOAT);
  float *TrainDataWeights = (float *)(PyArray_DATA(fTrainDataWeights));

  for (int i = 0; i < fNrowsTraining; i++) {
    const TMVA::Event *e = Data()->GetTrainingEvent(i);
    for (UInt_t j = 0; j < fNvars; j++) {
      TrainData[j + i * fNvars] = e->GetValue(j);
    }
    if (e->GetClass() == TMVA::Types::kSignal) TrainDataClasses[i] = TMVA::Types::kSignal;
    else TrainDataClasses[i] = TMVA::Types::kBackground;

    TrainDataWeights[i] = e->GetWeight();
  }

}

//_______________________________________________________________________
void MethodPyRFOneVsRest::Train()
{

  //NOTE: max_features must have 3 defferents variables int, float and string
  if (max_features == "auto" || max_features == "sqrt" || max_features == "log2"){
    max_features = Form("'%s'", max_features.Data());
  }
  PyObject *pomax_features = Eval(max_features);
  PyObject *pomax_depth = Eval(max_depth);
  PyObject *pomax_leaf_nodes = Eval(max_leaf_nodes);
  PyObject *porandom_state = Eval(random_state);
  PyObject *poclass_weight = Eval(class_weight);

  PyObject *argsRF = Py_BuildValue("(isOiifOOiiiOiiO)", n_estimators, criterion.Data(), pomax_depth, min_samples_split, \
                                  min_samples_leaf, min_weight_fraction_leaf, pomax_features, pomax_leaf_nodes, \
                                  bootstrap, oob_score, n_jobs, porandom_state, verbose, warm_start, poclass_weight);
  Py_DECREF(pomax_depth);
  PyObject_Print(argsRF, stdout, 0);
  std::cout << std::endl;

  // USING SINGLE fModule VARIABLES FOR RF AND OVR ------------------------------------------------------------------------------------------
  PyObject *pRFDict = PyModule_GetDict(fModule);
  PyObject *fRFClassifierClass = PyDict_GetItemString(pRFDict, "RandomForestClassifier");
  //    Log() << kFATAL <<"Train =" <<n_jobs<<Endl;

  // Create an instance of the class
  if (PyCallable_Check(fRFClassifierClass)) {
    //instance
    fClassifier = PyObject_CallObject(fRFClassifierClass, argsRF);
    PyObject_Print(fClassifier, stdout, 0);
    Py_DECREF(argsRF);

    // Feeding RandomForest Classifier to OneVsRest Classifier for multiclass classification.
    PyObject *argsOVR = Py_BuildValue("(Oi)", fClassifier, n_jobsOVR);
    PyObject_Print(argsOVR, stdout, 0);
    std::cout << std::endl;


    PyObject *pOVRName = PyUnicode_FromString("sklearn.multiclass");
  	fModule = PyImport_Import(pOVRName);
  	if(!fModule) {
      Log() << kFATAL << "Can't import sklearn.multiclass" << Endl;
      Log() << Endl;
 	}
  	Py_DECREF(pOVRName);
    PyObject *pOVRDict = PyModule_GetDict(fModule);
    PyObject *fOVRClassifierClass = PyDict_GetItemString(pOVRDict, "OneVsRestClassifier");

    if(PyCallable_Check(fOVRClassifierClass)){
      fClassifier = PyObject_CallObject(fOVRClassifierClass, argsOVR);

      PyObject_Print(fClassifier, stdout, 0);
      Py_DECREF(argsOVR);
    } else {
      PyErr_Print();
      Py_DECREF(pOVRDict);
      Py_DECREF(fOVRClassifierClass);
      Log() << kFATAL << "Can't call function OneVsRestClassifier" << Endl;
      Log() << Endl;
    }

    // OneVsRestClassifier using training dataset features and labels to train.
    fClassifier = PyObject_CallMethod(fClassifier, const_cast<char *>("fit"), const_cast<char *>("(OO)"), fTrainData, fTrainDataClasses);

    if(!fClassifier)
    {
      Log() << kFATAL << "Can't create classifier object from OneVsRestClassifier" << Endl;
      Log() << Endl;  
    }
    if (IsModelPersistence())
    {
      TString path = GetWeightFileDir() + "/PyOVRModel.PyData";
      Log() << Endl;
      Log() << gTools().Color("bold") << "--- Saving State File In:" << gTools().Color("reset") << path << Endl;
      Log() << Endl;
      Serialize(path,fClassifier);
    }

  } else {
    PyErr_Print();
    Py_DECREF(pRFDict);
    Py_DECREF(fRFClassifierClass);
    Log() << kFATAL << "Can't call function RandomForestClassifier" << Endl;
    Log() << Endl;
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------

  // USING SEPERATE fModule  VARIABLES FOR RF AND OVR ------------------------------------------------------------------------------------------
  // PyObject *pRFDict = PyModule_GetDict(fRFModule);
  // PyObject *fRFClassifierClass = PyDict_GetItemString(pRFDict, "RandomForestClassifier");
  // //    Log() << kFATAL <<"Train =" <<n_jobs<<Endl;

  // // Create an instance of the class
  // if (PyCallable_Check(fRFClassifierClass)) {
  //   //instance
  //   fRFClassifier = PyObject_CallObject(fRFClassifierClass, argsRF);
  //   PyObject_Print(fRFClassifier, stdout, 0);
  //   Py_DECREF(argsRF);

  //   // Feeding RandomForest Classifier to OneVsRest Classifier for multiclass classification.
  //   PyObject *argsOVR = Py_BuildValue("(Oi)", fRFClassifier, n_jobsOVR);
  //   PyObject_Print(argsOVR, stdout, 0);
  //   std::cout << std::endl;

  //   PyObject *pOVRDict = PyModule_GetDict(fOVRModule);
  //   PyObject *fOVRClassifierClass = PyDict_GetItemString(pOVRDict, "OneVsRestClassifier");

  //   if(PyCallable_Check(fOVRClassifierClass)){
  //     fOVRClassifier = PyObject_CallObject(fOVRClassifierClass, argsOVR);

  //     PyObject_Print(fOVRClassifier, stdout, 0);
  //     Py_DECREF(argsOVR);
  //   } else {
  //     PyErr_Print();
  //     Py_DECREF(pOVRDict);
  //     Py_DECREF(fOVRClassifierClass);
  //     Log() << kFATAL << "Can't call function OneVsRestClassifier" << Endl;
  //     Log() << Endl;
  //   }

  //   // OneVsRestClassifier using training dataset features and labels to train.
  //   fOVRClassifier = PyObject_CallMethod(fOVRClassifier, const_cast<char *>("fit"), const_cast<char *>("(OO)"), fTrainData, fTrainDataClasses);

  //   if(!fOVRClassifier)
  //   {
  //     Log() << kFATAL << "Can't create classifier object from OneVsRestClassifier" << Endl;
  //     Log() << Endl;  
  //   }
  //   if (IsModelPersistence())
  //   {
  //     TString path = GetWeightFileDir() + "/PyOVRModel.PyData";
  //     Log() << Endl;
  //     Log() << gTools().Color("bold") << "--- Saving State File In:" << gTools().Color("reset") << path << Endl;
  //     Log() << Endl;
  //     Serialize(path,fOVRClassifier);
  //   }

  // } else {
  //   PyErr_Print();
  //   Py_DECREF(pRFDict);
  //   Py_DECREF(fRFClassifierClass);
  //   Log() << kFATAL << "Can't call function RandomForestClassifier" << Endl;
  //   Log() << Endl;
  // }
  // -------------------------------------------------------------------------------------------------------------------------------------------

  // fClassifier = PyObject_CallMethod(fClassifier, const_cast<char *>("fit"), const_cast<char *>("(OOO)"), fTrainData, fTrainDataClasses, fTrainDataWeights);

  // if(!fClassifier)
  // {
  //   Log() << kFATAL << "Can't create classifier object from RandomForestClassifier" << Endl;
  //   Log() << Endl;  
  // }
  // if (IsModelPersistence())
  // {
  //   TString path = GetWeightFileDir() + "/PyRFModel.PyData";
  //   Log() << Endl;
  //   Log() << gTools().Color("bold") << "--- Saving State File In:" << gTools().Color("reset") << path << Endl;
  //   Log() << Endl;
  //   Serialize(path,fClassifier);
  // }
}

//_______________________________________________________________________
void MethodPyRFOneVsRest::TestClassification()
{
  MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodPyRFOneVsRest::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
  // cannot determine error
  NoErrorCalc(errLower, errUpper);

  if (IsModelPersistence()) ReadModelFromFile();

  Double_t mvaValue;
  const TMVA::Event *e = Data()->GetEvent();
  UInt_t nvars = e->GetNVariables();
  int dims[2];
  dims[0] = 1;
  dims[1] = nvars;
  PyArrayObject *pEvent= (PyArrayObject *)PyArray_FromDims(2, dims, NPY_FLOAT);
  float *pValue = (float *)(PyArray_DATA(pEvent));

  for (UInt_t i = 0; i < nvars; i++) pValue[i] = e->GetValue(i);
   
  PyArrayObject *result = (PyArrayObject *)PyObject_CallMethod(fClassifier, const_cast<char *>("predict_proba"), const_cast<char *>("(O)"), pEvent);
  double *proba = (double *)(PyArray_DATA(result));
  mvaValue = proba[0]; //getting signal prob
  Py_DECREF(result);
  Py_DECREF(pEvent);
  return mvaValue;
}

//_______________________________________________________________________
void MethodPyRFOneVsRest::ReadModelFromFile()
{
  if (!PyIsInitialized()) {
    PyInitialize();
  }

  TString path = GetWeightFileDir() + "/PyOVRModel.PyData";
  Log() << Endl;
  Log() << gTools().Color("bold") << "--- Loading State File From:" << gTools().Color("reset") << path << Endl;
  Log() << Endl;
  UnSerialize(path,&fClassifier);
  if(!fClassifier)
  {
    Log() << kFATAL << "Can't load OneVsRestRandomForestClassifier from Serialized data." << Endl;
    Log() << Endl;
  }
}

//_______________________________________________________________________
void MethodPyRFOneVsRest::GetHelpMessage() const
{
  // get help message text
  //
  // typical length of text line:
  //         "|--------------------------------------------------------------|"
  Log() << Endl;
  Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
  Log() << Endl;
  Log() << "Decision Trees and Rule-Based Models " << Endl;
  Log() << Endl;
  Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
  Log() << Endl;
  Log() << Endl;
  Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
  Log() << Endl;
  Log() << "<None>" << Endl;
}

