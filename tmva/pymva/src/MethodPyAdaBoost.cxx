// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyAdaBoost                                                      *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      AdaBoost      Classifier from Scikit learn                               *
 *                                                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/

#include <Python.h> // Needs to be included first to avoid redefinition of _POSIX_C_SOURCE
#include "TMVA/MethodPyAdaBoost.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include "TMVA/IMethod.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDF.h"
#include "TMVA/Ranking.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Timer.h"
#include "TMVA/VariableTransformBase.h"
#include "TMVA/Results.h"

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

REGISTER_METHOD(PyAdaBoost)

ClassImp(MethodPyAdaBoost);

//_______________________________________________________________________
MethodPyAdaBoost::MethodPyAdaBoost(const TString &jobName,
                                   const TString &methodTitle,
                                   DataSetInfo &dsi,
                                   const TString &theOption) :
   PyMethodBase(jobName, Types::kPyAdaBoost, methodTitle, dsi, theOption),
   fBaseEstimator("None"),
   fNestimators(50),
   fLearningRate(1.0),
   fAlgorithm("SAMME.R"),
   fRandomState("None")
{
}

//_______________________________________________________________________
MethodPyAdaBoost::MethodPyAdaBoost(DataSetInfo &theData,
                                   const TString &theWeightFile) :
   PyMethodBase(Types::kPyAdaBoost, theData, theWeightFile),
   fBaseEstimator("None"),
   fNestimators(50),
   fLearningRate(1.0),
   fAlgorithm("SAMME.R"),
   fRandomState("None")
{
}

//_______________________________________________________________________
MethodPyAdaBoost::~MethodPyAdaBoost(void)
{
}

//_______________________________________________________________________
Bool_t MethodPyAdaBoost::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass && numberClasses >= 2) return kTRUE;
   return kFALSE;
}

//_______________________________________________________________________
void MethodPyAdaBoost::DeclareOptions()
{
   MethodBase::DeclareCompatibilityOptions();

   DeclareOptionRef(fBaseEstimator, "BaseEstimator", "object, optional (default=DecisionTreeClassifier)\
      The base estimator from which the boosted ensemble is built.\
      Support for sample weighting is required, as well as proper `classes_`\
      and `n_classes_` attributes.");

   DeclareOptionRef(fNestimators, "NEstimators", "integer, optional (default=50)\
      The maximum number of estimators at which boosting is terminated.\
      In case of perfect fit, the learning procedure is stopped early.");

   DeclareOptionRef(fLearningRate, "LearningRate", "float, optional (default=1.)\
      Learning rate shrinks the contribution of each classifier by\
      ``learning_rate``. There is a trade-off between ``learning_rate`` and\
      ``n_estimators``.");

   DeclareOptionRef(fAlgorithm, "Algorithm", "{'SAMME', 'SAMME.R'}, optional (default='SAMME.R')\
      If 'SAMME.R' then use the SAMME.R real boosting algorithm.\
      ``base_estimator`` must support calculation of class probabilities.\
      If 'SAMME' then use the SAMME discrete boosting algorithm.\
      The SAMME.R algorithm typically converges faster than SAMME,\
      achieving a lower test error with fewer boosting iterations.");

   DeclareOptionRef(fRandomState, "RandomState", "int, RandomState instance or None, optional (default=None)\
      If int, random_state is the seed used by the random number generator;\
      If RandomState instance, random_state is the random number generator;\
      If None, the random number generator is the RandomState instance used\
      by `np.random`.");

   DeclareOptionRef(fFilenameClassifier, "FilenameClassifier",
      "Store trained classifier in this file");
}

//_______________________________________________________________________
// Check options and load them to local python namespace
void MethodPyAdaBoost::ProcessOptions()
{
   pBaseEstimator = Eval(fBaseEstimator);
   if (!pBaseEstimator) {
      Log() << kFATAL << Form("BaseEstimator = %s ... that does not work!", fBaseEstimator.Data())
            << " The options are Object or None." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "baseEstimator", pBaseEstimator);

   if (fNestimators <= 0) {
      Log() << kFATAL << "NEstimators <=0 ... that does not work!" << Endl;
   }
   pNestimators = Eval(Form("%i", fNestimators));
   PyDict_SetItemString(fLocalNS, "nEstimators", pNestimators);

   if (fLearningRate <= 0) {
      Log() << kFATAL << "LearningRate <=0 ... that does not work!" << Endl;
   }
   pLearningRate = Eval(Form("%f", fLearningRate));
   PyDict_SetItemString(fLocalNS, "learningRate", pLearningRate);

   if (fAlgorithm != "SAMME" && fAlgorithm != "SAMME.R") {
      Log() << kFATAL << Form("Algorithm = %s ... that does not work!", fAlgorithm.Data())
            << " The options are SAMME of SAMME.R." << Endl;
   }
   pAlgorithm = Eval(Form("'%s'", fAlgorithm.Data()));
   PyDict_SetItemString(fLocalNS, "algorithm", pAlgorithm);

   pRandomState = Eval(fRandomState);
   if (!pRandomState) {
      Log() << kFATAL << Form(" RandomState = %s... that does not work !! ", fRandomState.Data())
            << "If int, random_state is the seed used by the random number generator;"
            << "If RandomState instance, random_state is the random number generator;"
            << "If None, the random number generator is the RandomState instance used by `np.random`." << Endl;
   }
   PyDict_SetItemString(fLocalNS, "randomState", pRandomState);

   // If no filename is given, set default
   if(fFilenameClassifier.IsNull()) {
      fFilenameClassifier = GetWeightFileDir() + "/PyAdaBoostModel_" + GetName() + ".PyData";
   }
}

//_______________________________________________________________________
void MethodPyAdaBoost::Init()
{
   TMVA::Internal::PyGILRAII raii;
   _import_array(); //require to use numpy arrays

   // Check options and load them to local python namespace
   ProcessOptions();

   // Import module for ada boost classifier
   PyRunString("import sklearn.ensemble");

   // Get data properties
   fNvars = GetNVariables();
   fNoutputs = DataInfo().GetNClasses();
}

//_______________________________________________________________________
void MethodPyAdaBoost::Train()
{
   // Load training data (data, classes, weights) to python arrays
   int fNrowsTraining = Data()->GetNTrainingEvents(); //every row is an event, a class type and a weight
   npy_intp dimsData[2];
   dimsData[0] = fNrowsTraining;
   dimsData[1] = fNvars;
   PyArrayObject *  fTrainData = (PyArrayObject *)PyArray_SimpleNew(2, dimsData, NPY_FLOAT);
   PyDict_SetItemString(fLocalNS, "trainData", (PyObject*)fTrainData);
   float *TrainData = (float *)(PyArray_DATA(fTrainData));

   npy_intp dimsClasses = (npy_intp) fNrowsTraining;
   PyArrayObject *  fTrainDataClasses = (PyArrayObject *)PyArray_SimpleNew(1, &dimsClasses, NPY_FLOAT);
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
   PyRunString("classifier = sklearn.ensemble.AdaBoostClassifier(base_estimator=baseEstimator, n_estimators=nEstimators, learning_rate=learningRate, algorithm=algorithm, random_state=randomState)",
      "Failed to setup classifier");

   // Fit classifier
   // NOTE: We dump the output to a variable so that the call does not pollute stdout
   PyRunString("dump = classifier.fit(trainData, trainDataClasses, trainDataWeights)", "Failed to train classifier");

   // Store classifier
   fClassifier = PyDict_GetItemString(fLocalNS, "classifier");
   if(fClassifier == 0) {
      Log() << kFATAL << "Can't create classifier object from AdaBoostClassifier" << Endl;
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
void MethodPyAdaBoost::TestClassification()
{
   MethodBase::TestClassification();
}

//_______________________________________________________________________
std::vector<Double_t> MethodPyAdaBoost::GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress)
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
Double_t MethodPyAdaBoost::GetMvaValue(Double_t *errLower, Double_t *errUpper)
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
std::vector<Float_t>& MethodPyAdaBoost::GetMulticlassValues()
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

   return classValues;
}

//_______________________________________________________________________
void MethodPyAdaBoost::ReadModelFromFile()
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
const Ranking* MethodPyAdaBoost::CreateRanking()
{
   // Get feature importance from classifier as an array with length equal
   // number of variables, higher value signals a higher importance
   PyArrayObject* pRanking = (PyArrayObject*) PyObject_GetAttrString(fClassifier, "feature_importances_");
   // The python object is null if the base estimator does not support
   // variable ranking. Then, return NULL, which disables ranking.
   if(pRanking == 0) return NULL;

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
void MethodPyAdaBoost::GetHelpMessage() const
{
   // typical length of text line:
   //         "|--------------------------------------------------------------|"
   Log() << "An AdaBoost classifier is a meta-estimator that begins by fitting" << Endl;
   Log() << "a classifier on the original dataset and then fits additional copies" << Endl;
   Log() << "of the classifier on the same dataset but where the weights of incorrectly" << Endl;
   Log() << "classified instances are adjusted such that subsequent classifiers focus" << Endl;
   Log() << "more on difficult cases." << Endl;
   Log() << Endl;
   Log() << "Check out the scikit-learn documentation for more information." << Endl;
}
