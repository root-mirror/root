// @(#)root/tmva/pymva $Id$
// Author: Stefan Wunsch, 2016

#include <Python.h>
#include "TMVA/MethodPyKeras.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "TMVA/Types.h"
#include "TMVA/Config.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/Results.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/VariableTransformBase.h"
#include "TMVA/Tools.h"
#include "TMVA/Timer.h"
#include "TSystem.h"

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

REGISTER_METHOD(PyKeras)

ClassImp(MethodPyKeras);

MethodPyKeras::MethodPyKeras(const TString &jobName, const TString &methodTitle, DataSetInfo &dsi, const TString &theOption)
   : PyMethodBase(jobName, Types::kPyKeras, methodTitle, dsi, theOption) {
   fNumEpochs = 10;
   fBatchSize = 100;
   fVerbose = 1;
   fContinueTraining = false;
   fSaveBestOnly = true;
   fTriesEarlyStopping = -1;
   fLearningRateSchedule = ""; // empty string deactivates learning rate scheduler
   fFilenameTrainedModel = ""; // empty string sets output model filename to default (in weights/)
   fTensorBoard = "";          // empty string deactivates TensorBoard callback
}

MethodPyKeras::MethodPyKeras(DataSetInfo &theData, const TString &theWeightFile)
    : PyMethodBase(Types::kPyKeras, theData, theWeightFile) {
   fNumEpochs = 10;
   fNumThreads = 0;
   fBatchSize = 100;
   fVerbose = 1;
   fContinueTraining = false;
   fSaveBestOnly = true;
   fTriesEarlyStopping = -1;
   fLearningRateSchedule = ""; // empty string deactivates learning rate scheduler
   fFilenameTrainedModel = ""; // empty string sets output model filename to default (in weights/)
   fTensorBoard = "";          // empty string deactivates TensorBoard callback
}

MethodPyKeras::~MethodPyKeras() {
}

Bool_t MethodPyKeras::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t) {
   if (type == Types::kRegression) return kTRUE;
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass && numberClasses >= 2) return kTRUE;
   return kFALSE;
}

///////////////////////////////////////////////////////////////////////////////

void MethodPyKeras::DeclareOptions() {
   DeclareOptionRef(fFilenameModel, "FilenameModel", "Filename of the initial Keras model");
   DeclareOptionRef(fFilenameTrainedModel, "FilenameTrainedModel", "Filename of the trained output Keras model");
   DeclareOptionRef(fBatchSize, "BatchSize", "Training batch size");
   DeclareOptionRef(fNumEpochs, "NumEpochs", "Number of training epochs");
   DeclareOptionRef(fNumThreads, "NumThreads", "Number of CPU threads (only for Tensorflow backend)");
   DeclareOptionRef(fGpuOptions, "GpuOptions", "GPU options for tensorflow, such as allow_growth");
   DeclareOptionRef(fUseTFKeras, "tf.keras", "Use tensorflow from Keras");
   DeclareOptionRef(fVerbose, "Verbose", "Keras verbosity during training");
   DeclareOptionRef(fContinueTraining, "ContinueTraining", "Load weights from previous training");
   DeclareOptionRef(fSaveBestOnly, "SaveBestOnly", "Store only weights with smallest validation loss");
   DeclareOptionRef(fTriesEarlyStopping, "TriesEarlyStopping", "Number of epochs with no improvement in validation loss after which training will be stopped. The default or a negative number deactivates this option.");
   DeclareOptionRef(fLearningRateSchedule, "LearningRateSchedule", "Set new learning rate during training at specific epochs, e.g., \"50,0.01;70,0.005\"");
   DeclareOptionRef(fTensorBoard, "TensorBoard",
                    "Write a log during training to visualize and monitor the training performance with TensorBoard");

   DeclareOptionRef(fNumValidationString = "20%", "ValidationSize", "Part of the training data to use for validation. "
                    "Specify as 0.2 or 20% to use a fifth of the data set as validation set. "
                    "Specify as 100 to use exactly 100 events. (Default: 20%)");
   DeclareOptionRef(fUserCodeName = "", "UserCode",
                    "Optional python code provided by the user to be executed before loading the Keras model");
}

////////////////////////////////////////////////////////////////////////////////
/// Validation of the ValidationSize option. Allowed formats are 20%, 0.2 and
/// 100 etc.
///    - 20% and 0.2 selects 20% of the training set as validation data.
///    - 100 selects 100 events as the validation data.
///
/// @return number of samples in validation set
///
UInt_t TMVA::MethodPyKeras::GetNumValidationSamples()
{
   Int_t nValidationSamples = 0;
   UInt_t trainingSetSize = GetEventCollection(Types::kTraining).size();

   // Parsing + Validation
   // --------------------
   if (fNumValidationString.EndsWith("%")) {
      // Relative spec. format 20%
      TString intValStr = TString(fNumValidationString.Strip(TString::kTrailing, '%'));

      if (intValStr.IsFloat()) {
         Double_t valSizeAsDouble = fNumValidationString.Atof() / 100.0;
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         Log() << kFATAL << "Cannot parse number \"" << fNumValidationString
               << "\". Expected string like \"20%\" or \"20.0%\"." << Endl;
      }
   } else if (fNumValidationString.IsFloat()) {
      Double_t valSizeAsDouble = fNumValidationString.Atof();

      if (valSizeAsDouble < 1.0) {
         // Relative spec. format 0.2
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         // Absolute spec format 100 or 100.0
         nValidationSamples = valSizeAsDouble;
      }
   } else {
      Log() << kFATAL << "Cannot parse number \"" << fNumValidationString << "\". Expected string like \"0.2\" or \"100\"."
            << Endl;
   }

   // Value validation
   // ----------------
   if (nValidationSamples < 0) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString << "\" is negative." << Endl;
   }

   if (nValidationSamples == 0) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString << "\" is zero." << Endl;
   }

   if (nValidationSamples >= (Int_t)trainingSetSize) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString
            << "\" is larger than or equal in size to training set (size=\"" << trainingSetSize << "\")." << Endl;
   }

   return nValidationSamples;
}

/// Function processing the options
/// This is called only when creating the method before training not when
/// readinf from XML file. Called from MethodBase::ProcessSetup
/// that is called from Factory::BookMethod
void MethodPyKeras::ProcessOptions() {

   // Set default filename for trained model if option is not used
   if (fFilenameTrainedModel.IsNull()) {
      fFilenameTrainedModel = GetWeightFileDir() + "/TrainedModel_" + GetName() + ".h5";
   }

   // Setup model, either the initial model from `fFilenameModel` or
   // the trained model from `fFilenameTrainedModel`
   if (fContinueTraining) Log() << kINFO << "Continue training with trained model" << Endl;
   SetupKerasModel(fContinueTraining);
}

void MethodPyKeras::SetupKerasModel(bool loadTrainedModel) {

   // initialize first Keras. This is done only here when class has
   // all state variable set from options or read from XML file
   // Import Keras

   if (fUseTFKeras)
      Log() << kINFO << "Setting up tf.keras" << Endl;
   else
      Log() << kINFO << "Setting up keras with " << gSystem->Getenv("KERAS_BACKEND") << " backend" << Endl;

   bool useTFBackend = kFALSE;
   bool kerasIsCompatible = kTRUE;
   bool kerasIsPresent = kFALSE;

   if (!fUseTFKeras) {
      auto ret  = PyRun_String("import keras", Py_single_input, fGlobalNS, fLocalNS);
      // need importing also in global namespace
      if (ret != nullptr) ret = PyRun_String("import keras", Py_single_input, fGlobalNS, fGlobalNS);
      if (ret != nullptr)
         kerasIsPresent = kTRUE;
      if (kerasIsPresent) {
         // check compatibility with tensorflow
         if (GetKerasBackend() == kTensorFlow ) {
            useTFBackend = kTRUE;

            PyRunString("keras_major_version = int(keras.__version__.split('.')[0])");
            PyRunString("keras_minor_version = int(keras.__version__.split('.')[1])");
            PyObject *pyKerasMajorVersion = PyDict_GetItemString(fLocalNS, "keras_major_version");
            PyObject *pyKerasMinorVersion = PyDict_GetItemString(fLocalNS, "keras_minor_version");
            int kerasMajorVersion = PyLong_AsLong(pyKerasMajorVersion);
            int kerasMinorVersion = PyLong_AsLong(pyKerasMinorVersion);
            Log() << kINFO << "Using Keras version " << kerasMajorVersion << "." << kerasMinorVersion << Endl;
            // only version 2.3 is latest multi-backend version.
            // version 2.4 is just tf.keras and should not be used in standalone and will not work in this workflow
            // see https://github.com/keras-team/keras/releases/tag/2.4.0
            // for example variable  keras.backend.tensorflow_backend will not exist anymore in keras 2.4
            kerasIsCompatible = (kerasMajorVersion >= 2 && kerasMinorVersion == 3);

         }
      } else {
         // Keras is not found. try tyo use tf.keras
         Log() << kINFO << "Keras is not found. Trying using tf.keras" << Endl;
         fUseTFKeras = 1;
      }
   }

   // import Tensoprflow (if requested or because is keras backend)
   if (fUseTFKeras || useTFBackend) {
      auto ret = PyRun_String("import tensorflow as tf", Py_single_input, fGlobalNS, fLocalNS);
      if (ret != nullptr) ret = PyRun_String("import tensorflow as tf", Py_single_input, fGlobalNS, fGlobalNS);
      if (ret == nullptr) {
         Log() << kFATAL << "Importing TensorFlow failed" << Endl;
      }
      // check tensorflow version
      PyRunString("tf_major_version = int(tf.__version__.split('.')[0])");
      PyObject *pyTfVersion = PyDict_GetItemString(fLocalNS, "tf_major_version");
      int tfVersion = PyLong_AsLong(pyTfVersion);
      Log() << kINFO << "Using TensorFlow version " << tfVersion << Endl;

      if (tfVersion < 2) {
         if (fUseTFKeras == 1) {
            Log() << kWARNING << "Using TensorFlow version 1.x which does not contain tf.keras - use then TensorFlow as Keras backend" << Endl;
            fUseTFKeras = kFALSE;
            // case when Keras was not found
            if (!kerasIsPresent) {
               Log() << kFATAL << "Keras is not present and not a suitable TensorFlow version is found " << Endl;
               return;
            }
         }
      }
      else {
         // using version larger than 2.0 - can use tf.keras
         if (!kerasIsCompatible) {
            Log() << kWARNING << "The Keras version is not compatible with TensorFlow 2. Use instead tf.keras" << Endl;
            fUseTFKeras = 1;
         }
      }

      // if keras 2.3 and tensorflow 2 are found. Use tf.keras or keras ?
      // at the moment default is tf.keras=false to keep compatibility
      // but this might change in future releases
      if (fUseTFKeras) {
         Log() << kINFO << "Use Keras version from TensorFlow : tf.keras" << Endl;
         fKerasString = "tf.keras";
         PyRunString("K = tf.keras.backend");
         PyRun_String("K = tf.keras.backend", Py_single_input, fGlobalNS, fGlobalNS);
      }
      else {
         Log() << kINFO << "Use TensorFlow as Keras backend" << Endl;
         fKerasString = "keras";
         PyRunString("from keras.backend import tensorflow_backend as K");
         PyRun_String("from keras.backend import tensorflow_backend as K", Py_single_input, fGlobalNS, fGlobalNS);
      }

      // extra options for tensorflow
      // use different naming in tf2 for ConfigProto and Session
      TString configProto = (tfVersion >= 2) ? "tf.compat.v1.ConfigProto" : "tf.ConfigProto";
      TString session = (tfVersion >= 2) ? "tf.compat.v1.Session" : "tf.Session";

      // in case specify number of threads
      int num_threads = fNumThreads;
      if (num_threads > 0) {
         Log() << kINFO << "Setting the CPU number of threads =  " << num_threads << Endl;

         PyRunString(
            TString::Format("session_conf = %s(intra_op_parallelism_threads=%d,inter_op_parallelism_threads=%d)",
                            configProto.Data(), num_threads, num_threads));
      } else
         PyRunString(TString::Format("session_conf = %s()", configProto.Data()));

      // applying GPU options such as allow_growth=True to avoid allocating all memory on GPU
      // that prevents running later TMVA-GPU
      // Also new Nvidia RTX cards (e.g. RTX 2070)  require this option
      if (!fGpuOptions.IsNull()) {
         TObjArray *optlist = fGpuOptions.Tokenize(",");
         for (int item = 0; item < optlist->GetEntries(); ++item) {
            Log() << kINFO << "Applying GPU option:  gpu_options." << optlist->At(item)->GetName() << Endl;
            PyRunString(TString::Format("session_conf.gpu_options.%s", optlist->At(item)->GetName()));
         }
      }
      PyRunString(TString::Format("sess = %s(config=session_conf)", session.Data()));

      if (tfVersion < 2) {
         PyRunString("K.set_session(sess)");
      } else {
         PyRunString("tf.compat.v1.keras.backend.set_session(sess)");
      }
   }
   // case not using a Tensorflow backend
   else {
      fKerasString = "keras";
      if (fNumThreads > 0)
         Log() << kWARNING << "Cannot set the given " << fNumThreads << " threads when not using tensorflow as  backend"
               << Endl;
      if (!fGpuOptions.IsNull()) {
         Log() << kWARNING << "Cannot set the given GPU option " << fGpuOptions
               << " when not using tensorflow as  backend" << Endl;
      }
   }

   /*
    * Load Keras model from file
    */

   Log() << kINFO << " Loading Keras Model " << Endl;

   PyRunString("load_model_custom_objects=None");



   if (!fUserCodeName.IsNull()) {
      Log() << kINFO << " Executing user initialization code from  " << fUserCodeName << Endl;


      // run some python code provided by user for model initialization if needed
      TString cmd = "exec(open('" + fUserCodeName + "').read())";
      TString errmsg = "Error executing the provided user code";
      PyRunString(cmd, errmsg);

      PyRunString("print('custom objects for loading model : ',load_model_custom_objects)");
   }

   // Load initial model or already trained model
   TString filenameLoadModel;
   if (loadTrainedModel) {
      filenameLoadModel = fFilenameTrainedModel;
   }
   else {
      filenameLoadModel = fFilenameModel;
   }

   PyRunString("model = " + fKerasString + ".models.load_model('" + filenameLoadModel +
                     "', custom_objects=load_model_custom_objects)", "Failed to load Keras model from file: " + filenameLoadModel);

   Log() << kINFO << "Loaded model from file: " << filenameLoadModel << Endl;


   /*
    * Init variables and weights
    */

   // Get variables, classes and target numbers
   fNVars = GetNVariables();
   if (GetAnalysisType() == Types::kClassification || GetAnalysisType() == Types::kMulticlass) fNOutputs = DataInfo().GetNClasses();
   else if (GetAnalysisType() == Types::kRegression) fNOutputs = DataInfo().GetNTargets();
   else Log() << kFATAL << "Selected analysis type is not implemented" << Endl;

   // Init evaluation (needed for getMvaValue)
   fVals = new float[fNVars]; // holds values used for classification and regression
   npy_intp dimsVals[2] = {(npy_intp)1, (npy_intp)fNVars};
   PyArrayObject* pVals = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsVals, NPY_FLOAT, (void*)fVals);
   PyDict_SetItemString(fLocalNS, "vals", (PyObject*)pVals);

   fOutput.resize(fNOutputs); // holds classification probabilities or regression output
   npy_intp dimsOutput[2] = {(npy_intp)1, (npy_intp)fNOutputs};
   PyArrayObject* pOutput = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsOutput, NPY_FLOAT, (void*)&fOutput[0]);
   PyDict_SetItemString(fLocalNS, "output", (PyObject*)pOutput);

   // Mark the model as setup
   fModelIsSetup = true;
}

/// Initialization function called from MethodBase::SetupMethod()
/// Note that option string are not yet filled with their values.
/// This is done before ProcessOption method or after reading from XML file
void MethodPyKeras::Init() {

   TMVA::Internal::PyGILRAII raii;

   if (!PyIsInitialized()) {
      Log() << kFATAL << "Python is not initialized" << Endl;
   }
   _import_array(); // required to use numpy arrays

   // NOTE: sys.argv has to be cleared because otherwise TensorFlow breaks
   PyRunString("import sys; sys.argv = ['']", "Set sys.argv failed");

   // Set flag that model is not setup
   fModelIsSetup = false;
}

void MethodPyKeras::Train() {

   if(!fModelIsSetup) Log() << kFATAL << "Model is not setup for training" << Endl;

   /*
    * Load training data to numpy array
    */

   UInt_t nAllEvents = Data()->GetNTrainingEvents();
   UInt_t nValEvents = GetNumValidationSamples();
   UInt_t nTrainingEvents = nAllEvents - nValEvents;

   Log() << kINFO << "Split TMVA training data in " << nTrainingEvents << " training events and "
         << nValEvents << " validation events" << Endl;

   float* trainDataX = new float[nTrainingEvents*fNVars];
   float* trainDataY = new float[nTrainingEvents*fNOutputs];
   float* trainDataWeights = new float[nTrainingEvents];
   for (UInt_t i=0; i<nTrainingEvents; i++) {
      const TMVA::Event* e = GetTrainingEvent(i);
      // Fill variables
      for (UInt_t j=0; j<fNVars; j++) {
         trainDataX[j + i*fNVars] = e->GetValue(j);
      }
      // Fill targets
      // NOTE: For classification, convert class number in one-hot vector,
      // e.g., 1 -> [0, 1] or 0 -> [1, 0] for binary classification
      if (GetAnalysisType() == Types::kClassification || GetAnalysisType() == Types::kMulticlass) {
         for (UInt_t j=0; j<fNOutputs; j++) {
            trainDataY[j + i*fNOutputs] = 0;
         }
         trainDataY[e->GetClass() + i*fNOutputs] = 1;
      }
      else if (GetAnalysisType() == Types::kRegression) {
         for (UInt_t j=0; j<fNOutputs; j++) {
            trainDataY[j + i*fNOutputs] = e->GetTarget(j);
         }
      }
      else Log() << kFATAL << "Can not fill target vector because analysis type is not known" << Endl;
      // Fill weights
      // NOTE: If no weight branch is given, this defaults to ones for all events
      trainDataWeights[i] = e->GetWeight();
   }

   npy_intp dimsTrainX[2] = {(npy_intp)nTrainingEvents, (npy_intp)fNVars};
   npy_intp dimsTrainY[2] = {(npy_intp)nTrainingEvents, (npy_intp)fNOutputs};
   npy_intp dimsTrainWeights[1] = {(npy_intp)nTrainingEvents};
   PyArrayObject* pTrainDataX = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsTrainX, NPY_FLOAT, (void*)trainDataX);
   PyArrayObject* pTrainDataY = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsTrainY, NPY_FLOAT, (void*)trainDataY);
   PyArrayObject* pTrainDataWeights = (PyArrayObject*)PyArray_SimpleNewFromData(1, dimsTrainWeights, NPY_FLOAT, (void*)trainDataWeights);
   PyDict_SetItemString(fLocalNS, "trainX", (PyObject*)pTrainDataX);
   PyDict_SetItemString(fLocalNS, "trainY", (PyObject*)pTrainDataY);
   PyDict_SetItemString(fLocalNS, "trainWeights", (PyObject*)pTrainDataWeights);

   /*
    * Load validation data to numpy array
    */

   // NOTE: from TMVA,  we get the validation data as a subset of all the training data
   // we will not use test data for validation. They will be used for the real testing


   float* valDataX = new float[nValEvents*fNVars];
   float* valDataY = new float[nValEvents*fNOutputs];
   float* valDataWeights = new float[nValEvents];
   //validation events follows the trainig one in the TMVA training vector
   for (UInt_t i=0; i< nValEvents ; i++) {
      UInt_t ievt = nTrainingEvents + i; // TMVA event index
      const TMVA::Event* e = GetTrainingEvent(ievt);
      // Fill variables
      for (UInt_t j=0; j<fNVars; j++) {
         valDataX[j + i*fNVars] = e->GetValue(j);
      }
      // Fill targets
      if (GetAnalysisType() == Types::kClassification || GetAnalysisType() == Types::kMulticlass) {
         for (UInt_t j=0; j<fNOutputs; j++) {
            valDataY[j + i*fNOutputs] = 0;
         }
         valDataY[e->GetClass() + i*fNOutputs] = 1;
      }
      else if (GetAnalysisType() == Types::kRegression) {
         for (UInt_t j=0; j<fNOutputs; j++) {
            valDataY[j + i*fNOutputs] = e->GetTarget(j);
         }
      }
      else Log() << kFATAL << "Can not fill target vector because analysis type is not known" << Endl;
      // Fill weights
      valDataWeights[i] = e->GetWeight();
   }

   npy_intp dimsValX[2] = {(npy_intp)nValEvents, (npy_intp)fNVars};
   npy_intp dimsValY[2] = {(npy_intp)nValEvents, (npy_intp)fNOutputs};
   npy_intp dimsValWeights[1] = {(npy_intp)nValEvents};
   PyArrayObject* pValDataX = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsValX, NPY_FLOAT, (void*)valDataX);
   PyArrayObject* pValDataY = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsValY, NPY_FLOAT, (void*)valDataY);
   PyArrayObject* pValDataWeights = (PyArrayObject*)PyArray_SimpleNewFromData(1, dimsValWeights, NPY_FLOAT, (void*)valDataWeights);
   PyDict_SetItemString(fLocalNS, "valX", (PyObject*)pValDataX);
   PyDict_SetItemString(fLocalNS, "valY", (PyObject*)pValDataY);
   PyDict_SetItemString(fLocalNS, "valWeights", (PyObject*)pValDataWeights);

   /*
    * Train Keras model
    */
   Log() << kINFO << "Training Model Summary" << Endl;
   PyRunString("model.summary()");

   // Setup parameters

   PyObject* pBatchSize = PyLong_FromLong(fBatchSize);
   PyObject* pNumEpochs = PyLong_FromLong(fNumEpochs);
   PyObject* pVerbose = PyLong_FromLong(fVerbose);
   PyDict_SetItemString(fLocalNS, "batchSize", pBatchSize);
   PyDict_SetItemString(fLocalNS, "numEpochs", pNumEpochs);
   PyDict_SetItemString(fLocalNS, "verbose", pVerbose);

   // Setup training callbacks
   PyRunString("callbacks = []");

   // Callback: Save only weights with smallest validation loss
   if (fSaveBestOnly) {
      PyRunString("callbacks.append(" + fKerasString +".callbacks.ModelCheckpoint('"+fFilenameTrainedModel+"', monitor='val_loss', verbose=verbose, save_best_only=True, mode='auto'))", "Failed to setup training callback: SaveBestOnly");
      Log() << kINFO << "Option SaveBestOnly: Only model weights with smallest validation loss will be stored" << Endl;
   }

   // Callback: Stop training early if no improvement in validation loss is observed
   if (fTriesEarlyStopping>=0) {
      TString tries;
      tries.Form("%i", fTriesEarlyStopping);
      PyRunString("callbacks.append(" + fKerasString + ".callbacks.EarlyStopping(monitor='val_loss', patience="+tries+", verbose=verbose, mode='auto'))", "Failed to setup training callback: TriesEarlyStopping");
      Log() << kINFO << "Option TriesEarlyStopping: Training will stop after " << tries << " number of epochs with no improvement of validation loss" << Endl;
   }

   // Callback: Learning rate scheduler
   if (fLearningRateSchedule!="") {
      // Setup a python dictionary with the desired learning rate steps
      PyRunString("strScheduleSteps = '"+fLearningRateSchedule+"'\n"
                  "schedulerSteps = {}\n"
                  "for c in strScheduleSteps.split(';'):\n"
                  "    x = c.split(',')\n"
                  "    schedulerSteps[int(x[0])] = float(x[1])\n",
                  "Failed to setup steps for scheduler function from string: "+fLearningRateSchedule,
                   Py_file_input);
      // Set scheduler function as piecewise function with given steps
      PyRunString("def schedule(epoch, model=model, schedulerSteps=schedulerSteps):\n"
                  "    if epoch in schedulerSteps: return float(schedulerSteps[epoch])\n"
                  "    else: return float(model.optimizer.lr.get_value())\n",
                  "Failed to setup scheduler function with string: "+fLearningRateSchedule,
                  Py_file_input);
      // Setup callback
      PyRunString("callbacks.append(" + fKerasString + ".callbacks.LearningRateScheduler(schedule))",
                  "Failed to setup training callback: LearningRateSchedule");
      Log() << kINFO << "Option LearningRateSchedule: Set learning rate during training: " << fLearningRateSchedule << Endl;
   }

   // Callback: TensorBoard
   if (fTensorBoard != "") {
      TString logdir = TString("'") + fTensorBoard + TString("'");
      PyRunString(
         "callbacks.append(" + fKerasString + ".callbacks.TensorBoard(log_dir=" + logdir +
            ", histogram_freq=0, batch_size=batchSize, write_graph=True, write_grads=False, write_images=False))",
         "Failed to setup training callback: TensorBoard");
      Log() << kINFO << "Option TensorBoard: Log files for training monitoring are stored in: " << logdir << Endl;
   }

   // Train model
   PyRunString("history = model.fit(trainX, trainY, sample_weight=trainWeights, batch_size=batchSize, epochs=numEpochs, verbose=verbose, validation_data=(valX, valY, valWeights), callbacks=callbacks)",
               "Failed to train model");


   std::vector<float> fHistory; // Hold training history  (val_acc or loss etc)
   fHistory.resize(fNumEpochs); // holds training loss or accuracy output
   npy_intp dimsHistory[1] = { (npy_intp)fNumEpochs};
   PyArrayObject* pHistory = (PyArrayObject*)PyArray_SimpleNewFromData(1, dimsHistory, NPY_FLOAT, (void*)&fHistory[0]);
   PyDict_SetItemString(fLocalNS, "HistoryOutput", (PyObject*)pHistory);

   // Store training history data
   Int_t iHis=0;
   PyRunString("number_of_keys=len(history.history.keys())");
   PyObject* PyNkeys=PyDict_GetItemString(fLocalNS, "number_of_keys");
   int nkeys=PyLong_AsLong(PyNkeys);
   for (iHis=0; iHis<nkeys; iHis++) {

      PyRunString(TString::Format("copy_string=str(list(history.history.keys())[%d])",iHis));
      PyObject* stra=PyDict_GetItemString(fLocalNS, "copy_string");
      if(!stra) break;
#if PY_MAJOR_VERSION < 3   // for Python2
      const char *stra_name = PyBytes_AsString(stra);
      // need to add string delimiter for Python2
      TString sname = TString::Format("'%s'",stra_name);
      const char * name = sname.Data();
#else   // for Python3
      PyObject* repr = PyObject_Repr(stra);
      PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
      const char *name = PyBytes_AsString(str);
#endif

      Log() << kINFO << "Getting training history for item:" << iHis << " name = " << name << Endl;
      PyRunString(TString::Format("for i,p in enumerate(history.history[%s]):\n   HistoryOutput[i]=p\n",name),
                  TString::Format("Failed to get %s from training history",name));
      for (size_t i=0; i<fHistory.size(); i++)
         fTrainHistory.AddValue(name,i+1,fHistory[i]);

   }
//#endif

   /*
    * Store trained model to file (only if option 'SaveBestOnly' is NOT activated,
    * because we do not want to override the best model checkpoint)
    */

   if (!fSaveBestOnly) {
      PyRunString("model.save('"+fFilenameTrainedModel+"', overwrite=True)",
                  "Failed to save trained model: "+fFilenameTrainedModel);
      Log() << kINFO << "Trained model written to file: " << fFilenameTrainedModel << Endl;
   }

   /*
    * Clean-up
    */

   delete[] trainDataX;
   delete[] trainDataY;
   delete[] trainDataWeights;
   delete[] valDataX;
   delete[] valDataY;
   delete[] valDataWeights;
}

void MethodPyKeras::TestClassification() {
    MethodBase::TestClassification();
}

Double_t MethodPyKeras::GetMvaValue(Double_t *errLower, Double_t *errUpper) {
   // Cannot determine error
   NoErrorCalc(errLower, errUpper);

   // Check whether the model is setup
   // NOTE: unfortunately this is needed because during evaluation ProcessOptions is not called again
   if (!fModelIsSetup) {
      // Setup the trained model
      SetupKerasModel(true);
   }

   // Get signal probability (called mvaValue here)
   const TMVA::Event* e = GetEvent();
   for (UInt_t i=0; i<fNVars; i++) fVals[i] = e->GetValue(i);
   PyRunString("for i,p in enumerate(model.predict(vals)): output[i]=p\n",
               "Failed to get predictions");

   return fOutput[TMVA::Types::kSignal];
}

std::vector<Double_t> MethodPyKeras::GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress) {
   // Check whether the model is setup
   // NOTE: Unfortunately this is needed because during evaluation ProcessOptions is not called again
   if (!fModelIsSetup) {
      // Setup the trained model
      SetupKerasModel(true);
   }

   // Load data to numpy array
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

   float* data = new float[nEvents*fNVars];
   for (UInt_t i=0; i<nEvents; i++) {
      Data()->SetCurrentEvent(i);
      const TMVA::Event *e = GetEvent();
      for (UInt_t j=0; j<fNVars; j++) {
         data[j + i*fNVars] = e->GetValue(j);
      }
   }

   npy_intp dimsData[2] = {(npy_intp)nEvents, (npy_intp)fNVars};
   PyArrayObject* pDataMvaValues = (PyArrayObject*)PyArray_SimpleNewFromData(2, dimsData, NPY_FLOAT, (void*)data);
   if (pDataMvaValues==0) Log() << "Failed to load data to Python array" << Endl;

   // Get prediction for all events
   PyObject* pModel = PyDict_GetItemString(fLocalNS, "model");
   if (pModel==0) Log() << kFATAL << "Failed to get model Python object" << Endl;
   PyArrayObject* pPredictions = (PyArrayObject*) PyObject_CallMethod(pModel, (char*)"predict", (char*)"O", pDataMvaValues);
   if (pPredictions==0) Log() << kFATAL << "Failed to get predictions" << Endl;
   delete[] data;

   // Load predictions to double vector
   // NOTE: The signal probability is given at the output
   std::vector<double> mvaValues(nEvents);
   float* predictionsData = (float*) PyArray_DATA(pPredictions);
   for (UInt_t i=0; i<nEvents; i++) {
      mvaValues[i] = (double) predictionsData[i*fNOutputs + TMVA::Types::kSignal];
   }

   if (logProgress) {
      Log() << kINFO
            << "Elapsed time for evaluation of " << nEvents <<  " events: "
            << timer.GetElapsedTime() << "       " << Endl;
   }


   return mvaValues;
}

std::vector<Float_t>& MethodPyKeras::GetRegressionValues() {
   // Check whether the model is setup
   // NOTE: unfortunately this is needed because during evaluation ProcessOptions is not called again
   if (!fModelIsSetup){
      // Setup the model and load weights
      SetupKerasModel(true);
   }

   // Get regression values
   const TMVA::Event* e = GetEvent();
   for (UInt_t i=0; i<fNVars; i++) fVals[i] = e->GetValue(i);
   PyRunString("for i,p in enumerate(model.predict(vals)): output[i]=p\n",
               "Failed to get predictions");

   // Use inverse transformation of targets to get final regression values
   Event * eTrans = new Event(*e);
   for (UInt_t i=0; i<fNOutputs; ++i) {
      eTrans->SetTarget(i,fOutput[i]);
   }

   const Event* eTrans2 = GetTransformationHandler().InverseTransform(eTrans);
   for (UInt_t i=0; i<fNOutputs; ++i) {
      fOutput[i] = eTrans2->GetTarget(i);
   }

   return fOutput;
}

std::vector<Float_t>& MethodPyKeras::GetMulticlassValues() {
   // Check whether the model is setup
   // NOTE: unfortunately this is needed because during evaluation ProcessOptions is not called again
   if (!fModelIsSetup){
      // Setup the model and load weights
      SetupKerasModel(true);
   }

   // Get class probabilites
   const TMVA::Event* e = GetEvent();
   for (UInt_t i=0; i<fNVars; i++) fVals[i] = e->GetValue(i);
   PyRunString("for i,p in enumerate(model.predict(vals)): output[i]=p\n",
               "Failed to get predictions");

   return fOutput;
}

void MethodPyKeras::ReadModelFromFile() {
}

void MethodPyKeras::GetHelpMessage() const {
// typical length of text line:
//          "|--------------------------------------------------------------|"
   Log() << Endl;
   Log() << "Keras is a high-level API for the Theano and Tensorflow packages." << Endl;
   Log() << "This method wraps the training and predictions steps of the Keras" << Endl;
   Log() << "Python package for TMVA, so that dataloading, preprocessing and" << Endl;
   Log() << "evaluation can be done within the TMVA system. To use this Keras" << Endl;
   Log() << "interface, you have to generate a model with Keras first. Then," << Endl;
   Log() << "this model can be loaded and trained in TMVA." << Endl;
   Log() << Endl;
}

MethodPyKeras::EBackendType MethodPyKeras::GetKerasBackend()  {
   // get the keras backend

   // in case we use tf.keras backend is tensorflow
   if (UseTFKeras())  return kTensorFlow;

   // check first if using tensorflow backend
   PyRunString("keras_backend_is_set =  keras.backend.backend() == \"tensorflow\"");
   PyObject * keras_backend = PyDict_GetItemString(fLocalNS,"keras_backend_is_set");
   if (keras_backend  != nullptr && keras_backend == Py_True)
      return kTensorFlow;

   PyRunString("keras_backend_is_set =  keras.backend.backend() == \"theano\"");
   keras_backend = PyDict_GetItemString(fLocalNS,"keras_backend_is_set");
   if (keras_backend != nullptr && keras_backend == Py_True)
      return kTheano;

   PyRunString("keras_backend_is_set =  keras.backend.backend() == \"cntk\"");
   keras_backend = PyDict_GetItemString(fLocalNS,"keras_backend_is_set");
   if (keras_backend != nullptr && keras_backend == Py_True)
      return kCNTK;

   return kUndefined;
}

TString MethodPyKeras::GetKerasBackendName()  {
   // get the keras backend name
   EBackendType type = GetKerasBackend();
   if (type == kTensorFlow) return "TensorFlow";
   if (type == kTheano) return "Theano";
   if (type == kCNTK) return "CNTK";
   return "Undefined";
}
