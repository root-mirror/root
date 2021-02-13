# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building PyMVA tests
# @author Stefan Wunsch
############################################################################

project(pymva-tests)

set(Libraries Core MathCore TMVA PyMVA)
include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})

# Look for needed python modules
find_python_module(torch QUIET)
find_python_module(keras QUIET)
find_python_module(theano QUIET)
find_python_module(tensorflow QUIET)
find_python_module(sklearn QUIET)

# Enable tests based on available python modules
if(PY_TORCH_FOUND)
   configure_file(generatePyTorchModelClassification.py generatePyTorchModelClassification.py COPYONLY)
   configure_file(generatePyTorchModelMulticlass.py generatePyTorchModelMulticlass.py COPYONLY)
   configure_file(generatePyTorchModelRegression.py generatePyTorchModelRegression.py COPYONLY)
   # Test PyTorch: Binary classification
   ROOT_EXECUTABLE(testPyTorchClassification testPyTorchClassification.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Torch-Classification COMMAND testPyTorchClassification)

   # Test PyTorch: Regression
   ROOT_EXECUTABLE(testPyTorchRegression testPyTorchRegression.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Torch-Regression COMMAND testPyTorchRegression)

   # Test PyTorch: Multi-class classification
   ROOT_EXECUTABLE(testPyTorchMulticlass testPyTorchMulticlass.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Torch-Multiclass COMMAND testPyTorchMulticlass)
endif(PY_TORCH_FOUND)

if((PY_KERAS_FOUND AND PY_THEANO_FOUND) OR (PY_KERAS_FOUND AND PY_TENSORFLOW_FOUND))
   # Test PyKeras: Binary classification
   ROOT_EXECUTABLE(testPyKerasClassification testPyKerasClassification.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Keras-Classification COMMAND testPyKerasClassification)

   # Test PyKeras: Regression
   ROOT_EXECUTABLE(testPyKerasRegression testPyKerasRegression.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Keras-Regression COMMAND testPyKerasRegression)

   # Test PyKeras: Multi-class classification
   ROOT_EXECUTABLE(testPyKerasMulticlass testPyKerasMulticlass.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-Keras-Multiclass COMMAND testPyKerasMulticlass)
endif((PY_KERAS_FOUND AND PY_THEANO_FOUND) OR (PY_KERAS_FOUND AND PY_TENSORFLOW_FOUND))

if(PY_SKLEARN_FOUND)
   # Test PyRandomForest: Classification
   ROOT_EXECUTABLE(testPyRandomForestClassification testPyRandomForestClassification.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-RandomForest-Classification COMMAND testPyRandomForestClassification)

   # Test PyRandomForest: Multi-class classification
   ROOT_EXECUTABLE(testPyRandomForestMulticlass testPyRandomForestMulticlass.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-RandomForest-Multiclass COMMAND testPyRandomForestMulticlass)

   # Test PyGTB: Classification
   ROOT_EXECUTABLE(testPyGTBClassification testPyGTBClassification.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-GTB-Classification COMMAND testPyGTBClassification)

   # Test PyGTB: Multi-class classification
   ROOT_EXECUTABLE(testPyGTBMulticlass testPyGTBMulticlass.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-GTB-Multiclass COMMAND testPyGTBMulticlass)

   # Test PyAdaBoost: Classification
   ROOT_EXECUTABLE(testPyAdaBoostClassification testPyAdaBoostClassification.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-AdaBoost-Classification COMMAND testPyAdaBoostClassification)

   # Test PyAdaBoost: Multi-class classification
   ROOT_EXECUTABLE(testPyAdaBoostMulticlass testPyAdaBoostMulticlass.C
      LIBRARIES ${Libraries})
   ROOT_ADD_TEST(PyMVA-AdaBoost-Multiclass COMMAND testPyAdaBoostMulticlass)
endif(PY_SKLEARN_FOUND)