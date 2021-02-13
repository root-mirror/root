# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT tmva package
# @author Pere Mato, CERN
############################################################################
if (tmva-cpu)
  set(TMVA_CPU_HEADERS
      TMVA/DNN/Architectures/Cpu/Blas.h
  )
#else()
#  set(installoptions ${installoptions} FILTER "Cpu")
endif()

if(tmva-gpu)
  set(TMVA_GPU_HEADERS
       TMVA/DNN/Architectures/Cuda.h
       #  TMVA/DNN/Architectures/Cuda/CudaBuffers.h
       #  TMVA/DNN/Architectures/Cuda/CudaMatrix.h
       #  TMVA/DNN/Architectures/Cuda/CudaTensor.h
  )

else()
  set(installoptions ${installoptions} FILTER "Cuda")
endif()

# Enable parts dependent RDataFrame
if(dataframe)
    set(TMVA_EXTRA_HEADERS
        TMVA/RTensor.hxx
        TMVA/RTensorUtils.hxx
        TMVA/RStandardScaler.hxx
        TMVA/RReader.hxx
        TMVA/RInferenceUtils.hxx
        TMVA/RBDT.hxx
    )
    set(TMVA_EXTRA_SOURCES
        RBDT.cxx
    )
    list(APPEND TMVA_EXTRA_DEPENDENCIES ROOTDataFrame ROOTVecOps)
endif()

if (imt)
  list(APPEND TMVA_EXTRA_DEPENDENCIES Imt)
endif(imt)

ROOT_STANDARD_LIBRARY_PACKAGE(TMVA
  HEADERS
    TMVA/BDTEventWrapper.h
    TMVA/BinarySearchTree.h
    TMVA/BinarySearchTreeNode.h
    TMVA/BinaryTree.h
    TMVA/CCPruner.h
    TMVA/CCTreeWrapper.h
    TMVA/Classification.h
    TMVA/ClassifierFactory.h
    TMVA/ClassInfo.h
    TMVA/Config.h
    TMVA/Configurable.h
    TMVA/ConvergenceTest.h
    TMVA/CostComplexityPruneTool.h
    TMVA/CrossEntropy.h
    TMVA/CrossValidation.h
    TMVA/CvSplit.h
    TMVA/DataInputHandler.h
    TMVA/DataLoader.h
    TMVA/DataSetFactory.h
    TMVA/DataSet.h
    TMVA/DataSetInfo.h
    TMVA/DataSetManager.h
    TMVA/DecisionTree.h
    TMVA/DecisionTreeNode.h
    TMVA/Envelope.h
    TMVA/Event.h
    TMVA/ExpectedErrorPruneTool.h
    TMVA/Executor.h
    TMVA/Factory.h
    TMVA/FitterBase.h
    TMVA/GeneticAlgorithm.h
    TMVA/GeneticFitter.h
    TMVA/GeneticGenes.h
    TMVA/GeneticPopulation.h
    TMVA/GeneticRange.h
    TMVA/GiniIndex.h
    TMVA/GiniIndexWithLaplace.h
    TMVA/HyperParameterOptimisation.h
    TMVA/IFitterTarget.h
    TMVA/IMethod.h
    TMVA/Interval.h
    TMVA/IPruneTool.h
    TMVA/KDEKernel.h
    TMVA/LDA.h
    TMVA/LogInterval.h
    TMVA/LossFunction.h
    TMVA/MCFitter.h
    TMVA/MethodANNBase.h
    TMVA/MethodBase.h
    TMVA/MethodBayesClassifier.h
    TMVA/MethodBDT.h
    TMVA/MethodBoost.h
    TMVA/MethodCategory.h
    TMVA/MethodCFMlpANN_def.h
    TMVA/MethodCFMlpANN.h
    TMVA/MethodCFMlpANN_Utils.h
    TMVA/MethodCompositeBase.h
    TMVA/MethodCrossValidation.h
    TMVA/MethodCuts.h
    TMVA/MethodDL.h
    TMVA/MethodDNN.h
    TMVA/MethodDT.h
    TMVA/MethodFDA.h
    TMVA/MethodFisher.h
    TMVA/MethodHMatrix.h
    TMVA/MethodKNN.h
    TMVA/MethodLD.h
    TMVA/MethodLikelihood.h
    TMVA/MethodMLP.h
    TMVA/MethodPDEFoam.h
    TMVA/MethodPDERS.h
    TMVA/MethodRuleFit.h
    TMVA/MethodSVM.h
    TMVA/MethodTMlpANN.h
    TMVA/MinuitFitter.h
    TMVA/MinuitWrapper.h
    TMVA/MisClassificationError.h
    TMVA/ModulekNN.h
    TMVA/Monitoring.h
    TMVA/MsgLogger.h
    TMVA/NeuralNet.h
    TMVA/Node.h
    TMVA/NodekNN.h
    TMVA/OptimizeConfigParameters.h
    TMVA/Option.h
    TMVA/OptionMap.h
    TMVA/Pattern.h
    TMVA/PDEFoamCell.h
    TMVA/PDEFoamDecisionTreeDensity.h
    TMVA/PDEFoamDecisionTree.h
    TMVA/PDEFoamDensityBase.h
    TMVA/PDEFoamDiscriminantDensity.h
    TMVA/PDEFoamDiscriminant.h
    TMVA/PDEFoamEventDensity.h
    TMVA/PDEFoamEvent.h
    TMVA/PDEFoam.h
    TMVA/PDEFoamKernelBase.h
    TMVA/PDEFoamKernelGauss.h
    TMVA/PDEFoamKernelLinN.h
    TMVA/PDEFoamKernelTrivial.h
    TMVA/PDEFoamMultiTarget.h
    TMVA/PDEFoamTargetDensity.h
    TMVA/PDEFoamTarget.h
    TMVA/PDEFoamVect.h
    TMVA/PDF.h
    TMVA/QuickMVAProbEstimator.h
    TMVA/Ranking.h
    TMVA/Reader.h
    TMVA/RegressionVariance.h
    TMVA/ResultsClassification.h
    TMVA/Results.h
    TMVA/ResultsMulticlass.h
    TMVA/ResultsRegression.h
    TMVA/ROCCalc.h
    TMVA/ROCCurve.h
    TMVA/RootFinder.h
    TMVA/RuleCut.h
    TMVA/RuleEnsemble.h
    TMVA/RuleFitAPI.h
    TMVA/RuleFit.h
    TMVA/RuleFitParams.h
    TMVA/Rule.h
    TMVA/SdivSqrtSplusB.h
    TMVA/SeparationBase.h
    TMVA/SimulatedAnnealingFitter.h
    TMVA/SimulatedAnnealing.h
    TMVA/SVEvent.h
    TMVA/SVKernelFunction.h
    TMVA/SVKernelMatrix.h
    TMVA/SVWorkingSet.h
    TMVA/TActivationChooser.h
    TMVA/TActivation.h
    TMVA/TActivationIdentity.h
    TMVA/TActivationRadial.h
    TMVA/TActivationReLU.h
    TMVA/TActivationSigmoid.h
    TMVA/TActivationTanh.h
    TMVA/Timer.h
    TMVA/TNeuron.h
    TMVA/TNeuronInputAbs.h
    TMVA/TNeuronInputChooser.h
    TMVA/TNeuronInput.h
    TMVA/TNeuronInputSqSum.h
    TMVA/TNeuronInputSum.h
    TMVA/Tools.h
    TMVA/TrainingHistory.h
    TMVA/TransformationHandler.h
    TMVA/TSpline1.h
    TMVA/TSpline2.h
    TMVA/TSynapse.h
    TMVA/Types.h
    TMVA/VariableDecorrTransform.h
    TMVA/VariableGaussTransform.h
    TMVA/VariableIdentityTransform.h
    TMVA/VariableImportance.h
    TMVA/VariableInfo.h
    TMVA/VariableNormalizeTransform.h
    TMVA/VariablePCATransform.h
    TMVA/VariableRearrangeTransform.h
    TMVA/VariableTransformBase.h
    TMVA/VariableTransform.h
    TMVA/VarTransformHandler.h
    TMVA/Version.h
    TMVA/Volume.h
    TMVA/TreeInference/PythonHelpers.hxx
    TMVA/TreeInference/BranchlessTree.hxx
    TMVA/TreeInference/Forest.hxx
    TMVA/TreeInference/Objectives.hxx
    #  TMVA/DNN/Adadelta.h
    #  TMVA/DNN/Adagrad.h
    #  TMVA/DNN/Adam.h
    #  TMVA/DNN/BatchNormLayer.h
    #  TMVA/DNN/DataLoader.h
    #  TMVA/DNN/DeepNet.h
    #  TMVA/DNN/DenseLayer.h
    #  TMVA/DNN/DLMinimizers.h
    #  TMVA/DNN/Functions.h
    #  TMVA/DNN/GeneralLayer.h
    #  TMVA/DNN/ReshapeLayer.h
    #  TMVA/DNN/RMSProp.h
    #  TMVA/DNN/SGD.h
    #  TMVA/DNN/TensorDataLoader.h
    #  TMVA/DNN/CNN/ContextHandles.h
    #  TMVA/DNN/CNN/ConvLayer.h
    #  TMVA/DNN/CNN/MaxPoolLayer.h
    #  TMVA/DNN/RNN/RNNLayer.h
    #   TMVA/DNN/Architectures/Reference.h
    #  TMVA/DNN/Architectures/Reference/DataLoader.h
    #  TMVA/DNN/Architectures/Reference/TensorDataLoader.h
    #   TMVA/DNN/Architectures/Cpu.h
    #  TMVA/DNN/Architectures/Cpu/CpuBuffer.h
    #  TMVA/DNN/Architectures/Cpu/CpuMatrix.h
    #  TMVA/DNN/Architectures/Cpu/CpuTensor.h
    #  ${TMVA_CPU_HEADERS}
    #  ${TMVA_GPU_HEADERS}
    ${TMVA_EXTRA_HEADERS}
  SOURCES
    src/BDTEventWrapper.cxx
    src/BinarySearchTree.cxx
    src/BinarySearchTreeNode.cxx
    src/BinaryTree.cxx
    src/CCPruner.cxx
    src/CCTreeWrapper.cxx
    src/Classification.cxx
    src/ClassifierFactory.cxx
    src/ClassInfo.cxx
    src/Config.cxx
    src/Configurable.cxx
    src/ConvergenceTest.cxx
    src/CostComplexityPruneTool.cxx
    src/CrossEntropy.cxx
    src/CrossValidation.cxx
    src/CvSplit.cxx
    src/DataInputHandler.cxx
    src/DataLoader.cxx
    src/DataSet.cxx
    src/DataSetFactory.cxx
    src/DataSetInfo.cxx
    src/DataSetManager.cxx
    src/DecisionTree.cxx
    src/DecisionTreeNode.cxx
    src/Envelope.cxx
    src/Event.cxx
    src/ExpectedErrorPruneTool.cxx
    src/Factory.cxx
    src/FitterBase.cxx
    src/GeneticAlgorithm.cxx
    src/GeneticFitter.cxx
    src/GeneticGenes.cxx
    src/GeneticPopulation.cxx
    src/GeneticRange.cxx
    src/GiniIndex.cxx
    src/GiniIndexWithLaplace.cxx
    src/HyperParameterOptimisation.cxx
    src/IFitterTarget.cxx
    src/IMethod.cxx
    src/Interval.cxx
    src/KDEKernel.cxx
    src/LDA.cxx
    src/LogInterval.cxx
    src/LossFunction.cxx
    src/MCFitter.cxx
    src/MethodANNBase.cxx
    src/MethodBase.cxx
    src/MethodBayesClassifier.cxx
    src/MethodBDT.cxx
    src/MethodBoost.cxx
    src/MethodCategory.cxx
    src/MethodCFMlpANN.cxx
    src/MethodCFMlpANN_Utils.cxx
    src/MethodCompositeBase.cxx
    src/MethodCrossValidation.cxx
    src/MethodCuts.cxx
    src/MethodDL.cxx
    src/MethodDNN.cxx
    src/MethodDT.cxx
    src/MethodFDA.cxx
    src/MethodFisher.cxx
    src/MethodHMatrix.cxx
    src/MethodKNN.cxx
    src/MethodLD.cxx
    src/MethodLikelihood.cxx
    src/MethodMLP.cxx
    src/MethodPDEFoam.cxx
    src/MethodPDERS.cxx
    src/MethodPlugins.cxx
    src/MethodRuleFit.cxx
    src/MethodSVM.cxx
    src/MethodTMlpANN.cxx
    src/MinuitFitter.cxx
    src/MinuitWrapper.cxx
    src/MisClassificationError.cxx
    src/ModulekNN.cxx
    src/MsgLogger.cxx
    src/NeuralNet.cxx
    src/Node.cxx
    src/OptimizeConfigParameters.cxx
    src/Option.cxx
    src/OptionMap.cxx
    src/PDEFoamCell.cxx
    src/PDEFoam.cxx
    src/PDEFoamDecisionTree.cxx
    src/PDEFoamDecisionTreeDensity.cxx
    src/PDEFoamDensityBase.cxx
    src/PDEFoamDiscriminant.cxx
    src/PDEFoamDiscriminantDensity.cxx
    src/PDEFoamEvent.cxx
    src/PDEFoamEventDensity.cxx
    src/PDEFoamKernelBase.cxx
    src/PDEFoamKernelGauss.cxx
    src/PDEFoamKernelLinN.cxx
    src/PDEFoamKernelTrivial.cxx
    src/PDEFoamMultiTarget.cxx
    src/PDEFoamTarget.cxx
    src/PDEFoamTargetDensity.cxx
    src/PDEFoamVect.cxx
    src/PDF.cxx
    src/QuickMVAProbEstimator.cxx
    src/Ranking.cxx
    src/Reader.cxx
    src/RegressionVariance.cxx
    src/ResultsClassification.cxx
    src/Results.cxx
    src/ResultsMulticlass.cxx
    src/ResultsRegression.cxx
    src/ROCCalc.cxx
    src/ROCCurve.cxx
    src/RootFinder.cxx
    src/RuleCut.cxx
    src/Rule.cxx
    src/RuleEnsemble.cxx
    src/RuleFitAPI.cxx
    src/RuleFit.cxx
    src/RuleFitParams.cxx
    src/SdivSqrtSplusB.cxx
    src/SeparationBase.cxx
    src/SimulatedAnnealing.cxx
    src/SimulatedAnnealingFitter.cxx
    src/SVEvent.cxx
    src/SVKernelFunction.cxx
    src/SVKernelMatrix.cxx
    src/SVWorkingSet.cxx
    src/TActivationChooser.cxx
    src/TActivation.cxx
    src/TActivationIdentity.cxx
    src/TActivationRadial.cxx
    src/TActivationReLU.cxx
    src/TActivationSigmoid.cxx
    src/TActivationTanh.cxx
    src/Timer.cxx
    src/TNeuron.cxx
    src/TNeuronInputAbs.cxx
    src/TNeuronInputChooser.cxx
    src/TNeuronInput.cxx
    src/TNeuronInputSqSum.cxx
    src/TNeuronInputSum.cxx
    src/Tools.cxx
    src/TrainingHistory.cxx
    src/TransformationHandler.cxx
    src/TSpline1.cxx
    src/TSpline2.cxx
    src/TSynapse.cxx
    src/Types.cxx
    src/VariableDecorrTransform.cxx
    src/VariableGaussTransform.cxx
    src/VariableIdentityTransform.cxx
    src/VariableImportance.cxx
    src/VariableInfo.cxx
    src/VariableNormalizeTransform.cxx
    src/VariablePCATransform.cxx
    src/VariableRearrangeTransform.cxx
    src/VariableTransformBase.cxx
    src/VariableTransform.cxx
    src/VarTransformHandler.cxx
    src/Volume.cxx
    src/DNN/Architectures/Reference.cxx
    src/DNN/Architectures/Reference/DataLoader.cxx
    src/DNN/Architectures/Reference/TensorDataLoader.cxx
    src/DNN/Architectures/Cpu.cxx
    src/DNN/Architectures/Cpu/CpuBuffer.cxx
    src/DNN/Architectures/Cpu/CpuMatrix.cxx
    ${TMVA_EXTRA_SOURCES}
  DEPENDENCIES
    TreePlayer
    Tree
    Hist
    Matrix
    Minuit
    MLP
    MathCore
    Core
    RIO
    XMLIO
    ${TMVA_EXTRA_DEPENDENCIES}
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  INSTALL_OPTIONS
    ${installoptions}
)

if(MSVC)
  target_compile_definitions(TMVA PRIVATE _USE_MATH_DEFINES)
endif()

if(vdt OR builtin_vdt)
  target_include_directories(TMVA PRIVATE ${VDT_INCLUDE_DIRS})
endif()

if(tmva-cpu)
  target_include_directories(TMVA SYSTEM PRIVATE ${TBB_INCLUDE_DIRS})
  target_link_libraries(TMVA PRIVATE ${TBB_LIBRARIES})
  set_target_properties(TMVA PROPERTIES COMPILE_FLAGS "${TBB_CXXFLAGS}")

  if(BLAS_FOUND)
    target_link_libraries(TMVA PRIVATE ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
  elseif(GSL_FOUND)
    target_compile_definitions(TMVA PRIVATE -DDNN_USE_CBLAS)
    target_include_directories(TMVA SYSTEM PRIVATE ${GSL_INCLUDE_DIR})
    target_link_libraries(TMVA PRIVATE ${GSL_CBLAS_LIBRARY})
    if(builtin_gsl)
      add_dependencies(TMVA GSL)
    endif()
  else()
    message(FATAL_ERROR "tmva-cpu enabled but neither BLAS nor GSL BLAS were found")
  endif()
endif()

if(tmva-gpu)
  if (tmva-cudnn)
    message(STATUS "Using Cuda+cuDNN for TMVA Deep Learning on GPU")
    target_sources(TMVA PRIVATE
                   src/DNN/Architectures/Cuda.cu
                   src/DNN/Architectures/Cuda/CudaBuffers.cxx
                   src/DNN/Architectures/Cuda/CudaMatrix.cu
                   src/DNN/Architectures/Cuda/CudaTensor.cu
                   src/DNN/Architectures/Cudnn/TensorDataLoader.cxx
                   src/DNN/Architectures/Cudnn.cu)
    target_include_directories(TMVA PRIVATE ${CUDA_INCLUDE_DIRS} ${CUDNN_INCLUDE_DIRS} )
    target_link_libraries(TMVA PRIVATE ${CUDA_CUBLAS_LIBRARIES} ${CUDNN_LIBRARIES})
  else()
    message(STATUS "cuDNN not found or disabled - use only Cuda+Cublas for TMVA Deep Learning on GPU")
    target_include_directories(TMVA PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(TMVA PRIVATE ${CUDA_CUBLAS_LIBRARIES})
    target_sources(TMVA PRIVATE
                   src/DNN/Architectures/Cuda.cu
                   src/DNN/Architectures/Cuda/CudaBuffers.cxx
                   src/DNN/Architectures/Cuda/CudaMatrix.cu
                   src/DNN/Architectures/Cuda/CudaTensor.cu )
   endif()
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
ROOT_ADD_TEST_SUBDIRECTORY(test/crossvalidation)
ROOT_ADD_TEST_SUBDIRECTORY(test/DNN)
ROOT_ADD_TEST_SUBDIRECTORY(test/Method)
ROOT_ADD_TEST_SUBDIRECTORY(test/ROC)
ROOT_ADD_TEST_SUBDIRECTORY(test/envelope)
ROOT_ADD_TEST_SUBDIRECTORY(test/DNN/CNN)
ROOT_ADD_TEST_SUBDIRECTORY(test/DNN/RNN)
ROOT_ADD_TEST_SUBDIRECTORY(test/DNN/LSTM)
ROOT_ADD_TEST_SUBDIRECTORY(test/DNN/GRU)
