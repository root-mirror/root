#ifdef __CINT__

#include "RConfigure.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ nestedclass;


#pragma link C++ class TMVA::Configurable+;
#pragma link C++ class TMVA::Event+;
#pragma link C++ class TMVA::kNN::Event+;
#pragma link C++ class TMVA::Factory+;

#pragma link C++ class TMVA::Envelope+;
#pragma link C++ class TMVA::OptionMap+;
#pragma link C++ class TMVA::VariableImportance+;
#pragma link C++ class TMVA::CrossValidation+;
#pragma link C++ class TMVA::CrossValidationFoldResult+;
#pragma link C++ class TMVA::CvSplit+;
#pragma link C++ class TMVA::CvSplitKFolds + ;
#pragma link C++ class TMVA::HyperParameterOptimisation+;

#pragma link C++ class TMVA::Experimental::Classification + ;
#pragma link C++ class TMVA::Experimental::ClassificationResult + ;

//required to enable serialization on DataLoader for paralellism.
#pragma link C++ class TMVA::OptionBase+;
#pragma link C++ class TMVA::Results+;
#pragma link C++ class TMVA::ResultsClassification+;
#pragma link C++ class TMVA::ResultsMulticlass+;
#pragma link C++ class TMVA::ResultsRegression+;
#pragma link C++ class TMVA::DataLoader+;
#pragma link C++ class TMVA::TreeInfo+;
#pragma link C++ class TMVA::VariableInfo+;
#pragma link C++ class TMVA::ClassInfo+;
#pragma link C++ class TMVA::DataInputHandler+;
#pragma link C++ class TMVA::DataSet+;
#pragma link C++ class TMVA::DataSetInfo+;
#pragma link C++ class TMVA::DataSetManager+;
#pragma link C++ class TMVA::DataSetFactory+;

// the classifiers
#pragma link C++ class TMVA::MethodBase+;
#pragma link C++ class TMVA::MethodCompositeBase+;
#pragma link C++ class TMVA::MethodANNBase+;
#pragma link C++ class TMVA::MethodTMlpANN+;
#pragma link C++ class TMVA::MethodRuleFit+;
#pragma link C++ class TMVA::MethodCuts+;
#pragma link C++ class TMVA::MethodFisher+;
#pragma link C++ class TMVA::MethodKNN+;
#pragma link C++ class TMVA::MethodCFMlpANN+;
#pragma link C++ class TMVA::MethodCFMlpANN_Utils+;
#pragma link C++ class TMVA::MethodLikelihood+;
#pragma link C++ class TMVA::MethodHMatrix+;
#pragma link C++ class TMVA::MethodPDERS+;
#pragma link C++ class TMVA::MethodBDT+;
#pragma link C++ class TMVA::MethodDT+;
#pragma link C++ class TMVA::MethodSVM+;
#pragma link C++ class TMVA::MethodBayesClassifier+;
#pragma link C++ class TMVA::MethodFDA+;
#pragma link C++ class TMVA::MethodMLP+;
#pragma link C++ class TMVA::MethodBoost+;
#pragma link C++ class TMVA::MethodPDEFoam+;
#pragma link C++ class TMVA::MethodLD+;
#pragma link C++ class TMVA::MethodCategory+;
#pragma link C++ class TMVA::MethodDNN+;
#pragma link C++ class TMVA::MethodCrossValidation+;
#pragma link C++ class TMVA::MethodDL+;

#ifdef R__HAS_DATAFRAME
// BDT inference
#pragma link C++ class TMVA::Experimental::RBDT<TMVA::Experimental::BranchlessForest<float>>;
#pragma link C++ class TMVA::Experimental::RBDT<TMVA::Experimental::BranchlessJittedForest<float>>;
#endif
#endif
