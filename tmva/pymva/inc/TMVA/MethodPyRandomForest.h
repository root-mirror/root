// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyRandomForest                                                  *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      scikit-learn Package RandomForestClassifier  method based on python       *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodPyRandomForest
#define ROOT_TMVA_MethodPyRandomForest

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodPyRandomForest                                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/PyMethodBase.h"
#include <vector>

namespace TMVA {

   class Factory;  // DSMTEST
   class Reader;   // DSMTEST
   class DataSetManager;  // DSMTEST
   class Types;
   class MethodPyRandomForest : public PyMethodBase {

   public :
      // constructors
      MethodPyRandomForest(const TString &jobName,
                           const TString &methodTitle,
                           DataSetInfo &theData,
                           const TString &theOption = "");

      MethodPyRandomForest(DataSetInfo &dsi,
                           const TString &theWeightFile);

      ~MethodPyRandomForest(void);
      void Train();

      // options treatment
      void Init();
      void DeclareOptions();
      void ProcessOptions();

      // create ranking
      const Ranking *CreateRanking();

      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

      // performs classifier testing
      virtual void TestClassification();

      // Get class probabilities of given event
      Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0);
      std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);
      std::vector<Float_t>& GetMulticlassValues();

      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /* parent */) const {} // = 0;
      virtual void ReadWeightsFromXML(void * /* wghtnode */) {} // = 0;
      virtual void ReadWeightsFromStream(std::istream &) {} //= 0; // backward compatibility

      void ReadModelFromFile();

   private :
      DataSetManager *fDataSetManager;     // DSMTEST
      friend class Factory;                   // DSMTEST
      friend class Reader;                    // DSMTEST

   protected:
      std::vector<Double_t> mvaValues;
      std::vector<Float_t> classValues;

      UInt_t fNvars; // number of variables
      UInt_t fNoutputs; // number of outputs
      TString fFilenameClassifier; // Path to serialized classifier (default in `weights` folder)

      // RandomForest options

      PyObject* pNestimators;
      Int_t fNestimators; //integer, optional (default=10)
      //The number of trees in the forest.

      PyObject* pCriterion;
      TString fCriterion; //string, optional (default="gini")
      //The function to measure the quality of a split. Supported criteria are
      //"gini" for the Gini impurity and "entropy" for the information gain.
      //Note: this parameter is tree-specific.

      PyObject* pMaxDepth;
      TString fMaxDepth; //integer or None, optional (default=None)
      //The maximum depth of the tree. If None, then nodes are expanded until
      //all leaves are pure or until all leaves contain less than `fMinSamplesSplit`.

      PyObject* pMinSamplesSplit;
      Int_t fMinSamplesSplit; //integer, optional (default=2)
      //The minimum number of samples required to split an internal node.

      PyObject* pMinSamplesLeaf;
      Int_t fMinSamplesLeaf; //integer, optional (default=1)
      //The minimum number of samples in newly created leaves.  A split is
      //discarded if after the split, one of the leaves would contain less then
      //``min_samples_leaf`` samples.
      //Note: this parameter is tree-specific.

      PyObject* pMinWeightFractionLeaf;
      Double_t fMinWeightFractionLeaf; //float, optional (default=0.)
      //The minimum weighted fraction of the input samples required to be at a
      //leaf node.
      //Note: this parameter is tree-specific.

      PyObject* pMaxFeatures;
      TString fMaxFeatures; //int, float, string or None, optional (default="auto")
      //The number of features to consider when looking for the best split:
      //- If int, then consider `max_features` features at each split.
      //- If float, then `max_features` is a percentage and
      //`int(max_features * n_features)` features are considered at each split.
      //- If "auto", then `max_features=sqrt(n_features)`.
      //- If "sqrt", then `max_features=sqrt(n_features)`.
      //- If "log2", then `max_features=log2(n_features)`.
      //- If None, then `max_features=n_features`.
      //     Note: the search for a split does not stop until at least one
      //     valid partition of the node samples is found, even if it requires to
      //     effectively inspect more than ``max_features`` features.
      //     Note: this parameter is tree-specific.

      PyObject* pMaxLeafNodes;
      TString fMaxLeafNodes; //int or None, optional (default=None)
      //Grow trees with ``max_leaf_nodes`` in best-first fashion.
      //Best nodes are defined as relative reduction in impurity.
      //If None then unlimited number of leaf nodes.
      //If not None then ``max_depth`` will be ignored.

      PyObject* pBootstrap;
      Bool_t fBootstrap; //boolean, optional (default=True)
      //Whether bootstrap samples are used when building trees.

      PyObject* pOobScore;
      Bool_t fOobScore; //Whether to use out-of-bag samples to estimate
      //the generalization error.

      PyObject* pNjobs;
      Int_t fNjobs; // integer, optional (default=1)
      //The number of jobs to run in parallel for both `fit` and `predict`.
      //If -1, then the number of jobs is set to the number of cores.

      PyObject* pRandomState;
      TString fRandomState; //int, RandomState instance or None, optional (default=None)
      //If int, random_state is the seed used by the random number generator;
      //If RandomState instance, random_state is the random number generator;
      //If None, the random number generator is the RandomState instance used
      //by `np.random`.

      PyObject* pVerbose;
      Int_t fVerbose; //Controls the verbosity of the tree building process.

      PyObject* pWarmStart;
      Bool_t fWarmStart; //bool, optional (default=False)
      //When set to ``True``, reuse the solution of the previous call to fit
      //and add more estimators to the ensemble, otherwise, just fit a whole
      //new forest.

      PyObject* pClassWeight;
      TString fClassWeight; //dict, list of dicts, "auto", "subsample" or None, optional
      //Weights associated with classes in the form ``{class_label: weight}``.
      //If not given, all classes are supposed to have weight one. For
      //multi-output problems, a list of dicts can be provided in the same
      //order as the columns of y.
      //The "auto" mode uses the values of y to automatically adjust
      //weights inversely proportional to class frequencies in the input data.
      //The "subsample" mode is the same as "auto" except that weights are
      //computed based on the bootstrap sample for every tree grown.
      //For multi-output, the weights of each column of y will be multiplied.
      //Note that these weights will be multiplied with sample_weight (passed
      //through the fit method) if sample_weight is specified.

      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodPyRandomForest, 0)
   };

} // namespace TMVA

#endif // ROOT_TMVA_MethodPyRandomForest
