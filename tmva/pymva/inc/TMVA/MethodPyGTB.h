// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyGTB                                                           *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      scikit-learn Package GradientBoostingClassifier  method based on python   *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodPyGTB
#define ROOT_TMVA_MethodPyGTB

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodPyGTB                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/PyMethodBase.h"
#include <vector>

namespace TMVA {

   class Factory;
   class Reader;
   class DataSetManager;
   class Types;
   class MethodPyGTB : public PyMethodBase {

   public :
      MethodPyGTB(const TString &jobName,
                  const TString &methodTitle,
                  DataSetInfo &theData,
                  const TString &theOption = "");
      MethodPyGTB(DataSetInfo &dsi,
                  const TString &theWeightFile);
      ~MethodPyGTB(void);

      void Train();
      void Init();
      void DeclareOptions();
      void ProcessOptions();

      const Ranking *CreateRanking();

      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

      virtual void TestClassification();

      Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0);
      std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);
      std::vector<Float_t>& GetMulticlassValues();

      virtual void ReadModelFromFile();

      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /* parent */ ) const {} // = 0;
      virtual void ReadWeightsFromXML(void * /*wghtnode*/) {} // = 0;
      virtual void ReadWeightsFromStream(std::istream &) {} //= 0; backward compatibility

   private :
      DataSetManager *fDataSetManager;
      friend class Factory;
      friend class Reader;

   protected:
      std::vector<Double_t> mvaValues;
      std::vector<Float_t> classValues;

      UInt_t fNvars; // number of variables
      UInt_t fNoutputs; // number of outputs
      TString fFilenameClassifier; // Path to serialized classifier (default in `weights` folder)

      //GTB options

      PyObject* pLoss;
      TString fLoss; // {'deviance', 'exponential'}, optional (default='deviance')
      //loss function to be optimized. 'deviance' refers to
      //deviance (= logistic regression) for classification
      //with probabilistic outputs. For loss 'exponential' gradient
      //boosting recovers the AdaBoost algorithm.

      PyObject* pLearningRate;
      Double_t fLearningRate; //float, optional (default=0.1)
      //learning rate shrinks the contribution of each tree by `learning_rate`.
      //There is a trade-off between learning_rate and n_estimators.

      PyObject* pNestimators;
      Int_t fNestimators; //integer, optional (default=10)
      //The number of trees in the forest.

      PyObject* pSubsample;
      Double_t fSubsample; //float, optional (default=1.0)
      //The fraction of samples to be used for fitting the individual base
      //learners. If smaller than 1.0 this results in Stochastic Gradient
      //Boosting. `subsample` interacts with the parameter `n_estimators`.
      //Choosing `subsample < 1.0` leads to a reduction of variance
      //and an increase in bias.

      PyObject* pMinSamplesSplit;
      Int_t fMinSamplesSplit; // integer, optional (default=2)
      //The minimum number of samples required to split an internal node.

      PyObject* pMinSamplesLeaf;
      Int_t fMinSamplesLeaf; //integer, optional (default=1)
      //The minimum number of samples required to be at a leaf node.

      PyObject* pMinWeightFractionLeaf;
      Double_t fMinWeightFractionLeaf; //float, optional (default=0.)
      //The minimum weighted fraction of the input samples required to be at a leaf node.

      PyObject* pMaxDepth;
      Int_t fMaxDepth; //integer, optional (default=3)
      //maximum depth of the individual regression estimators. The maximum
      //depth limits the number of nodes in the tree. Tune this parameter
      //for best performance; the best value depends on the interaction
      //of the input variables.
      //Ignored if ``max_leaf_nodes`` is not None.

      PyObject* pInit;
      TString fInit; //BaseEstimator, None, optional (default=None)
      //An estimator object that is used to compute the initial
      //predictions. ``init`` has to provide ``fit`` and ``predict``.
      //If None it uses ``loss.init_estimator``.

      PyObject* pRandomState;
      TString fRandomState; //int, RandomState instance or None, optional (default=None)
      //If int, random_state is the seed used by the random number generator;
      //If RandomState instance, random_state is the random number generator;
      //If None, the random number generator is the RandomState instance used
      //by `np.random`.

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

      PyObject* pVerbose;
      Int_t fVerbose; //Controls the verbosity of the tree building process.

      PyObject* pMaxLeafNodes;
      TString fMaxLeafNodes; //int or None, optional (default=None)
      //Grow trees with ``max_leaf_nodes`` in best-first fashion.
      //Best nodes are defined as relative reduction in impurity.
      //If None then unlimited number of leaf nodes.
      //If not None then ``max_depth`` will be ignored.

      PyObject* pWarmStart;
      Bool_t fWarmStart; //bool, optional (default=False)
      //When set to ``True``, reuse the solution of the previous call to fit
      //and add more estimators to the ensemble, otherwise, just fit a whole
      //new forest.

      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodPyGTB, 0)
   };

} // namespace TMVA

#endif // ROOT_TMVA_PyMethodGTB
