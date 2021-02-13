// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPyAdaBoost                                                      *
 * Web    : http://oproject.org                                                   *
 *                                                                                *
 * Description:                                                                   *
 *      scikit-learn package AdaBoostClassifier method based on python            *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodPyAdaBoost
#define ROOT_TMVA_MethodPyAdaBoost

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodPyAdaBoost                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/PyMethodBase.h"

#include "TString.h"
#include <vector>

namespace TMVA {

   class Factory;
   class Reader;
   class DataSetManager;
   class Types;
   class MethodPyAdaBoost : public PyMethodBase {

   public :
      MethodPyAdaBoost(const TString &jobName,
                       const TString &methodTitle,
                       DataSetInfo &theData,
                       const TString &theOption = "");

      MethodPyAdaBoost(DataSetInfo &dsi,
                       const TString &theWeightFile);

      ~MethodPyAdaBoost();

      void Train();

      void Init();
      void DeclareOptions();
      void ProcessOptions();

      // create ranking
      const Ranking *CreateRanking();

      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

      // performs classifier testing
      virtual void TestClassification();

      Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0);
      std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);
      std::vector<Float_t>& GetMulticlassValues();

      virtual void ReadModelFromFile();

      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /*parent */ ) const {} // = 0;
      virtual void ReadWeightsFromXML(void * /*wghtnode*/ ) {} // = 0;
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

      //AdaBoost options

      PyObject* pBaseEstimator;
      TString fBaseEstimator; //object, optional (default=DecisionTreeClassifier)
      //The base estimator from which the boosted ensemble is built.
      //Support for sample weighting is required, as well as proper `classes_`
      //and `n_classes_` attributes.

      PyObject* pNestimators;
      Int_t fNestimators; //integer, optional (default=10)
      //The number of trees in the forest.

      PyObject* pLearningRate;
      Double_t fLearningRate; //loat, optional (default=1.)
      //Learning rate shrinks the contribution of each classifier by
      //``learning_rate``. There is a trade-off between ``learning_rate`` and ``n_estimators``.

      PyObject* pAlgorithm;
      TString fAlgorithm; //{'SAMME', 'SAMME.R'}, optional (default='SAMME.R')
      //If 'SAMME.R' then use the SAMME.R real boosting algorithm.
      //``base_estimator`` must support calculation of class probabilities.
      //If 'SAMME' then use the SAMME discrete boosting algorithm.
      //The SAMME.R algorithm typically converges faster than SAMME,
      //achieving a lower test error with fewer boosting iterations.

      PyObject* pRandomState;
      TString fRandomState; //int, RandomState instance or None, optional (default=None)
      //If int, random_state is the seed used by the random number generator;
      //If RandomState instance, random_state is the random number generator;
      //If None, the random number generator is the RandomState instance used by `np.random`.

      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodPyAdaBoost, 0)
   };

} // namespace TMVA

#endif // ROOT_TMVA_MethodPyAdaBoost
