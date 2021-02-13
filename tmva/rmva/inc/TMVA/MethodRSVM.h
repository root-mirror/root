// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata,Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodRSVM                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package RSVM  method based on ROOTR                                    *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodRSVM
#define ROOT_TMVA_RMethodRSVM

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodRSVM                                                          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/RMethodBase.h"
#include <vector>

namespace TMVA {

   class Factory;  // DSMTEST
   class Reader;   // DSMTEST
   class DataSetManager;  // DSMTEST
   class Types;
   class MethodRSVM : public RMethodBase {

   public :

      // constructors
      MethodRSVM(const TString &jobName,
                 const TString &methodTitle,
                 DataSetInfo &theData,
                 const TString &theOption = "");

      MethodRSVM(DataSetInfo &dsi,
                 const TString &theWeightFile);


      ~MethodRSVM(void);
      void     Train();
      // options treatment
      void     Init();
      void     DeclareOptions();
      void     ProcessOptions();
      // create ranking
      const Ranking *CreateRanking()
      {
         return NULL;  // = 0;
      }


      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

      // performs classifier testing
      virtual void     TestClassification();


      Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0);

      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo(void * /*parent*/) const {}  // = 0;
      virtual void ReadWeightsFromXML(void * /*wghtnode*/) {} // = 0;
      virtual void ReadWeightsFromStream(std::istream &) {} //= 0;       // backward compatibility
      void ReadModelFromFile();

      // signal/background classification response for all current set of data 
      virtual std::vector<Double_t> GetMvaValues(Long64_t firstEvt = 0, Long64_t lastEvt = -1, Bool_t logProgress = false);

   private :
      DataSetManager    *fDataSetManager;     // DSMTEST
      friend class Factory;                   // DSMTEST
      friend class Reader;                    // DSMTEST
   protected:
      UInt_t fMvaCounter;
      std::vector<Float_t> fProbResultForTrainSig;
      std::vector<Float_t> fProbResultForTestSig;

      //Booking options
      Bool_t fScale;//A logical vector indicating the variables to be scaled. If
      //‘scale’ is of length 1, the value is recycled as many times
      //as needed.  Per default, data are scaled internally (both ‘x’
      //and ‘y’ variables) to zero mean and unit variance. The center
      //and scale values are returned and used for later predictions.
      TString fType;//‘svm’ can be used as a classification machine, as a
      //regression machine, or for novelty detection.  Depending of
      //whether ‘y’ is a factor or not, the default setting for
      //‘type’ is ‘C-classification’ or ‘eps-regression’,
      //respectively, but may be overwritten by setting an explicit value.
      //Valid options are:
      // - ‘C-classification’
      // - ‘nu-classification’
      // - ‘one-classification’ (for novelty detection)
      // - ‘eps-regression’
      // - ‘nu-regression’
      TString fKernel;//the kernel used in training and predicting. You might
      //consider changing some of the following parameters, depending on the kernel type.
      //linear: u'*v
      //polynomial: (gamma*u'*v + coef0)^degree
      //radial basis: exp(-gamma*|u-v|^2)
      //sigmoid: tanh(gamma*u'*v + coef0)
      Int_t fDegree;//parameter needed for kernel of type ‘polynomial’ (default: 3)
      Float_t fGamma;//parameter needed for all kernels except ‘linear’ (default: 1/(data dimension))
      Float_t fCoef0;//parameter needed for kernels of type ‘polynomial’ and ‘sigmoid’ (default: 0)
      Float_t fCost;//cost of constraints violation (default: 1)-it is the
      //‘C’-constant of the regularization term in the Lagrange formulation.
      Float_t fNu;//parameter needed for ‘nu-classification’, ‘nu-regression’, and ‘one-classification’
      Float_t fCacheSize;//cache memory in MB (default 40)
      Float_t fTolerance;//tolerance of termination criterion (default: 0.001)
      Float_t fEpsilon;//epsilon in the insensitive-loss function (default: 0.1)
      Bool_t fShrinking;//option whether to use the shrinking-heuristics (default: ‘TRUE’)
      Float_t fCross;//if a integer value k>0 is specified, a k-fold cross
      //validation on the training data is performed to assess the
      //quality of the model: the accuracy rate for classification
      //and the Mean Squared Error for regression
      Bool_t fProbability;//logical indicating whether the model should allow for probability predictions.
      Bool_t fFitted;//logical indicating whether the fitted values should be computed and included in the model or not (default: ‘TRUE’)

      static Bool_t IsModuleLoaded;
      ROOT::R::TRFunctionImport svm;
      ROOT::R::TRFunctionImport predict;
      ROOT::R::TRFunctionImport asfactor;
      ROOT::R::TRObject *fModel;
      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRSVM, 0)
   };
} // namespace TMVA
#endif
