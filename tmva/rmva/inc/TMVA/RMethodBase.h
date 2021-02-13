// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata,Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodBase                                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method based on ROOTR                      *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodBase
#define ROOT_TMVA_RMethodBase

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodBase                                                          //
//                                                                      //
// Virtual base class for all TMVA method based on ROOTR                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/MethodBase.h"

#include <TRInterface.h>

#include <vector>
#include <string>

class TGraph;
class TTree;
class TDirectory;
class TSpline;
class TH1F;
class TH1D;

namespace TMVA {

   class Ranking;
   class PDF;
   class TSpline1;
   class MethodCuts;
   class MethodBoost;
   class DataSetInfo;

   class RMethodBase : public MethodBase {

      friend class Factory;
   protected:
      ROOT::R::TRInterface &r;
   public:

      // default constructur
      RMethodBase(const TString &jobName,
                  Types::EMVA methodType,
                  const TString &methodTitle,
                  DataSetInfo &dsi,
                  const TString &theOption = "", ROOT::R::TRInterface &_r = ROOT::R::TRInterface::Instance());

      // constructor used for Testing + Application of the MVA, only (no training),
      // using given weight file
      RMethodBase(Types::EMVA methodType,
                  DataSetInfo &dsi,
                  const TString &weightFile, ROOT::R::TRInterface &_r = ROOT::R::TRInterface::Instance());

      // default destructur
      virtual ~RMethodBase() {};
      virtual void     Train() = 0;
      // options treatment
      virtual void     Init()           = 0;
      virtual void     DeclareOptions() = 0;
      virtual void     ProcessOptions() = 0;
      // create ranking
      virtual const Ranking *CreateRanking() = 0;

      virtual Double_t GetMvaValue(Double_t *errLower = 0, Double_t *errUpper = 0) = 0;

      Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets) = 0;
   protected:
      // the actual "weights"
      virtual void AddWeightsXMLTo(void *parent) const = 0;
      virtual void ReadWeightsFromXML(void *wghtnode) = 0;
      virtual void ReadWeightsFromStream(std::istream &) = 0;        // backward compatibility
      virtual void ReadWeightsFromStream(TFile &) {}                 // backward compatibility


      void LoadData();//Read data from Data() Aand DataInfo() to Dataframes and Vectors
   protected:
      ROOT::R::TRDataFrame fDfTrain;//signal and backgrd
      ROOT::R::TRDataFrame fDfTest;
      TVectorD             fWeightTrain;
      TVectorD             fWeightTest;
      std::vector<std::string> fFactorTrain;
      std::vector<std::string> fFactorTest;
      ROOT::R::TRDataFrame fDfSpectators;

   private:
      ClassDef(RMethodBase, 0) // Virtual base class for all TMVA method

   };
} // namespace TMVA

#endif


