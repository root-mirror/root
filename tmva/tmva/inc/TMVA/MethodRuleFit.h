// @(#)root/tmva $Id$
// Author: Fredrik Tegenfeldt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodRuleFit                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Friedman's RuleFit method                                                 *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U.                                                             *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodRuleFit
#define ROOT_TMVA_MethodRuleFit

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodRuleFit                                                        //
//                                                                      //
// J Friedman's RuleFit method                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/MethodBase.h"
#include "TMatrixDfwd.h"
#include "TVectorD.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/RuleFit.h"
#include <vector>

namespace TMVA {

   class SeparationBase;

   class MethodRuleFit : public MethodBase {

   public:

      MethodRuleFit( const TString& jobName,
                     const TString& methodTitle,
                     DataSetInfo& theData,
                     const TString& theOption = "");

      MethodRuleFit( DataSetInfo& theData,
                     const TString& theWeightFile);

      virtual ~MethodRuleFit( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ );

      // training method
      void Train( void );

      using MethodBase::ReadWeightsFromStream;

      // write weights to file
      void AddWeightsXMLTo     ( void* parent ) const;

      // read weights from file
      void ReadWeightsFromStream( std::istream& istr );
      void ReadWeightsFromXML   ( void* wghtnode );

      // calculate the MVA value
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0 );

      // write method specific histos to target file
      void WriteMonitoringHistosToFile( void ) const;

      // ranking of input variables
      const Ranking* CreateRanking();

      Bool_t                                   UseBoost()           const   { return fUseBoost; }

      // accessors
      RuleFit*                                 GetRuleFitPtr()              { return &fRuleFit; }
      const RuleFit*                           GetRuleFitConstPtr() const   { return &fRuleFit; }
      TDirectory*                              GetMethodBaseDir()   const   { return BaseDir(); }
      const std::vector<TMVA::Event*>&         GetTrainingEvents()  const   { return fEventSample; }
      const std::vector<TMVA::DecisionTree*>&  GetForest()          const   { return fForest; }
      Int_t                                    GetNTrees()          const   { return fNTrees; }
      Double_t                                 GetTreeEveFrac()     const   { return fTreeEveFrac; }
      const SeparationBase*                    GetSeparationBaseConst() const { return fSepType; }
      SeparationBase*                          GetSeparationBase()  const   { return fSepType; }
      TMVA::DecisionTree::EPruneMethod         GetPruneMethod()     const   { return fPruneMethod; }
      Double_t                                 GetPruneStrength()   const   { return fPruneStrength; }
      Double_t                                 GetMinFracNEve()     const   { return fMinFracNEve; }
      Double_t                                 GetMaxFracNEve()     const   { return fMaxFracNEve; }
      Int_t                                    GetNCuts()           const   { return fNCuts; }
      //
      Int_t                                    GetGDNPathSteps()    const   { return fGDNPathSteps; }
      Double_t                                 GetGDPathStep()      const   { return fGDPathStep; }
      Double_t                                 GetGDErrScale()      const   { return fGDErrScale; }
      Double_t                                 GetGDPathEveFrac()   const   { return fGDPathEveFrac; }
      Double_t                                 GetGDValidEveFrac()  const   { return fGDValidEveFrac; }
      //
      Double_t                                 GetLinQuantile()     const   { return fLinQuantile; }

      const TString                            GetRFWorkDir()       const   { return fRFWorkDir; }
      Int_t                                    GetRFNrules()        const   { return fRFNrules; }
      Int_t                                    GetRFNendnodes()     const   { return fRFNendnodes; }

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      void MakeClassRuleCuts( std::ostream& ) const;

      void MakeClassLinear( std::ostream& ) const;

      // get help message text
      void GetHelpMessage() const;

      // initialize rulefit
      void Init( void );

      // copy all training events into a stl::vector
      void InitEventSample( void );

      // initialize monitor ntuple
      void InitMonitorNtuple();

      void TrainTMVARuleFit();
      void TrainJFRuleFit();

   private:

      // check variable range and set var to lower or upper if out of range
      template<typename T>
         inline Bool_t VerifyRange( MsgLogger& mlog, const char *varstr, T& var, const T& vmin, const T& vmax );

      template<typename T>
         inline Bool_t VerifyRange( MsgLogger& mlog, const char *varstr, T& var, const T& vmin, const T& vmax, const T& vdef );

      template<typename T>
         inline Int_t VerifyRange( const T& var, const T& vmin, const T& vmax );

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      RuleFit                      fRuleFit;        // RuleFit instance
      std::vector<TMVA::Event *>   fEventSample;    // the complete training sample
      Double_t                     fSignalFraction; // scalefactor for bkg events to modify initial s/b fraction in training data

      // ntuple
      TTree                       *fMonitorNtuple;  // pointer to monitor rule ntuple
      Double_t                     fNTImportance;   // ntuple: rule importance
      Double_t                     fNTCoefficient;  // ntuple: rule coefficient
      Double_t                     fNTSupport;      // ntuple: rule support
      Int_t                        fNTNcuts;        // ntuple: rule number of cuts
      Int_t                        fNTNvars;        // ntuple: rule number of vars
      Double_t                     fNTPtag;         // ntuple: rule P(tag)
      Double_t                     fNTPss;          // ntuple: rule P(tag s, true s)
      Double_t                     fNTPsb;          // ntuple: rule P(tag s, true b)
      Double_t                     fNTPbs;          // ntuple: rule P(tag b, true s)
      Double_t                     fNTPbb;          // ntuple: rule P(tag b, true b)
      Double_t                     fNTSSB;          // ntuple: rule S/(S+B)
      Int_t                        fNTType;         // ntuple: rule type (+1->signal, -1->bkg)

      // options
      TString                      fRuleFitModuleS;// which rulefit module to use
      Bool_t                       fUseRuleFitJF;  // if true interface with J.Friedmans RuleFit module
      TString                      fRFWorkDir;     // working directory from Friedmans module
      Int_t                        fRFNrules;      // max number of rules (only Friedmans module)
      Int_t                        fRFNendnodes;   // max number of rules (only Friedmans module)
      std::vector<DecisionTree *>  fForest;        // the forest
      Int_t                        fNTrees;        // number of trees in forest
      Double_t                     fTreeEveFrac;   // fraction of events used for training each tree
      SeparationBase              *fSepType;       // the separation used in node splitting
      Double_t                     fMinFracNEve;   // min fraction of number events
      Double_t                     fMaxFracNEve;   // ditto max
      Int_t                        fNCuts;         // grid used in cut applied in node splitting
      TString                      fSepTypeS;        // forest generation: separation type - see DecisionTree
      TString                      fPruneMethodS;    // forest generation: prune method - see DecisionTree
      TMVA::DecisionTree::EPruneMethod fPruneMethod; // forest generation: method used for pruning - see DecisionTree
      Double_t                     fPruneStrength;   // forest generation: prune strength - see DecisionTree
      TString                      fForestTypeS;     // forest generation: how the trees are generated
      Bool_t                       fUseBoost;        // use boosted events for forest generation
      //
      Double_t                     fGDPathEveFrac; //  GD path: fraction of subsamples used for the fitting
      Double_t                     fGDValidEveFrac; // GD path: fraction of subsamples used for the fitting
      Double_t                     fGDTau;          // GD path: def threshold fraction [0..1]
      Double_t                     fGDTauPrec;      // GD path: precision of estimated tau
      Double_t                     fGDTauMin;       // GD path: min threshold fraction [0..1]
      Double_t                     fGDTauMax;       // GD path: max threshold fraction [0..1]
      UInt_t                       fGDTauScan;      // GD path: number of points to scan
      Double_t                     fGDPathStep;     // GD path: step size in path
      Int_t                        fGDNPathSteps;   // GD path: number of steps
      Double_t                     fGDErrScale;     // GD path: stop
      Double_t                     fMinimp;         // rule/linear: minimum importance
      //
      TString                      fModelTypeS;     // rule ensemble: which model (rule,linear or both)
      Double_t                     fRuleMinDist;    // rule min distance - see RuleEnsemble
      Double_t                     fLinQuantile;    // quantile cut to remove outliers - see RuleEnsemble

      ClassDef(MethodRuleFit,0);  // Friedman's RuleFit method
   };

} // namespace TMVA


//_______________________________________________________________________
template<typename T>
inline Int_t TMVA::MethodRuleFit::VerifyRange( const T& var, const T& vmin, const T& vmax )
{
   // check range and return +1 if above, -1 if below or 0 if inside
   if (var>vmax) return  1;
   if (var<vmin) return -1;
   return 0;
}

//_______________________________________________________________________
template<typename T>
inline Bool_t TMVA::MethodRuleFit::VerifyRange( TMVA::MsgLogger& mlog, const char *varstr, T& var, const T& vmin, const T& vmax )
{
   // verify range and print out message
   // if outside range, set to closest limit
   Int_t dir = TMVA::MethodRuleFit::VerifyRange(var,vmin,vmax);
   Bool_t modif=kFALSE;
   if (dir==1) {
      modif = kTRUE;
      var=vmax;
   }
   if (dir==-1) {
      modif = kTRUE;
      var=vmin;
   }
   if (modif) {
      mlog << kWARNING << "Option <" << varstr << "> " << (dir==1 ? "above":"below") << " allowed range. Reset to new value = " << var << Endl;
   }
   return modif;
}

//_______________________________________________________________________
template<typename T>
inline Bool_t TMVA::MethodRuleFit::VerifyRange( TMVA::MsgLogger& mlog, const char *varstr, T& var, const T& vmin, const T& vmax, const T& vdef )
{
   // verify range and print out message
   // if outside range, set to given default value
   Int_t dir = TMVA::MethodRuleFit::VerifyRange(var,vmin,vmax);
   Bool_t modif=kFALSE;
   if (dir!=0) {
      modif = kTRUE;
      var=vdef;
   }
   if (modif) {
      mlog << kWARNING << "Option <" << varstr << "> " << (dir==1 ? "above":"below") << " allowed range. Reset to default value = " << var << Endl;
   }
   return modif;
}


#endif // MethodRuleFit_H
