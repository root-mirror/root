// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RuleEnsemble                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      A class generating an ensemble of rules                                   *
 *      Input:  a forest of decision trees                                        *
 *      Output: an ensemble of rules                                              *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-KP Heidelberg, Ger. *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U.                                                             *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_RuleEnsemble
#define ROOT_TMVA_RuleEnsemble

#include "TMath.h"

#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/Rule.h"
#include "TMVA/Types.h"

#include <vector>

class TH1F;

namespace TMVA {

   class MethodBase;
   class RuleFit;
   class MethodRuleFit;
   class RuleEnsemble;
   class MsgLogger;

   std::ostream& operator<<( std::ostream& os, const RuleEnsemble& event );

   class RuleEnsemble {

      // output operator for a RuleEnsemble
      friend std::ostream& operator<< ( std::ostream& os, const RuleEnsemble& rules );

   public:

      enum ELearningModel { kFull=0, kRules=1, kLinear=2 };

      // main constructor
      RuleEnsemble( RuleFit* rf );

      // copy constructor
      RuleEnsemble( const RuleEnsemble& other );

      // empty constructor
      RuleEnsemble();

      // destructor
      virtual ~RuleEnsemble();

      // initialize
      void Initialize( const RuleFit* rf );

      // set message type
      void SetMsgType( EMsgType t );

      // makes the model - calls MakeRules() and MakeLinearTerms()
      void MakeModel();

      // generates the rules from a given forest of decision trees
      void MakeRules( const std::vector< const TMVA::DecisionTree *>& forest );

      // make the linear terms
      void MakeLinearTerms();

      // select linear model
      void SetModelLinear() { fLearningModel = kLinear; }

      // select rule model
      void SetModelRules()  { fLearningModel = kRules; }

      // select full (linear+rules) model
      void SetModelFull()   { fLearningModel = kFull; }

      // set rule collection (if not created by MakeRules())
      void SetRules( const std::vector< TMVA::Rule *> & rules );

      // set RuleFit ptr
      void SetRuleFit( const RuleFit *rf ) { fRuleFit = rf; }

      // set coefficients
      void  SetCoefficients( const std::vector< Double_t >& v );
      void  SetCoefficient( UInt_t i, Double_t v )                  { if (i<fRules.size()) fRules[i]->SetCoefficient(v); }
      //
      void  SetOffset(Double_t v=0.0)                               { fOffset=v; }
      void  AddOffset(Double_t v)                                   { fOffset+=v; }
      void  SetLinCoefficients( const std::vector< Double_t >& v )  { fLinCoefficients = v; }
      void  SetLinCoefficient( UInt_t i, Double_t v )               { fLinCoefficients[i] = v; }
      void  SetLinDM( const std::vector<Double_t>   & xmin ) { fLinDM   = xmin; }
      void  SetLinDP( const std::vector<Double_t>   & xmax ) { fLinDP   = xmax; }
      void  SetLinNorm( const std::vector<Double_t> & norm ) { fLinNorm = norm; }

      Double_t CalcLinNorm( Double_t stdev ) { return ( stdev>0 ? fAverageRuleSigma/stdev : 1.0 ); }

      // clear coefficients
      void  ClearCoefficients( Double_t val=0 )    { for (UInt_t i=0; i<fRules.size(); i++)           fRules[i]->SetCoefficient(val); }
      void  ClearLinCoefficients( Double_t val=0 ) { for (UInt_t i=0; i<fLinCoefficients.size(); i++) fLinCoefficients[i]=val; }
      void  ClearLinNorm( Double_t val=1.0 )       { for (UInt_t i=0; i<fLinNorm.size(); i++)         fLinNorm[i]=val; }

      // set maximum allowed distance between equal rules
      void SetRuleMinDist(Double_t d)          { fRuleMinDist = d; }

      // set minimum rule importance - used by CleanupRules()
      void SetImportanceCut(Double_t minimp=0) { fImportanceCut=minimp; }

      // set the quantile for linear terms
      void SetLinQuantile(Double_t q)          { fLinQuantile=q; }

      // set average sigma for rules
      void SetAverageRuleSigma(Double_t v) { if (v>0.5) v=0.5; fAverageRuleSigma = v; fAverageSupport = 0.5*(1.0+TMath::Sqrt(1.0-4.0*v*v)); }

      // Calculate the number of possible rules from a given tree
      Int_t CalcNRules( const TMVA::DecisionTree* dtree );
      // Recursively search for end-nodes; used by CalcNRules()
      void  FindNEndNodes( const TMVA::Node* node, Int_t& nendnodes );

      // set current event to be used
      void SetEvent( const Event & e ) { fEvent = &e; fEventCacheOK = kFALSE; }

      // fill cached values of rule/linear respons
      void UpdateEventVal();

      // fill binary rule respons for all events (or selected subset)
      void MakeRuleMap(const std::vector<const TMVA::Event *> *events=0, UInt_t ifirst=0, UInt_t ilast=0);

      // clear rule map
      void ClearRuleMap() { fRuleMap.clear(); fRuleMapEvents=0; }

      // evaluates the event using the ensemble of rules
      // the following uses fEventCache, that is per event saved in cache
      Double_t EvalEvent() const;
      Double_t EvalEvent( const Event & e );

      // same as previous but using other model coefficients
      Double_t EvalEvent( Double_t ofs,
                          const std::vector<Double_t> & coefs,
                          const std::vector<Double_t> & lincoefs) const;
      Double_t EvalEvent( const Event & e,
                          Double_t ofs,
                          const std::vector<Double_t> & coefs,
                          const std::vector<Double_t> & lincoefs);

      // same as above but using the event index
      // these will use fRuleMap - MUST call MakeRuleMap() before - no check...
      Double_t EvalEvent( UInt_t evtidx ) const;
      Double_t EvalEvent( UInt_t evtidx,
                          Double_t ofs,
                          const std::vector<Double_t> & coefs,
                          const std::vector<Double_t> & lincoefs) const;

      // evaluate the linear term using event by reference
      //      Double_t EvalLinEvent( UInt_t vind ) const;
      Double_t EvalLinEvent() const;
      Double_t EvalLinEvent( const std::vector<Double_t> & coefs ) const;
      Double_t EvalLinEvent( const Event &e );
      Double_t EvalLinEvent( const Event &e, UInt_t vind );
      Double_t EvalLinEvent( const Event &e, const std::vector<Double_t> & coefs );

      // idem but using evtidx - must call MakeRuleMap() first
      Double_t EvalLinEvent( UInt_t evtidx ) const;
      Double_t EvalLinEvent( UInt_t evtidx, const std::vector<Double_t> & coefs ) const;
      Double_t EvalLinEvent( UInt_t evtidx, UInt_t vind ) const;
      Double_t EvalLinEvent( UInt_t evtidx, UInt_t vind, Double_t coefs ) const;

      // evaluate linear terms used to fill fEventLinearVal
      Double_t EvalLinEventRaw( UInt_t vind, const Event &e, Bool_t norm ) const;
      Double_t EvalLinEventRaw( UInt_t vind, UInt_t evtidx,  Bool_t norm ) const;

      // calculate p(y=1|x) for a given event using the linear terms
      Double_t PdfLinear( Double_t & nsig, Double_t & ntot ) const;

      // calculate p(y=1|x) for a given event using the rules
      Double_t PdfRule( Double_t & nsig, Double_t & ntot ) const;

      // calculate F* = 2*p(y=1|x) - 1
      Double_t FStar() const;
      Double_t FStar(const TMVA::Event & e );

      // set reference importance for all model objects
      void SetImportanceRef(Double_t impref);

      // calculates the support for all rules given the set of events
      void CalcRuleSupport();

      // calculates rule importance
      void CalcImportance();

      // calculates rule importance
      Double_t CalcRuleImportance();

      // calculates linear importance
      Double_t CalcLinImportance();

      // calculates variable importance
      void CalcVarImportance();

      // remove rules of low importance
      void CleanupRules();

      // remove linear terms of low importance
      void CleanupLinear();

      // remove similar rules
      void RemoveSimilarRules();

      // get rule statistics
      void RuleStatistics();

      // get rule response stats
      void RuleResponseStats();

      // copy operator
      void operator=( const RuleEnsemble& other ) { Copy( other ); }

      // calculate sum of the squared coefficients
      Double_t CoefficientRadius();

      // fill the vector with the coefficients
      void GetCoefficients( std::vector< Double_t >& v );

      // accessors
      const MethodRuleFit*                   GetMethodRuleFit()   const;
      const MethodBase*                      GetMethodBase()      const;
      const RuleFit*                         GetRuleFit()         const { return fRuleFit; }
      //
      const std::vector<const TMVA::Event *>*     GetTrainingEvents()  const;
      const Event*                    GetTrainingEvent(UInt_t i) const;
      const Event*                    GetEvent() const { return fEvent; }
      //
      Bool_t                          DoLinear()             const { return (fLearningModel==kFull) || (fLearningModel==kLinear); }
      Bool_t                          DoRules()              const { return (fLearningModel==kFull) || (fLearningModel==kRules); }
      Bool_t                          DoOnlyRules()          const { return (fLearningModel==kRules); }
      Bool_t                          DoOnlyLinear()         const { return (fLearningModel==kLinear); }
      Bool_t                          DoFull()               const { return (fLearningModel==kFull); }
      ELearningModel                  GetLearningModel()     const { return fLearningModel; }
      Double_t                        GetImportanceCut()     const { return fImportanceCut; }
      Double_t                        GetImportanceRef()     const { return fImportanceRef; }
      Double_t                        GetOffset()            const { return fOffset; }
      UInt_t                          GetNRules()            const { return (DoRules() ? fRules.size():0); }
      const std::vector<TMVA::Rule*>& GetRulesConst()        const { return fRules; }
      std::vector<TMVA::Rule*>&       GetRules()                   { return fRules; }
      const std::vector< Double_t >&  GetLinCoefficients()   const { return fLinCoefficients; }
      const std::vector< Double_t >&  GetLinNorm()           const { return fLinNorm; }
      const std::vector< Double_t >&  GetLinImportance()     const { return fLinImportance; }
      const std::vector< Double_t >&  GetVarImportance()     const { return fVarImportance; }
      UInt_t                          GetNLinear()           const { return (DoLinear() ? fLinNorm.size():0); }
      Double_t                        GetLinQuantile()       const { return fLinQuantile; }

      const Rule    *GetRulesConst(int i)        const { return fRules[i]; }
      Rule          *GetRules(int i)                   { return fRules[i]; }

      UInt_t         GetRulesNCuts(int i)        const { return fRules[i]->GetRuleCut()->GetNcuts(); }
      Double_t       GetRuleMinDist()            const { return fRuleMinDist; }
      Double_t       GetLinCoefficients(int i)   const { return fLinCoefficients[i]; }
      Double_t       GetLinNorm(int i)           const { return fLinNorm[i]; }
      Double_t       GetLinDM(int i)             const { return fLinDM[i]; }
      Double_t       GetLinDP(int i)             const { return fLinDP[i]; }
      Double_t       GetLinImportance(int i)     const { return fLinImportance[i]; }
      Double_t       GetVarImportance(int i)     const { return fVarImportance[i]; }
      Double_t       GetRulePTag(int i)          const { return fRulePTag[i]; }
      Double_t       GetRulePSS(int i)           const { return fRulePSS[i]; }
      Double_t       GetRulePSB(int i)           const { return fRulePSB[i]; }
      Double_t       GetRulePBS(int i)           const { return fRulePBS[i]; }
      Double_t       GetRulePBB(int i)           const { return fRulePBB[i]; }

      Bool_t         IsLinTermOK(int i)          const { return fLinTermOK[i]; }
      //
      Double_t       GetAverageSupport()             const { return fAverageSupport; }
      Double_t       GetAverageRuleSigma()           const { return fAverageRuleSigma; }
      Double_t       GetEventRuleVal(UInt_t i)       const { return (fEventRuleVal[i] ? 1.0:0.0); }
      Double_t       GetEventLinearVal(UInt_t i)     const { return fEventLinearVal[i]; }
      Double_t       GetEventLinearValNorm(UInt_t i) const { return fEventLinearVal[i]*fLinNorm[i]; }
      //
      const std::vector<UInt_t>  & GetEventRuleMap(UInt_t evtidx) const { return fRuleMap[evtidx]; }
      const TMVA::Event *GetRuleMapEvent(UInt_t evtidx) const { return (*fRuleMapEvents)[evtidx]; }
      Bool_t         IsRuleMapOK()               const { return fRuleMapOK; }

      // print rule generation info
      void  PrintRuleGen() const;

      // print the ensemble
      void  Print() const;

      // print the model in a cryptic way
      void  PrintRaw   ( std::ostream& os  ) const; // obsolete
      void* AddXMLTo   ( void* parent ) const;

      // read the model from input stream
      void  ReadRaw    ( std::istream& istr ); // obsolete
      void  ReadFromXML( void* wghtnode );


   private:

      // delete all rules
      void DeleteRules() { for (UInt_t i=0; i<fRules.size(); i++) delete fRules[i]; fRules.clear(); }

      // copy method
      void  Copy( RuleEnsemble const& other );

      // set all coeffs to default values
      void  ResetCoefficients();

      // make rules form one decision tree
      void  MakeRulesFromTree( const DecisionTree *dtree );

      // add a rule with the given end-node
      void  AddRule( const Node *node );

      // make a rule
      Rule *MakeTheRule( const Node *node );


      ELearningModel                fLearningModel;     // can be full (rules+linear), rules, linear
      Double_t                      fImportanceCut;     // minimum importance accepted
      Double_t                      fLinQuantile;       // quantile cut to remove outliers
      Double_t                      fOffset;            // offset in discriminator function
      std::vector< TMVA::Rule* >    fRules;             // vector of rules
      std::vector< Char_t >         fLinTermOK;         // flags linear terms with sufficient strong importance <-- stores boolean
      std::vector< Double_t >       fLinDP;             // delta+ in eq 24, ref 2
      std::vector< Double_t >       fLinDM;             // delta-
      std::vector< Double_t >       fLinCoefficients;   // linear coefficients, one per variable
      std::vector< Double_t >       fLinNorm;           // norm of ditto, see after eq 26 in ref 2
      std::vector< TH1F* >          fLinPDFB;           // pdfs for each variable, background
      std::vector< TH1F* >          fLinPDFS;           // pdfs for each variable, signal
      std::vector< Double_t >       fLinImportance;     // linear term importance
      std::vector< Double_t >       fVarImportance;     // one importance per input variable
      Double_t                      fImportanceRef;     // reference importance (max)
      Double_t                      fAverageSupport;    // average support (over all rules)
      Double_t                      fAverageRuleSigma;  // average rule sigma
      //
      std::vector< Double_t >       fRuleVarFrac;       // fraction of rules using a given variable - size of vector = n(variables)
      std::vector< Double_t >       fRulePSS;           // p(tag as S|S) - tagged as S if rule is SIG and the event is accepted
      std::vector< Double_t >       fRulePSB;           // p(tag as S|B)
      std::vector< Double_t >       fRulePBS;           // p(tag as B|S)
      std::vector< Double_t >       fRulePBB;           // p(tag as B|B)
      std::vector< Double_t >       fRulePTag;          // p(tag)
      Double_t                      fRuleFSig;          // N(sig)/N(sig)+N(bkg)
      Double_t                      fRuleNCave;         // N(cuts) average
      Double_t                      fRuleNCsig;         // idem sigma
      //
      Double_t                      fRuleMinDist;       // minimum rule distance
      UInt_t                        fNRulesGenerated;   // number of rules generated, before cleanup
      //
      const Event*                  fEvent;             // current event.
      Bool_t                        fEventCacheOK;      // true if rule/linear respons are updated
      std::vector<Char_t>           fEventRuleVal;      // the rule respons of current event <----- stores boolean
      std::vector<Double_t>         fEventLinearVal;    // linear respons
      //
      Bool_t                        fRuleMapOK;         // true if MakeRuleMap() has been called
      std::vector< std::vector<UInt_t> > fRuleMap;           // map of rule responses
      UInt_t                        fRuleMapInd0;       // start index
      UInt_t                        fRuleMapInd1;       // last index
      const std::vector<const TMVA::Event *> *fRuleMapEvents; // pointer to vector of events used
      //
      const RuleFit*                fRuleFit;           // pointer to rule fit object

      mutable MsgLogger*            fLogger;            //! message logger
      MsgLogger& Log() const { return *fLogger; }
   };
}

//_______________________________________________________________________
inline void TMVA::RuleEnsemble::UpdateEventVal()
{
   //
   // Update rule and linear respons using the current event
   //
   if (fEventCacheOK) return;
   //
   if (DoRules()) {
      UInt_t nrules = fRules.size();
      fEventRuleVal.resize(nrules,kFALSE);
      for (UInt_t r=0; r<nrules; r++) {
         fEventRuleVal[r] = fRules[r]->EvalEvent(*fEvent);
      }
   }
   if (DoLinear()) {
      UInt_t nlin = fLinTermOK.size();
      fEventLinearVal.resize(nlin,0);
      for (UInt_t r=0; r<nlin; r++) {
         fEventLinearVal[r] = EvalLinEventRaw(r,*fEvent,kFALSE); // not normalised!
      }
   }
   fEventCacheOK = kTRUE;
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent() const
{
   // evaluate current event

   Int_t nrules = fRules.size();
   Double_t rval=fOffset;
   Double_t linear=0;
   //
   // evaluate all rules
   // normally it should NOT use the normalized rules - the flag should be kFALSE
   //
   if (DoRules()) {
      for ( Int_t i=0; i<nrules; i++ ) {
         if (fEventRuleVal[i])
            rval += fRules[i]->GetCoefficient();
      }
   }
   //
   // Include linear part - the call below incorporates both coefficient and normalisation (fLinNorm)
   //
   if (DoLinear()) linear = EvalLinEvent();
   rval +=linear;

   return rval;
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent( Double_t ofs,
                                               const std::vector<Double_t> & coefs,
                                               const std::vector<Double_t> & lincoefs ) const
{
   // evaluate current event with given offset and coefs

   Int_t nrules    = fRules.size();
   Double_t rval   = ofs;
   Double_t linear = 0;
   //
   // evaluate all rules
   //
   if (DoRules()) {
      for ( Int_t i=0; i<nrules; i++ ) {
         if (fEventRuleVal[i])
            rval += coefs[i];
      }
   }
   //
   // Include linear part - the call below incorporates both coefficient and normalisation (fLinNorm)
   //
   if (DoLinear()) linear = EvalLinEvent(lincoefs);
   rval +=linear;

   return rval;
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent(const TMVA::Event & e)
{
   // evaluate event e
   SetEvent(e);
   UpdateEventVal();
   return EvalEvent();
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent(const TMVA::Event & e,
                                              Double_t ofs,
                                              const std::vector<Double_t> & coefs,
                                              const std::vector<Double_t> & lincoefs )
{
   // evaluate event e
   SetEvent(e);
   UpdateEventVal();
   return EvalEvent(ofs,coefs,lincoefs);
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent(UInt_t evtidx) const
{
   // evaluate event with index evtidx
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   //
   Double_t rval=fOffset;
   if (DoRules()) {
      UInt_t nrules = fRuleMap[evtidx].size();
      UInt_t rind;
      for (UInt_t ir = 0; ir<nrules; ir++) {
         rind = fRuleMap[evtidx][ir];
         rval += fRules[rind]->GetCoefficient();
      }
   }
   if (DoLinear()) {
      UInt_t nlin = fLinTermOK.size();
      for (UInt_t r=0; r<nlin; r++) {
         if (fLinTermOK[r]) {
            rval += fLinCoefficients[r] * EvalLinEventRaw(r,*(*fRuleMapEvents)[evtidx],kTRUE);
         }
      }
   }
   return rval;
}

//_____________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalEvent(UInt_t evtidx,
                                              Double_t ofs,
                                              const std::vector<Double_t> & coefs,
                                              const std::vector<Double_t> & lincoefs ) const
{
   // evaluate event with index evtidx and user given model coefficients
   //
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   Double_t rval=ofs;
   if (DoRules()) {
      UInt_t nrules = fRuleMap[evtidx].size();
      UInt_t rind;
      for (UInt_t ir = 0; ir<nrules; ir++) {
         rind = fRuleMap[evtidx][ir];
         rval += coefs[rind];
      }
   }
   if (DoLinear()) {
      rval += EvalLinEvent( evtidx, lincoefs );
   }
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEventRaw( UInt_t vind, const TMVA::Event & e, Bool_t norm) const
{
   // evaluate the event linearly (not normalized)

   Double_t val  = e.GetValue(vind);
   Double_t rval = TMath::Min( fLinDP[vind], TMath::Max( fLinDM[vind], val ) );
   if (norm) rval *= fLinNorm[vind];
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEventRaw( UInt_t vind, UInt_t evtidx, Bool_t norm) const
{
   // evaluate the event linearly (not normalized)

   Double_t val  = (*fRuleMapEvents)[evtidx]->GetValue(vind);
   Double_t rval = TMath::Min( fLinDP[vind], TMath::Max( fLinDM[vind], val ) );
   if (norm) rval *= fLinNorm[vind];
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent() const
{
   // evaluate event linearly

   Double_t rval=0;
   for (UInt_t v=0; v<fLinTermOK.size(); v++) {
      if (fLinTermOK[v])
         rval += fLinCoefficients[v]*fEventLinearVal[v]*fLinNorm[v];
   }
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent(const std::vector<Double_t> & coefs) const
{
   // evaluate event linearly using the given coefficients

   Double_t rval=0;
   for (UInt_t v=0; v<fLinTermOK.size(); v++) {
      if (fLinTermOK[v])
         rval += coefs[v]*fEventLinearVal[v]*fLinNorm[v];
   }
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( const TMVA::Event& e )
{
   // evaluate event linearly

   SetEvent(e);
   UpdateEventVal();
   return EvalLinEvent();
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( const TMVA::Event& e, UInt_t vind )
{
   // evaluate linear term vind

   SetEvent(e);
   UpdateEventVal();
   return GetEventLinearValNorm(vind);
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( const TMVA::Event& e, const std::vector<Double_t> & coefs )
{
   // evaluate event linearly using the given coefficients

   SetEvent(e);
   UpdateEventVal();
   return EvalLinEvent(coefs);
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( UInt_t evtidx, const std::vector<Double_t> & coefs ) const
{
   // evaluate event linearly using the given coefficients
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   Double_t rval=0;
   UInt_t nlin = fLinTermOK.size();
   for (UInt_t r=0; r<nlin; r++) {
      if (fLinTermOK[r]) {
         rval += coefs[r] * EvalLinEventRaw(r,*(*fRuleMapEvents)[evtidx],kTRUE);
      }
   }
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( UInt_t evtidx ) const
{
   // evaluate event linearly using the given coefficients
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   Double_t rval=0;
   UInt_t nlin = fLinTermOK.size();
   for (UInt_t r=0; r<nlin; r++) {
      if (fLinTermOK[r]) {
         rval += fLinCoefficients[r] * EvalLinEventRaw(r,*(*fRuleMapEvents)[evtidx],kTRUE);
      }
   }
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( UInt_t evtidx, UInt_t vind ) const
{
   // evaluate event linearly using the given coefficients
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   Double_t rval;
   rval = fLinCoefficients[vind] * EvalLinEventRaw(vind,*(*fRuleMapEvents)[evtidx],kTRUE);
   return rval;
}

//_______________________________________________________________________
inline Double_t TMVA::RuleEnsemble::EvalLinEvent( UInt_t evtidx, UInt_t vind, Double_t coefs ) const
{
   // evaluate event linearly using the given coefficients
   if ((evtidx<fRuleMapInd0) || (evtidx>fRuleMapInd1)) return 0;
   Double_t rval;
   rval = coefs * EvalLinEventRaw(vind,*(*fRuleMapEvents)[evtidx],kTRUE);
   return rval;
}

#endif
