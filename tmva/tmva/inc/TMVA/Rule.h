// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Rule                                                                  *
 *                                                                                *
 * Description:                                                                   *
 *      A class describing a 'rule'                                               *
 *      Each internal node of a tree defines a rule from all the parental nodes.  *
 *      A rule consists of at least 2 nodes.                                      *
 *      Input: a decision tree (in the constructor)                               *
 *             its coefficient                                                    *
 *                                                                                *
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

#ifndef ROOT_TMVA_Rule
#define ROOT_TMVA_Rule

#include "TMath.h"
#include <vector>
#include <iostream>

#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/RuleCut.h"

namespace TMVA {

   class RuleEnsemble;
   class MsgLogger;
   class Rule;

   std::ostream& operator<<( std::ostream& os, const Rule & rule );

   class Rule {

      // output operator for a Rule
      friend std::ostream& operator<< ( std::ostream& os, const Rule & rule );

   public:

      // main constructor
      Rule( RuleEnsemble *re, const std::vector< const TMVA::Node * > & nodes );

      // main constructor
      Rule( RuleEnsemble *re );

      // copy constructor
      Rule( const Rule & other ) { Copy( other ); }

      // empty constructor
      Rule();

      virtual ~Rule();

      // set message type
      void SetMsgType( EMsgType t );

      // set RuleEnsemble ptr
      void SetRuleEnsemble( const RuleEnsemble *re ) { fRuleEnsemble = re; }

      // set RuleCut ptr
      void SetRuleCut( RuleCut *rc )           { fCut = rc; }

      // set Rule norm
      void SetNorm(Double_t norm)       { fNorm = (norm>0 ? 1.0/norm:1.0); }

      // set coefficient
      void SetCoefficient(Double_t v)   { fCoefficient=v; }

      // set support
      void SetSupport(Double_t v)       { fSupport=v; fSigma = TMath::Sqrt(v*(1.0-v));}

      // set s/(s+b)
      void SetSSB(Double_t v)           { fSSB=v; }

      // set N(eve) accepted by rule
      void SetSSBNeve(Double_t v)       { fSSBNeve=v; }

      // set reference importance
      void SetImportanceRef(Double_t v) { fImportanceRef=(v>0 ? v:1.0); }

      // calculate importance
      void CalcImportance()             { fImportance = TMath::Abs(fCoefficient)*fSigma; }

      // get the relative importance
      Double_t GetRelImportance()  const { return fImportance/fImportanceRef; }

      // evaluate the Rule for the given Event using the coefficient
      //      inline Double_t EvalEvent( const Event& e, Bool_t norm ) const;

      // evaluate the Rule for the given Event, not using normalization or the coefficient
      inline Bool_t EvalEvent( const Event& e ) const;

      // test if two rules are equal
      Bool_t Equal( const Rule & other, Bool_t useCutValue, Double_t maxdist ) const;

      // get distance between two equal (ie apart from the cut values) rules
      Double_t RuleDist( const Rule & other, Bool_t useCutValue ) const;

      // returns true if the trained S/(S+B) of the last node is > 0.5
      Double_t GetSSB()       const { return fSSB; }
      Double_t GetSSBNeve()   const { return fSSBNeve; }
      Bool_t   IsSignalRule() const { return (fSSB>0.5); }

      // copy operator
      void operator=( const Rule & other )  { Copy( other ); }

      // identical operator
      Bool_t operator==( const Rule & other ) const;

      Bool_t operator<( const Rule & other ) const;

      // get number of variables used in Rule
      UInt_t GetNumVarsUsed() const { return fCut->GetNvars(); }

      // get number of cuts in Rule
      UInt_t GetNcuts() const { return fCut->GetNcuts(); }

      // check if variable is used by the rule
      Bool_t ContainsVariable(UInt_t iv) const;

      // accessors
      const RuleCut*      GetRuleCut()       const { return fCut; }
      const RuleEnsemble* GetRuleEnsemble()  const { return fRuleEnsemble; }
      Double_t            GetCoefficient()   const { return fCoefficient; }
      Double_t            GetSupport()       const { return fSupport; }
      Double_t            GetSigma()         const { return fSigma; }
      Double_t            GetNorm()          const { return fNorm; }
      Double_t            GetImportance()    const { return fImportance; }
      Double_t            GetImportanceRef() const { return fImportanceRef; }

      // print the rule using flogger
      void PrintLogger( const char *title=0 ) const;

      // print just the raw info, used for weight file generation
      void  PrintRaw   ( std::ostream& os  ) const; // obsolete
      void* AddXMLTo   ( void* parent ) const;

      void  ReadRaw    ( std::istream& os    ); // obsolete
      void  ReadFromXML( void* wghtnode );

   private:

      // set sigma - don't use this as non private!
      void SetSigma(Double_t v)         { fSigma=v; }

      // print info about the Rule
      void Print( std::ostream& os ) const;

      // copy from another rule
      void Copy( const Rule & other );

      // get the name of variable with index i
      const TString & GetVarName( Int_t i) const;

      RuleCut*             fCut;           // all cuts associated with the rule
      Double_t             fNorm;          // normalization - usually 1.0/t(k)
      Double_t             fSupport;       // s(k)
      Double_t             fSigma;         // t(k) = sqrt(s*(1-s))
      Double_t             fCoefficient;   // rule coeff. a(k)
      Double_t             fImportance;    // importance of rule
      Double_t             fImportanceRef; // importance ref
      const RuleEnsemble*  fRuleEnsemble;  // pointer to parent RuleEnsemble
      Double_t             fSSB;           // S/(S+B) for rule
      Double_t             fSSBNeve;       // N(events) reaching the last node in reevaluation

      mutable MsgLogger*   fLogger;        //! message logger
      MsgLogger& Log() const { return *fLogger; }

   };

} // end of TMVA namespace

//_______________________________________________________________________
inline Bool_t TMVA::Rule::EvalEvent( const TMVA::Event& e ) const
{
   // Checks if event is accepted by rule.
   // Return true if yes and false if not.
   //
   return fCut->EvalEvent(e);
}

#endif
