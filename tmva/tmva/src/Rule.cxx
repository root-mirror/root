// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Rule                                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      A class describing a 'rule'                                               *
 *      Each internal node of a tree defines a rule from all the parental nodes.  *
 *      A rule consists of at least 2 nodes.                                      *
 *      Input: a decision tree (in the constructor)                               *
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

/*! \class TMVA::Rule
\ingroup TMVA

Implementation of a rule.

A rule is simply a branch or a part of a branch in a tree.
It fulfills the following:

  - First node is the root node of the originating tree
  - Consists of a minimum of 2 nodes
  - A rule returns for a given event:
    - 0 : if the event fails at any node
    - 1 : otherwise
  - If the rule contains <2 nodes, it returns 0 SHOULD NOT HAPPEN!

The coefficient is found by either brute force or some sort of
intelligent fitting. See the RuleEnsemble class for more info.
*/

#include "TMVA/Rule.h"

#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodRuleFit.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/RuleCut.h"
#include "TMVA/RuleFit.h"
#include "TMVA/RuleEnsemble.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include <iomanip>

////////////////////////////////////////////////////////////////////////////////
/// the main constructor for a Rule

TMVA::Rule::Rule( RuleEnsemble *re,
                  const std::vector< const Node * >& nodes )
   : fCut           ( 0 )
   , fNorm          ( 1.0 )
   , fSupport       ( 0.0 )
   , fSigma         ( 0.0 )
   , fCoefficient   ( 0.0 )
   , fImportance    ( 0.0 )
   , fImportanceRef ( 1.0 )
   , fRuleEnsemble  ( re )
   , fSSB           ( 0 )
   , fSSBNeve       ( 0 )
   , fLogger( new MsgLogger("RuleFit") )
{
   //
   // input:
   //   nodes  - a vector of Node; from these all possible rules will be created
   //
   //

   fCut     = new RuleCut( nodes );
   fSSB     = fCut->GetPurity();
   fSSBNeve = fCut->GetCutNeve();
}

////////////////////////////////////////////////////////////////////////////////
/// the simple constructor

TMVA::Rule::Rule( RuleEnsemble *re )
   : fCut           ( 0 )
   , fNorm          ( 1.0 )
   , fSupport       ( 0.0 )
   , fSigma         ( 0.0 )
   , fCoefficient   ( 0.0 )
   , fImportance    ( 0.0 )
   , fImportanceRef ( 1.0 )
   , fRuleEnsemble  ( re )
   , fSSB           ( 0 )
   , fSSBNeve       ( 0 )
   , fLogger( new MsgLogger("RuleFit") )
{
}

////////////////////////////////////////////////////////////////////////////////
/// the simple constructor

TMVA::Rule::Rule()
   : fCut           ( 0 )
   , fNorm          ( 1.0 )
   , fSupport       ( 0.0 )
   , fSigma         ( 0.0 )
   , fCoefficient   ( 0.0 )
   , fImportance    ( 0.0 )
   , fImportanceRef ( 1.0 )
   , fRuleEnsemble  ( 0 )
   , fSSB           ( 0 )
   , fSSBNeve       ( 0 )
   , fLogger( new MsgLogger("RuleFit") )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Rule::~Rule()
{
   delete fCut;
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// check if variable in node

Bool_t TMVA::Rule::ContainsVariable(UInt_t iv) const
{
   Bool_t found    = kFALSE;
   Bool_t doneLoop = kFALSE;
   UInt_t nvars    = fCut->GetNvars();
   UInt_t i        = 0;
   //
   while (!doneLoop) {
      found = (fCut->GetSelector(i) == iv);
      i++;
      doneLoop = (found || (i==nvars));
   }
   return found;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Rule::SetMsgType( EMsgType t )
{
   fLogger->SetMinType(t);
}


////////////////////////////////////////////////////////////////////////////////
/// Compare two rules.
///
///  - useCutValue:
///    - true -> calculate a distance between the two rules based on the cut values
///              if the rule cuts are not equal, the distance is < 0 (-1.0)
///              return true if d<mindist
///    - false-> ignore mindist, return true if rules are equal, ignoring cut values
///  - mindist:     min distance allowed between rules; if < 0 => set useCutValue=false;

Bool_t TMVA::Rule::Equal( const Rule& other, Bool_t useCutValue, Double_t mindist ) const
{
   Bool_t rval=kFALSE;
   if (mindist<0) useCutValue=kFALSE;
   Double_t d = RuleDist( other, useCutValue );
   // cut value used - return true if 0<=d<mindist
   if (useCutValue) rval = ( (!(d<0)) && (d<mindist) );
   else rval = (!(d<0));
   // cut value not used, return true if <> -1
   return rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns:
///
///  * -1.0 : rules are NOT equal, i.e, variables and/or cut directions are wrong
///  * >=0: rules are equal apart from the cutvalue, returns \f$ d = \sqrt{\sum(c1-c2)^2} \f$
///
/// If not useCutValue, the distance is exactly zero if they are equal

Double_t TMVA::Rule::RuleDist( const Rule& other, Bool_t useCutValue ) const
{
   if (fCut->GetNvars()!=other.GetRuleCut()->GetNvars()) return -1.0; // check number of cuts
   //
   const UInt_t nvars  = fCut->GetNvars();
   //
   Int_t    sel;         // cut variable
   Double_t rms;         // rms of cut variable
   Double_t smin;        // distance between the lower range
   Double_t smax;        // distance between the upper range
   Double_t vminA,vmaxA; // min,max range of cut A (cut from this Rule)
   Double_t vminB,vmaxB; // idem from other Rule
   //
   // compare nodes
   // A 'distance' is assigned if the two rules has exactly the same set of cuts but with
   // different cut values.
   // The distance is given in number of sigmas
   //
   UInt_t   in     = 0;    // cut index
   Double_t sumdc2 = 0;    // sum of 'distances'
   Bool_t   equal  = true; // flag if cut are equal
   //
   const RuleCut *otherCut = other.GetRuleCut();
   while ((equal) && (in<nvars)) {
      // check equality in cut topology
      equal = ( (fCut->GetSelector(in) == (otherCut->GetSelector(in))) &&
                (fCut->GetCutDoMin(in) == (otherCut->GetCutDoMin(in))) &&
                (fCut->GetCutDoMax(in) == (otherCut->GetCutDoMax(in))) );
      // if equal topology, check cut values
      if (equal) {
         if (useCutValue) {
            sel   = fCut->GetSelector(in);
            vminA = fCut->GetCutMin(in);
            vmaxA = fCut->GetCutMax(in);
            vminB = other.GetRuleCut()->GetCutMin(in);
            vmaxB = other.GetRuleCut()->GetCutMax(in);
            // messy - but ok...
            rms = fRuleEnsemble->GetRuleFit()->GetMethodBase()->GetRMS(sel);
            smin=0;
            smax=0;
            if (fCut->GetCutDoMin(in))
               smin = ( rms>0 ? (vminA-vminB)/rms : 0 );
            if (fCut->GetCutDoMax(in))
               smax = ( rms>0 ? (vmaxA-vmaxB)/rms : 0 );
            sumdc2 += smin*smin + smax*smax;
            //            sumw   += 1.0/(rms*rms); // TODO: probably not needed
         }
      }
      in++;
   }
   if (!useCutValue) sumdc2 = (equal ? 0.0:-1.0); // ignore cut values
   else              sumdc2 = (equal ? sqrt(sumdc2) : -1.0);

   return sumdc2;
}

////////////////////////////////////////////////////////////////////////////////
/// comparison operator ==

Bool_t TMVA::Rule::operator==( const Rule& other ) const
{
   return this->Equal( other, kTRUE, 1e-3 );
}

////////////////////////////////////////////////////////////////////////////////
/// comparison operator <

Bool_t TMVA::Rule::operator<( const Rule& other ) const
{
   return (fImportance < other.GetImportance());
}

////////////////////////////////////////////////////////////////////////////////
/// std::ostream operator

std::ostream& TMVA::operator<< ( std::ostream& os, const Rule& rule )
{
   rule.Print( os );
   return os;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the name of a rule

const TString & TMVA::Rule::GetVarName( Int_t i ) const
{
   return fRuleEnsemble->GetMethodBase()->GetInputLabel(i);
}

////////////////////////////////////////////////////////////////////////////////
/// copy function

void TMVA::Rule::Copy( const Rule& other )
{
   if(this != &other) {
      SetRuleEnsemble( other.GetRuleEnsemble() );
      fCut = new RuleCut( *(other.GetRuleCut()) );
      fSSB     = other.GetSSB();
      fSSBNeve = other.GetSSBNeve();
      SetCoefficient(other.GetCoefficient());
      SetSupport( other.GetSupport() );
      SetSigma( other.GetSigma() );
      SetNorm( other.GetNorm() );
      CalcImportance();
      SetImportanceRef( other.GetImportanceRef() );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print function

void TMVA::Rule::Print( std::ostream& os ) const
{
   const UInt_t nvars = fCut->GetNvars();
   if (nvars<1) os << "     *** WARNING - <EMPTY RULE> ***" << std::endl; // TODO: Fix this, use fLogger
   //
   Int_t sel;
   Double_t valmin, valmax;
   //
   os << "    Importance  = " << Form("%1.4f", fImportance/fImportanceRef) << std::endl;
   os << "    Coefficient = " << Form("%1.4f", fCoefficient) << std::endl;
   os << "    Support     = " << Form("%1.4f", fSupport)  << std::endl;
   os << "    S/(S+B)     = " << Form("%1.4f", fSSB)  << std::endl;

   for ( UInt_t i=0; i<nvars; i++) {
      os << "    ";
      sel    = fCut->GetSelector(i);
      valmin = fCut->GetCutMin(i);
      valmax = fCut->GetCutMax(i);
      //
      os << Form("* Cut %2d",i+1) << " : " << std::flush;
      if (fCut->GetCutDoMin(i)) os << Form("%10.3g",valmin) << " < " << std::flush;
      else                      os << "             " << std::flush;
      os << GetVarName(sel) << std::flush;
      if (fCut->GetCutDoMax(i)) os << " < " << Form("%10.3g",valmax) << std::flush;
      else                      os << "             " << std::flush;
      os << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print function

void TMVA::Rule::PrintLogger(const char *title) const
{
   const UInt_t nvars = fCut->GetNvars();
   if (nvars<1) Log() << kWARNING << "BUG TRAP: EMPTY RULE!!!" << Endl;
   //
   Int_t sel;
   Double_t valmin, valmax;
   //
   if (title) Log() << kINFO << title;
   Log() << kINFO
         << "Importance  = " << Form("%1.4f", fImportance/fImportanceRef) << Endl;

   for ( UInt_t i=0; i<nvars; i++) {

      Log() << kINFO << "            ";
      sel    = fCut->GetSelector(i);
      valmin = fCut->GetCutMin(i);
      valmax = fCut->GetCutMax(i);
      //
      Log() << kINFO << Form("Cut %2d",i+1) << " : ";
      if (fCut->GetCutDoMin(i)) Log() << kINFO << Form("%10.3g",valmin) << " < ";
      else                      Log() << kINFO << "             ";
      Log() << kINFO << GetVarName(sel);
      if (fCut->GetCutDoMax(i)) Log() << kINFO << " < " << Form("%10.3g",valmax);
      else                      Log() << kINFO << "             ";
      Log() << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// extensive print function used to print info for the weight file

void TMVA::Rule::PrintRaw( std::ostream& os ) const
{
   Int_t dp = os.precision();
   const UInt_t nvars = fCut->GetNvars();
   os << "Parameters: "
      << std::setprecision(10)
      << fImportance << " "
      << fImportanceRef << " "
      << fCoefficient << " "
      << fSupport << " "
      << fSigma << " "
      << fNorm << " "
      << fSSB << " "
      << fSSBNeve << " "
      << std::endl;                                         \
   os << "N(cuts): " << nvars << std::endl; // mark end of nodes
   for ( UInt_t i=0; i<nvars; i++) {
      os << "Cut " << i << " : " << std::flush;
      os <<        fCut->GetSelector(i)
         << std::setprecision(10)
         << " " << fCut->GetCutMin(i)
         << " " << fCut->GetCutMax(i)
         << " " << (fCut->GetCutDoMin(i) ? "T":"F")
         << " " << (fCut->GetCutDoMax(i) ? "T":"F")
         << std::endl;
   }
   os << std::setprecision(dp);
}

////////////////////////////////////////////////////////////////////////////////

void* TMVA::Rule::AddXMLTo( void* parent ) const
{
   void* rule = gTools().AddChild( parent, "Rule" );
   const UInt_t nvars = fCut->GetNvars();

   gTools().AddAttr( rule, "Importance", fImportance    );
   gTools().AddAttr( rule, "Ref",        fImportanceRef );
   gTools().AddAttr( rule, "Coeff",      fCoefficient   );
   gTools().AddAttr( rule, "Support",    fSupport       );
   gTools().AddAttr( rule, "Sigma",      fSigma         );
   gTools().AddAttr( rule, "Norm",       fNorm          );
   gTools().AddAttr( rule, "SSB",        fSSB           );
   gTools().AddAttr( rule, "SSBNeve",    fSSBNeve       );
   gTools().AddAttr( rule, "Nvars",      nvars          );

   for (UInt_t i=0; i<nvars; i++) {
      void* cut = gTools().AddChild( rule, "Cut" );
      gTools().AddAttr( cut, "Selector", fCut->GetSelector(i) );
      gTools().AddAttr( cut, "Min",      fCut->GetCutMin(i) );
      gTools().AddAttr( cut, "Max",      fCut->GetCutMax(i) );
      gTools().AddAttr( cut, "DoMin",    (fCut->GetCutDoMin(i) ? "T":"F") );
      gTools().AddAttr( cut, "DoMax",    (fCut->GetCutDoMax(i) ? "T":"F") );
   }

   return rule;
}

////////////////////////////////////////////////////////////////////////////////
/// read rule from XML

void TMVA::Rule::ReadFromXML( void* wghtnode )
{
   TString nodeName = TString( gTools().GetName(wghtnode) );
   if (nodeName != "Rule") Log() << kFATAL << "<ReadFromXML> Unexpected node name: " << nodeName << Endl;

   gTools().ReadAttr( wghtnode, "Importance", fImportance    );
   gTools().ReadAttr( wghtnode, "Ref",        fImportanceRef );
   gTools().ReadAttr( wghtnode, "Coeff",      fCoefficient   );
   gTools().ReadAttr( wghtnode, "Support",    fSupport       );
   gTools().ReadAttr( wghtnode, "Sigma",      fSigma         );
   gTools().ReadAttr( wghtnode, "Norm",       fNorm          );
   gTools().ReadAttr( wghtnode, "SSB",        fSSB           );
   gTools().ReadAttr( wghtnode, "SSBNeve",    fSSBNeve       );

   UInt_t nvars;
   gTools().ReadAttr( wghtnode, "Nvars",      nvars          );
   if (fCut) delete fCut;
   fCut = new RuleCut();
   fCut->SetNvars( nvars );

   // read Cut
   void*    ch = gTools().GetChild( wghtnode );
   UInt_t   i = 0;
   UInt_t   ui;
   Double_t d;
   Char_t   c;
   while (ch) {
      gTools().ReadAttr( ch, "Selector", ui );
      fCut->SetSelector( i, ui );
      gTools().ReadAttr( ch, "Min",      d );
      fCut->SetCutMin  ( i, d );
      gTools().ReadAttr( ch, "Max",      d );
      fCut->SetCutMax  ( i, d );
      gTools().ReadAttr( ch, "DoMin",    c );
      fCut->SetCutDoMin( i, (c == 'T' ? kTRUE : kFALSE ) );
      gTools().ReadAttr( ch, "DoMax",    c );
      fCut->SetCutDoMax( i, (c == 'T' ? kTRUE : kFALSE ) );

      i++;
      ch = gTools().GetNextChild(ch);
   }

   // sanity check
   if (i != nvars) Log() << kFATAL << "<ReadFromXML> Mismatch in number of cuts: " << i << " != " << nvars << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// read function (format is the same as written by PrintRaw)

void TMVA::Rule::ReadRaw( std::istream& istr )
{
   TString dummy;
   UInt_t nvars;
   istr >> dummy
        >> fImportance
        >> fImportanceRef
        >> fCoefficient
        >> fSupport
        >> fSigma
        >> fNorm
        >> fSSB
        >> fSSBNeve;
   // coverity[tainted_data_argument]
   istr >> dummy >> nvars;
   Double_t cutmin,cutmax;
   UInt_t   sel,idum;
   Char_t   bA, bB;
   //
   if (fCut) delete fCut;
   fCut = new RuleCut();
   fCut->SetNvars( nvars );
   for ( UInt_t i=0; i<nvars; i++) {
      istr >> dummy >> idum; // get 'Node' and index
      istr >> dummy;         // get ':'
      istr >> sel >> cutmin >> cutmax >> bA >> bB;
      fCut->SetSelector(i,sel);
      fCut->SetCutMin(i,cutmin);
      fCut->SetCutMax(i,cutmax);
      fCut->SetCutDoMin(i,(bA=='T' ? kTRUE:kFALSE));
      fCut->SetCutDoMax(i,(bB=='T' ? kTRUE:kFALSE));
   }
}
