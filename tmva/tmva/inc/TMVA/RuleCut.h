// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Rule                                                                  *
 *                                                                                *
 * Description:                                                                   *
 *      A class describing a 'rule cut'                                           *
 *                                                                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U.                                                             *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/
#ifndef ROOT_TMVA_RuleCut
#define ROOT_TMVA_RuleCut

#include "TMVA/Event.h"

#include <vector>

namespace TMVA {

   class Node;
   class MsgLogger;

   class RuleCut {

   public:

      // main constructor
      RuleCut( const std::vector< const TMVA::Node * > & nodes );

      // copy constructor
   RuleCut( const RuleCut & other ) : fLogger(0) { Copy( other ); }

      // empty constructor
      RuleCut();

      // destructor
      virtual ~RuleCut();

      // evaluate an event
      inline Bool_t EvalEvent( const Event &eve );

      // get cut range for a given selector
      Bool_t GetCutRange(Int_t sel,Double_t &rmin, Double_t &rmax, Bool_t &dormin, Bool_t &dormax) const;

      // number of cuts
      UInt_t GetNcuts() const;

      // set members
      inline void SetNvars( UInt_t nc );
      void SetNeve( Double_t n )                   { fCutNeve     = n;   }
      void SetPurity( Double_t ssb )               { fPurity      = ssb; }
      void SetSelector( Int_t i, UInt_t s )        { fSelector[i] = s; }
      void SetCutMin( Int_t i, Double_t v )        { fCutMin[i]   = v; }
      void SetCutMax( Int_t i, Double_t v )        { fCutMax[i]   = v; }
      void SetCutDoMin( Int_t i, Bool_t v )        { fCutDoMin[i] = v; }
      void SetCutDoMax( Int_t i, Bool_t v )        { fCutDoMax[i] = v; }

      // accessors
      UInt_t   GetNvars()              const { return fSelector.size(); }
      UInt_t   GetSelector(Int_t is)   const { return fSelector[is]; }
      Double_t GetCutMin(Int_t is)     const { return fCutMin[is]; }
      Double_t GetCutMax(Int_t is)     const { return fCutMax[is]; }
      Char_t   GetCutDoMin(Int_t is)   const { return fCutDoMin[is]; }
      Char_t   GetCutDoMax(Int_t is)   const { return fCutDoMax[is]; }
      Double_t GetCutNeve()            const { return fCutNeve; }
      Double_t GetPurity()             const { return fPurity; }

   private:
      // copy
      inline void Copy( const RuleCut & other);

      // make the cuts from the array of nodes
      void MakeCuts( const std::vector< const TMVA::Node * > & nodes );

      std::vector<UInt_t>   fSelector; // array of selectors (expressions)
      std::vector<Double_t> fCutMin;   // array of lower limits
      std::vector<Double_t> fCutMax;   // array of upper limits
      std::vector<Char_t>   fCutDoMin; // array of usage flags for lower limits <--- stores boolean
      std::vector<Char_t>   fCutDoMax; // array of usage flags for upper limits <--- stores boolean
      Double_t              fCutNeve;  // N(events) after cut (possibly weighted)
      Double_t              fPurity;  // S/(S+B) on training data


      mutable MsgLogger*    fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }
   };
}

//_______________________________________________________________________
inline void TMVA::RuleCut::Copy( const TMVA::RuleCut & other )
{
   // copy from another
   if (&other != this) {
      for (UInt_t ns=0; ns<other.GetNvars(); ns++) {
         fSelector.push_back( other.GetSelector(ns) );
         fCutMin.push_back( other.GetCutMin(ns) );
         fCutMax.push_back( other.GetCutMax(ns) );
         fCutDoMin.push_back( other.GetCutDoMin(ns) );
         fCutDoMax.push_back( other.GetCutDoMax(ns) );
      }
      fCutNeve = other.GetCutNeve();
      fPurity = other.GetPurity();
   }
}

//_______________________________________________________________________
inline Bool_t TMVA::RuleCut::EvalEvent( const Event &eve )
{
   // evaluate event using the cut

   // Loop over all cuts
   Int_t    sel;
   Double_t val;
   Bool_t done=kFALSE;
   Bool_t minOK, cutOK=kFALSE;
   UInt_t nc=0;
   while (!done) {
      sel = fSelector[nc];
      val = eve.GetValue(sel);
      minOK = (fCutDoMin[nc] ? (val>fCutMin[nc]):kTRUE); // min cut ok
      cutOK = (minOK ? ((fCutDoMax[nc] ? (val<fCutMax[nc]):kTRUE)) : kFALSE); // cut ok
      nc++;
      done = ((!cutOK) || (nc==fSelector.size())); // done if 
   }
   //   return ( cutOK ? 1.0: 0.0 );
   return cutOK;
}

//_______________________________________________________________________
inline void TMVA::RuleCut::SetNvars( UInt_t nc )
{
   // set the number of cuts
   fSelector.clear();
   fCutMin.clear();
   fCutMax.clear();
   fCutDoMin.clear();
   fCutDoMax.clear();
   //
   fSelector.resize(nc);
   fCutMin.resize(nc);
   fCutMax.resize(nc);
   fCutDoMin.resize(nc);
   fCutDoMax.resize(nc);
}

#endif
