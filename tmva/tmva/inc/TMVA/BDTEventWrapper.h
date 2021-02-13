
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : BDTEventWrapper                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *                                                                                *
 *                                                                                *
 * Author: Doug Schouten (dschoute@sfu.ca)                                        *
 *                                                                                *
 * Copyright (c) 2007:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Texas at Austin, USA                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_BDTEventWrapper
#define ROOT_TMVA_BDTEventWrapper

#include "RtypesCore.h"
#include "Event.h"
#include "ThreadLocalStorage.h"

namespace TMVA {

   class BDTEventWrapper{

   public:

      BDTEventWrapper( const Event* );
      ~BDTEventWrapper();

      // Require '<' operator to use std::sort algorithms on collection of Events
      Bool_t operator <( const BDTEventWrapper& other ) const;

      // Set the accumulated weight, for sorted signal/background events
      /**
       * @param type - true for signal, false for background
       * @param weight - the total weight
       */
      void SetCumulativeWeight( Bool_t type, Double_t weight );

      // Get the accumulated weight
      /**
       * @param type - true for signal, false for background
       * @return the cumulative weight for sorted signal/background events
       */
      Double_t GetCumulativeWeight( Bool_t type ) const;

      // Set the index of the variable to compare on
      /**
       * @param iVar - index of the variable in fEvent to use
       */
      inline static void SetVarIndex( Int_t iVar ) { if (iVar >= 0) GetVarIndex() = iVar; }

      // Return the value of variable fVarIndex for this event
      /**
       * @return value of variable fVarIndex for this event
       */
      inline Double_t GetVal() const { return fEvent->GetValue(GetVarIndex()); }
      const Event* operator*() const { return fEvent; }

      inline Double_t GetVal(Int_t var) const { return fEvent->GetValue(var); }
   private:

      // This is a workaround for OSx where static thread_local data members are
      // not supported. The C++ solution would indeed be the following:
      static Int_t& GetVarIndex(){TTHREAD_TLS(Int_t) fVarIndex(0); return fVarIndex;}; // index of the variable to sort on

      const Event* fEvent;     // pointer to the event

      Double_t     fBkgWeight; // cumulative background weight for splitting
      Double_t     fSigWeight; // same for the signal weights
   };
}

inline Bool_t TMVA::BDTEventWrapper::operator<( const BDTEventWrapper& other ) const
{
   return GetVal() < other.GetVal();
}

#endif
