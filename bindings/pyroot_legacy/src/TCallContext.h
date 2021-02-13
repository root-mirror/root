#ifndef PYROOT_TCALLCONTEXT_H
#define PYROOT_TCALLCONTEXT_H

// Standard
#include <vector>


namespace PyROOT {

// general place holder for function parameters
   struct TParameter {
      union Value {
         Bool_t       fBool;
         Short_t      fShort;
         UShort_t     fUShort;
         Int_t        fInt;
         UInt_t       fUInt;
         Long_t       fLong;
         ULong_t      fULong;
         Long64_t     fLongLong;
         ULong64_t    fULongLong;
         Float_t      fFloat;
         Double_t     fDouble;
         LongDouble_t fLongDouble;
         void*        fVoidp;
      } fValue;
      void* fRef;
      char  fTypeCode;
   };

// extra call information
   struct TCallContext {
      TCallContext( std::vector< TParameter >::size_type sz = 0 ) : fArgs( sz ), fFlags( 0 ) {}

      enum ECallFlags {
         kNone           =    0,
         kIsSorted       =    1,   // if method overload priority determined
         kIsCreator      =    2,   // if method creates python-owned objects
         kIsConstructor  =    4,   // if method is a C++ constructor
         kUseHeuristics  =    8,   // if method applies heuristics memory policy
         kUseStrict      =   16,   // if method applies strict memory policy
         kManageSmartPtr =   32,   // if executor should manage smart pointers
         kReleaseGIL     =   64,   // if method should release the GIL
         kFast           =  128,   // if method should NOT handle signals
         kSafe           =  256    // if method should return on signals
      };

   // memory handling
      static ECallFlags sMemoryPolicy;
      static Bool_t SetMemoryPolicy( ECallFlags e );

   // signal safety
      static ECallFlags sSignalPolicy;
      static Bool_t SetSignalPolicy( ECallFlags e );

   // payload
      std::vector< TParameter > fArgs;
      UInt_t fFlags;
   };

   inline Bool_t IsSorted( UInt_t flags ) {
      return flags & TCallContext::kIsSorted;
   }

   inline Bool_t IsCreator( UInt_t flags ) {
      return flags & TCallContext::kIsCreator;
   }

   inline Bool_t IsConstructor( UInt_t flags ) {
      return flags & TCallContext::kIsConstructor;
   }

   inline Bool_t ManagesSmartPtr( TCallContext* ctxt ) {
      return ctxt->fFlags & TCallContext::kManageSmartPtr;
   }

   inline Bool_t ReleasesGIL( UInt_t flags ) {
      return flags & TCallContext::kReleaseGIL;
   }

   inline Bool_t ReleasesGIL( TCallContext* ctxt ) {
      return ctxt ? (ctxt->fFlags & TCallContext::kReleaseGIL) : kFALSE;
   }

   inline Bool_t UseStrictOwnership( TCallContext* ctxt ) {
      if ( ctxt && (ctxt->fFlags & TCallContext::kUseStrict) )
         return kTRUE;
      if ( ctxt && (ctxt->fFlags & TCallContext::kUseHeuristics) )
         return kFALSE;

      return TCallContext::sMemoryPolicy == TCallContext::kUseStrict;
   }

} // namespace PyROOT

#endif // !PYROOT_TCALLCONTEXT_H
