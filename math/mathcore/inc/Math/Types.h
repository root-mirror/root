#ifndef ROOT_Math_VecTypes
#define ROOT_Math_VecTypes

#include "RConfigure.h"

#ifdef R__HAS_VECCORE

#if defined(R__HAS_VC)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wconditional-uninitialized"
#endif

#include <Vc/Vc>
#pragma GCC diagnostic pop
#endif

#include <VecCore/VecCore>

namespace ROOT {

namespace Internal {
   using ScalarBackend = vecCore::backend::Scalar;
#ifdef VECCORE_ENABLE_VC
   using VectorBackend = vecCore::backend::VcVector;
#else
   using VectorBackend = vecCore::backend::Scalar;
#endif
}
   using Float_v  = typename Internal::VectorBackend::Float_v;
   using Double_v = typename Internal::VectorBackend::Double_v;
   using Int_v    = typename Internal::VectorBackend::Int_v;
   using Int32_v  = typename Internal::VectorBackend::Int32_v;
   using UInt_v   = typename Internal::VectorBackend::UInt_v;
   using UInt32_v = typename Internal::VectorBackend::UInt32_v;
}

#else // R__HAS_VECCORE

// We do not have explicit vectorisation support enabled. Fall back to regular ROOT types.

#include "RtypesCore.h"

namespace ROOT {
   using Float_v  = Float_t;
   using Double_v = Double_t;
   using Int_v    = Int_t;
   using Int32_v  = Int_t; // FIXME: Should we introduce Int32_t in RtypesCore.h?
   using UInt_v   = UInt_t;
   using UInt32_v = UInt_t; // FIXME: Should we introduce UInt32_t in RtypesCore.h?
}
#endif

#endif // ROOT_Math_VecTypes
