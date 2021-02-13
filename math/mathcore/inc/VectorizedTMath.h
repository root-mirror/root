#ifndef ROOT_VectorizedTMath
#define ROOT_VectorizedTMath

#include "RtypesCore.h"
#include "Math/Types.h"
#include "TMath.h"

#if defined(R__HAS_VECCORE) && defined(R__HAS_VC)

namespace TMath {
::ROOT::Double_v Log2(::ROOT::Double_v &x);
::ROOT::Double_v BreitWigner(::ROOT::Double_v &x, Double_t mean = 0, Double_t gamma = 1);
::ROOT::Double_v Gaus(::ROOT::Double_v &x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE);
::ROOT::Double_v LaplaceDist(::ROOT::Double_v &x, Double_t alpha = 0, Double_t beta = 1);
::ROOT::Double_v LaplaceDistI(::ROOT::Double_v &x, Double_t alpha = 0, Double_t beta = 1);
::ROOT::Double_v Freq(::ROOT::Double_v &x);
::ROOT::Double_v BesselI0_Split_More(::ROOT::Double_v &ax);
::ROOT::Double_v BesselI0_Split_Less(::ROOT::Double_v &x);
::ROOT::Double_v BesselI0(::ROOT::Double_v &x);
::ROOT::Double_v BesselI1_Split_More(::ROOT::Double_v &ax, ::ROOT::Double_v &x);
::ROOT::Double_v BesselI1_Split_Less(::ROOT::Double_v &x);
::ROOT::Double_v BesselI1(::ROOT::Double_v &x);
::ROOT::Double_v BesselJ0_Split1_More(::ROOT::Double_v &ax);
::ROOT::Double_v BesselJ0_Split1_Less(::ROOT::Double_v &x);
::ROOT::Double_v BesselJ0(::ROOT::Double_v &x);
::ROOT::Double_v BesselJ1_Split1_More(::ROOT::Double_v &ax, ::ROOT::Double_v &x);
::ROOT::Double_v BesselJ1_Split1_Less(::ROOT::Double_v &x);
::ROOT::Double_v BesselJ1(::ROOT::Double_v &x);
} // namespace TMath

#endif // VECCORE and VC exist check

#endif
