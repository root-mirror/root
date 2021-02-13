#ifndef ROOT_VECTYPE
#define ROOT_VECTYPE

#include "Math/Vector2D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4Dfwd.h"
#include "TrackMathCore.h"
#include <string>

template <class V>
struct VecType {
   static std::string name() { return "MathCoreVector"; }
};
template <>
struct VecType<ROOT::Math::XYVector> {
   static std::string name() { return "XYVector"; }
   static std::string name32() { return "ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<Double32_t> >"; }
};
template <>
struct VecType<ROOT::Math::Polar2DVector> {
   static std::string name() { return "Polar2DVector"; }
};
template <>
struct VecType<ROOT::Math::XYZVector> {
   static std::string name() { return "XYZVector"; }
   static std::string name32() { return "ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double32_t> >"; }
};
template <>
struct VecType<ROOT::Math::Polar3DVector> {
   static std::string name() { return "Polar3DVector"; }
};
template <>
struct VecType<ROOT::Math::RhoEtaPhiVector> {
   static std::string name() { return "RhoEtaPhiVector"; }
};
template <>
struct VecType<ROOT::Math::RhoZPhiVector> {
   static std::string name() { return "RhoZPhiVector"; }
};
template <>
struct VecType<ROOT::Math::PxPyPzEVector> {
   static std::string name() { return "PxPyPzEVector"; }
   static std::string name32() { return "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double32_t> >"; }
};
template <>
struct VecType<ROOT::Math::PtEtaPhiEVector> {
   static std::string name() { return "PtEtaPhiEVector"; }
};
template <>
struct VecType<ROOT::Math::PtEtaPhiMVector> {
   static std::string name() { return "PtEtaPhiMVector"; }
};
template <>
struct VecType<ROOT::Math::PxPyPzMVector> {
   static std::string name() { return "PxPyPzMVector"; }
};
template <>
struct VecType<ROOT::Math::SVector<double, 3>> {
   static std::string name() { return "SVector3"; }
   static std::string name32() { return "ROOT::Math::SVector<Double32_t,3>"; }
};
template <>
struct VecType<ROOT::Math::SVector<double, 4>> {
   static std::string name() { return "SVector4"; }
   static std::string name32() { return "ROOT::Math::SVector<Double32_t,4>"; }
};
template <>
struct VecType<TrackD> {
   static std::string name() { return "TrackD"; }
};
template <>
struct VecType<TrackD32> {
   static std::string name() { return "TrackD32"; }
};
template <>
struct VecType<TrackErrD> {
   static std::string name() { return "TrackErrD"; }
};
template <>
struct VecType<TrackErrD32> {
   static std::string name() { return "TrackErrD32"; }
};
template <>
struct VecType<VecTrack<TrackD>> {
   static std::string name() { return "VecTrackD"; }
};
template <>
struct VecType<VecTrack<TrackErrD>> {
   static std::string name() { return "VecTrackErrD"; }
};

#endif // ROOT_VECTYPE
