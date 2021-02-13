// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006 - 2019

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//==============================================================================
// LinkDef.h - REve objects and services.
//==============================================================================

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// REveVector
#pragma link C++ class   ROOT::Experimental::REveVectorT<Float_t>+;
#pragma link C++ class   ROOT::Experimental::REveVectorT<Double_t>+;
#pragma link C++ typedef ROOT::Experimental::REveVector;
#pragma link C++ typedef ROOT::Experimental::REveVectorF;
#pragma link C++ typedef ROOT::Experimental::REveVectorD;

#pragma link C++ class   ROOT::Experimental::REveVector4T<Float_t>+;
#pragma link C++ class   ROOT::Experimental::REveVector4T<Double_t>+;
#pragma link C++ typedef ROOT::Experimental::REveVector4;
#pragma link C++ typedef ROOT::Experimental::REveVector4F;
#pragma link C++ typedef ROOT::Experimental::REveVector4D;

#pragma link C++ class   ROOT::Experimental::REveVector2T<Float_t>+;
#pragma link C++ class   ROOT::Experimental::REveVector2T<Double_t>+;
#pragma link C++ typedef ROOT::Experimental::REveVector2;
#pragma link C++ typedef ROOT::Experimental::REveVector2F;
#pragma link C++ typedef ROOT::Experimental::REveVector2D;

// Operators for REveVectorXT<Float_t>
#pragma link C++ function operator+(const ROOT::Experimental::REveVectorT<Float_t>&, const ROOT::Experimental::REveVectorT<Float_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVectorT<Float_t>&, const ROOT::Experimental::REveVectorT<Float_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVectorT<Float_t>&, Float_t);
#pragma link C++ function operator*(Float_t, const ROOT::Experimental::REveVectorT<Float_t>&);
#pragma link C++ function operator+(const ROOT::Experimental::REveVector4T<Float_t>&, const ROOT::Experimental::REveVector4T<Float_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVector4T<Float_t>&, const ROOT::Experimental::REveVector4T<Float_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVector4T<Float_t>&, Float_t);
#pragma link C++ function operator*(Float_t, const ROOT::Experimental::REveVector4T<Float_t>&);
#pragma link C++ function operator+(const ROOT::Experimental::REveVector2T<Float_t>&, const ROOT::Experimental::REveVector2T<Float_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVector2T<Float_t>&, const ROOT::Experimental::REveVector2T<Float_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVector2T<Float_t>&, Float_t);
#pragma link C++ function operator*(Float_t, const ROOT::Experimental::REveVector2T<Float_t>&);
// Operators for REveVectorXT<Double_t>
#pragma link C++ function operator+(const ROOT::Experimental::REveVectorT<Double_t>&, const ROOT::Experimental::REveVectorT<Double_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVectorT<Double_t>&, const ROOT::Experimental::REveVectorT<Double_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVectorT<Double_t>&, Double_t);
#pragma link C++ function operator*(Double_t, const ROOT::Experimental::REveVectorT<Double_t>&);
#pragma link C++ function operator+(const ROOT::Experimental::REveVector4T<Double_t>&, const ROOT::Experimental::REveVector4T<Double_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVector4T<Double_t>&, const ROOT::Experimental::REveVector4T<Double_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVector4T<Double_t>&, Double_t);
#pragma link C++ function operator*(Double_t, const ROOT::Experimental::REveVector4T<Double_t>&);
#pragma link C++ function operator+(const ROOT::Experimental::REveVector2T<Double_t>&, const ROOT::Experimental::REveVector2T<Double_t>&);
#pragma link C++ function operator-(const ROOT::Experimental::REveVector2T<Double_t>&, const ROOT::Experimental::REveVector2T<Double_t>&);
#pragma link C++ function operator*(const ROOT::Experimental::REveVector2T<Double_t>&, Double_t);
#pragma link C++ function operator*(Double_t, const ROOT::Experimental::REveVector2T<Double_t>&);

// REvePathMark
#pragma link C++ class   ROOT::Experimental::REvePathMarkT<Float_t>+;
#pragma link C++ class   ROOT::Experimental::REvePathMarkT<Double_t>+;
#pragma link C++ typedef ROOT::Experimental::REvePathMark;
#pragma link C++ typedef ROOT::Experimental::REvePathMarkF;
#pragma link C++ typedef ROOT::Experimental::REvePathMarkD;

// REveTrans
#pragma link C++ class ROOT::Experimental::REveTrans-;

// REveUtil
#pragma link C++ class ROOT::Experimental::REveUtil+;
#pragma link C++ class ROOT::Experimental::REveException+;
#pragma link C++ class ROOT::Experimental::REveGeoManagerHolder+;
#pragma link C++ class ROOT::Experimental::REveRefCnt+;
#pragma link C++ class ROOT::Experimental::REveRefBackPtr+;

// REveManager
#pragma link C++ class ROOT::Experimental::REveManager+;
#pragma link C++ global ROOT::Experimental::gEve;
#pragma link C++ class ROOT::Experimental::REveManager::RRedrawDisabler+;
#pragma link C++ class ROOT::Experimental::REveManager::RExceptionHandler+;

// REveVSD
#pragma link C++ class ROOT::Experimental::REveMCTrack+;
#pragma link C++ class ROOT::Experimental::REveHit+;
#pragma link C++ class ROOT::Experimental::REveCluster+;

#pragma link C++ class   ROOT::Experimental::REveRecTrackT<Float_t>+;
#pragma link C++ class   ROOT::Experimental::REveRecTrackT<Double_t>+;
#pragma link C++ typedef ROOT::Experimental::REveRecTrack;
#pragma link C++ typedef ROOT::Experimental::REveRecTrackF;
#pragma link C++ typedef ROOT::Experimental::REveRecTrackD;

#pragma link C++ class ROOT::Experimental::REveRecKink+;
#pragma link C++ class ROOT::Experimental::REveRecV0+;
#pragma link C++ class ROOT::Experimental::REveRecCascade+;
#pragma link C++ class ROOT::Experimental::REveMCRecCrossRef+;

#pragma link C++ class ROOT::Experimental::REveVSD+;

// REveTreeTools
#pragma link C++ class ROOT::Experimental::REveSelectorToEventList+;
#pragma link C++ class ROOT::Experimental::REvePointSelectorConsumer+;
#pragma link C++ class ROOT::Experimental::REvePointSelector+;

#pragma link C++ class ROOT::Experimental::REveRenderData+;

// REveElement
#pragma link C++ class ROOT::Experimental::REveElement+;
#pragma link C++ class ROOT::Experimental::REveAunt+;
#pragma link C++ class ROOT::Experimental::REveAuntAsList+;

#pragma link C++ class std::list<ROOT::Experimental::REveElement*>+;
#pragma link C++ typedef ROOT::Experimental::REveElement::List_t;

#pragma link C++ class std::set<ROOT::Experimental::REveElement*>+;
#pragma link C++ typedef ROOT::Experimental::REveElement::Set_t;

// REveCompound
#pragma link C++ class ROOT::Experimental::REveCompound+;
#pragma link C++ class ROOT::Experimental::REveCompoundProjected+;

// REveSelection
#pragma link C++ class ROOT::Experimental::REveSelection+;
#pragma link C++ class ROOT::Experimental::REveSecondarySelectable+;

// 3D Viewers and Scenes
#pragma link C++ class ROOT::Experimental::REveScene+;
#pragma link C++ class ROOT::Experimental::REveSceneList+;
#pragma link C++ class ROOT::Experimental::REveSceneInfo+;
#pragma link C++ class ROOT::Experimental::REveViewer+;
#pragma link C++ class ROOT::Experimental::REveViewerList+;

// Data classes
#pragma link C++ class ROOT::Experimental::REveViewContext+;
#pragma link C++ class ROOT::Experimental::REveDataCollection+;
#pragma link C++ class ROOT::Experimental::REveDataItem+;
#pragma link C++ class ROOT::Experimental::REveDataItemList+;
#pragma link C++ class ROOT::Experimental::REveDataProxyBuilderBase+;
#pragma link C++ class ROOT::Experimental::REveDataSimpleProxyBuilder+;
#pragma link C++ class ROOT::Experimental::REveDataTable+;
#pragma link C++ class ROOT::Experimental::REveDataColumn+;

// Projections / non-linear transformations
#pragma link C++ class ROOT::Experimental::REveProjectable+;
#pragma link C++ class ROOT::Experimental::REveProjected+;
#pragma link C++ class ROOT::Experimental::REveProjection+;
#pragma link C++ class ROOT::Experimental::REveProjection::PreScaleEntry_t+;
#pragma link C++ class std::vector<ROOT::Experimental::REveProjection::PreScaleEntry_t>;
#pragma link C++ typedef ROOT::Experimental::REveProjection::vPreScale_t;
#pragma link C++ class ROOT::Experimental::REveRhoZProjection+;
#pragma link C++ class ROOT::Experimental::REveRPhiProjection+;
#pragma link C++ class ROOT::Experimental::REve3DProjection+;

#pragma link C++ class ROOT::Experimental::REveProjectionManager+;
// #pragma link C++ class ROOT::Experimental::REveProjectionAxes+;

// Generic configuration
// #pragma link C++ class REveParamList;
// #pragma link C++ class REveParamList::FloatConfig_t+;
// #pragma link C++ class REveParamList::IntConfig_t+;
// #pragma link C++ class REveParamList::BoolConfig_t+;

// REveTrack
#pragma link C++ class ROOT::Experimental::REveTrack+;
#pragma link C++ class ROOT::Experimental::REveTrackList+;
#pragma link C++ class ROOT::Experimental::REveTrackProjected+;
#pragma link C++ class ROOT::Experimental::REveTrackListProjected+;

// REveTrackPropagator
#pragma link C++ class ROOT::Experimental::REveTrackPropagator+;
#pragma link C++ class ROOT::Experimental::REveMagField+;
#pragma link C++ class ROOT::Experimental::REveMagFieldConst+;
#pragma link C++ class ROOT::Experimental::REveMagFieldDuo+;

// REvePointSet
#pragma link C++ class ROOT::Experimental::REvePointSet+;
#pragma link C++ class ROOT::Experimental::REvePointSetArray+;
#pragma link C++ class ROOT::Experimental::REvePointSetProjected+;

// REveLine
#pragma link C++ class ROOT::Experimental::REveLine+;
#pragma link C++ class ROOT::Experimental::REveLineProjected+;

// Shapes
#pragma link C++ class ROOT::Experimental::REveShape+;
#pragma link C++ class ROOT::Experimental::REvePolygonSetProjected+;
#pragma link C++ class ROOT::Experimental::REveGeoShape+;
#pragma link C++ class ROOT::Experimental::REveGeoShapeProjected+;
#pragma link C++ class ROOT::Experimental::REveGeoShapeExtract+;
#pragma link C++ class ROOT::Experimental::REveGeoPolyShape+;

// Not yet ported
// #pragma link C++ class ROOT::Experimental::REveGeoNode+;
// #pragma link C++ class ROOT::Experimental::REveGeoTopNode+;

// REveCalo
#pragma link C++ class ROOT::Experimental::REveCaloData+;
#pragma link C++ class ROOT::Experimental::REveCaloData::SliceInfo_t+;
#pragma link C++ class ROOT::Experimental::REveCaloDataVec;
#pragma link C++ class ROOT::Experimental::REveCaloDataHist+;
#pragma link C++ class ROOT::Experimental::REveCaloViz+;
#pragma link C++ class ROOT::Experimental::REveCalo3D+;
#pragma link C++ class ROOT::Experimental::REveCalo2D+;
#pragma link C++ class ROOT::Experimental::REveCaloLego+;

// REveJetCone
#pragma link C++ class ROOT::Experimental::REveJetCone+;
#pragma link C++ class ROOT::Experimental::REveJetConeProjected+;

// Ellipse
#pragma link C++ class ROOT::Experimental::REveEllipsoid+;
#pragma link C++ class ROOT::Experimental::REveEllipsoidProjected+;

// REveLineSet
#pragma link C++ class ROOT::Experimental::REveStraightLineSet+;
#pragma link C++ class ROOT::Experimental::REveStraightLineSetProjected+;

// Digits

// REveBox-es
#pragma link C++ class ROOT::Experimental::REveFrameBox+;
#pragma link C++ class ROOT::Experimental::REveBox+;
#pragma link C++ class ROOT::Experimental::REveBoxProjected+;

// Digits
#pragma link C++ class ROOT::Experimental::REveRGBAPalette+;
#pragma link C++ class ROOT::Experimental::REveDigitSet+;
#pragma link C++ class ROOT::Experimental::REveBoxSet+;

// REveChunkManager
#pragma link C++ class ROOT::Experimental::REveChunkManager+;
#pragma link C++ class ROOT::Experimental::REveChunkManager::iterator;

// Tables
#pragma link C++ class ROOT::Experimental::REveTableViewInfo;

// Geometry viewer
#pragma link C++ class ROOT::Experimental::REveGeomNodeBase+;
#pragma link C++ class ROOT::Experimental::REveGeomNode+;
#pragma link C++ class ROOT::Experimental::REveGeomVisible+;
#pragma link C++ class ROOT::Experimental::RGeomRenderInfo+;
#pragma link C++ class ROOT::Experimental::RGeomRawRenderInfo+;
#pragma link C++ class ROOT::Experimental::RGeomShapeRenderInfo;
#pragma link C++ class ROOT::Experimental::REveGeomDescription+;
#pragma link C++ class ROOT::Experimental::REveGeomDrawing+;
#pragma link C++ class ROOT::Experimental::REveGeomRequest+;
#pragma link C++ class ROOT::Experimental::REveGeomNodeInfo+;
#pragma link C++ class ROOT::Experimental::REveGeomConfig+;

#endif
