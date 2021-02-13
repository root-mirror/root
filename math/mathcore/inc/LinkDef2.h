// @(#)root/mathcore:$Id$

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace ROOT;
#pragma link C++ namespace ROOT::Math;

// for automatic loading
#ifdef MAKE_MAPS
#pragma link C++ class TMath;
//#pragma link C++ class ROOT::Math;
#endif

#pragma link C++ class std::vector<Double_t>::iterator;
#pragma link C++ class std::vector<Double_t>::const_iterator;
#pragma link C++ class std::vector<Double_t>::reverse_iterator;

#pragma link C++ global gRandom;

#pragma link C++ class TRandom+;
#pragma link C++ class TRandom1+;
#pragma link C++ class TRandom2+;
#pragma link C++ class TRandom3-;

#pragma link C++ class ROOT::Math::TRandomEngine+;
#pragma link C++ class ROOT::Math::LCGEngine+;
#pragma link C++ class ROOT::Math::MersenneTwisterEngine+;
#pragma link C++ class ROOT::Math::MixMaxEngine<240,0>+;
#pragma link C++ class ROOT::Math::MixMaxEngine<256,2>+;
#pragma link C++ class ROOT::Math::MixMaxEngine<17,1>+;
//#pragma link C++ class mixmax::mixmax_engine<240>+;
//#pragma link C++ class mixmax::mixmax_engine<256>+;
//#pragma link C++ class mixmax::mixmax_engine<17>+;
//#pragma link C++ struct mixmax::mixmax_engine<240>::rng_state_st+;
//#pragma link C++ struct mixmax::mixmax_engine<256>::rng_state_st+;
//#pragma link C++ struct mixmax::mixmax_engine<17>::rng_state_st+;
//#pragma link C++ struct mixmax::_Generator<ULong64_t,0,2305843009213693951>+;

#pragma link C++ class std::mersenne_twister_engine< uint_fast64_t, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29, 0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43, 6364136223846793005ULL >+;

#pragma link C++ class std::subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12>+;
#pragma link C++ class std::discard_block_engine<std::ranlux48_base, 389, 11>+;

#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<240,0>>+;
#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<256,0>>+;
#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<256,2>>+;
#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<256,4>>+;
#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<17,0>>+;
#pragma link C++ class TRandomGen<ROOT::Math::MixMaxEngine<17,1>>+;
#pragma link C++ class TRandomGen<ROOT::Math::RanluxppEngine2048>+;
#pragma link C++ class TRandomGen<ROOT::Math::StdEngine<std::mt19937_64>>+;
#pragma link C++ class TRandomGen<ROOT::Math::StdEngine<std::ranlux48>>+;


#pragma link C++ class ROOT::Math::StdRandomEngine+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::LCGEngine>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MersenneTwisterEngine>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<240,0>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<256,0>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<256,2>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<256,4>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<17,0>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<17,1>>+;
#pragma link C++ class ROOT::Math::Random<ROOT::Math::MixMaxEngine<17,2>>+;

// #pragma link C++ typedef ROOT::Math::RandomMT19937;
// #pragma link C++ typedef ROOT::Math::RandomMT64;
// #pragma link C++ typedef ROOT::Math::RandomRanlux48;
// #pragma link C++ typedef TRandomMixMax;
// #pragma link C++ typedef TRandomMixMax256;
// #pragma link C++ typedef TRandomMT64;
// #pragma link C++ typedef TRandomRanlux48;


// #pragma link C++ class TRandomNew3+;



#pragma link C++ class TStatistic+;


#pragma link C++ class TKDTree<Int_t, Double_t>+;
#pragma link C++ class TKDTree<Int_t, Float_t>+;
#pragma link C++ typedef TKDTreeID;
#pragma link C++ typedef TKDTreeIF;
#pragma link C++ class TKDTreeBinning-;


// ROOT::Math namespace
#pragma link C++ typedef ROOT::Math::IGenFunction;
#pragma link C++ typedef ROOT::Math::IMultiGenFunction;
#pragma link C++ typedef ROOT::Math::IGradFunction;
#pragma link C++ typedef ROOT::Math::IMultiGradFunction;

#pragma link C++ class ROOT::Math::IBaseFunctionOneDim+;
#pragma link C++ class ROOT::Math::IGradientOneDim+;
#pragma link C++ class ROOT::Math::IGradientFunctionOneDim+;
#pragma link C++ class ROOT::Math::IBaseParam+;

#pragma link C++ class ROOT::Math::IParametricFunctionOneDim+;
#pragma link C++ class ROOT::Math::IParametricGradFunctionOneDim+;

#pragma link C++ class ROOT::Math::IBaseFunctionMultiDim+;
#pragma link C++ class ROOT::Math::IGradientMultiDim+;
#pragma link C++ class ROOT::Math::IGradientFunctionMultiDim+;
#pragma link C++ class ROOT::Math::IParametricFunctionMultiDim+;
#pragma link C++ class ROOT::Math::IParametricGradFunctionMultiDim+;

#pragma link C++ class ROOT::Math::ParamFunctor+;
#pragma link C++ class ROOT::Math::Functor-;
#pragma link C++ class ROOT::Math::GradFunctor-;
#pragma link C++ class ROOT::Math::Functor1D-;
#pragma link C++ class ROOT::Math::GradFunctor1D-;

#pragma link C++ class ROOT::Math::Minimizer+;
#pragma link C++ class ROOT::Math::MinimizerOptions+;
#pragma link C++ class ROOT::Math::MinimTransformFunction-;
#pragma link C++ class ROOT::Math::MinimTransformVariable+;
#pragma link C++ class ROOT::Math::BasicMinimizer+;
#pragma link C++ class ROOT::Math::IntegratorOneDimOptions+;
#pragma link C++ class ROOT::Math::IntegratorMultiDimOptions+;
#pragma link C++ class ROOT::Math::BaseIntegratorOptions+;
#pragma link C++ class ROOT::Math::IOptions+;
#pragma link C++ class ROOT::Math::GenAlgoOptions+;
#pragma link C++ class ROOT::Math::IntegratorOneDim+;
#pragma link C++ class ROOT::Math::IntegratorMultiDim+;
#pragma link C++ class ROOT::Math::VirtualIntegrator+;
#pragma link C++ class ROOT::Math::VirtualIntegratorOneDim+;
#pragma link C++ class ROOT::Math::VirtualIntegratorMultiDim+;
#pragma link C++ class ROOT::Math::AdaptiveIntegratorMultiDim+;
#pragma link C++ typedef ROOT::Math::Integrator;

#pragma link C++ namespace ROOT::Math::IntegrationOneDim;
#pragma link C++ enum ROOT::Math::IntegrationOneDim::Type;
#pragma link C++ namespace ROOT::Math::IntegrationMultiDim;
// #pragma link C++ typedef ROOT::Math::IntegratorOneDim::Type;
// #pragma link C++ typedef ROOT::Math::IntegratorMultiDim::Type;


#pragma link C++ class ROOT::Math::BasicFitMethodFunction<ROOT::Math::IBaseFunctionMultiDim>+;
#ifndef _WIN32
#pragma link C++ class ROOT::Math::BasicFitMethodFunction<ROOT::Math::IGradientFunctionMultiDim>+;
#else
// problem due to virtual inheritance
#pragma link C++ class ROOT::Math::BasicFitMethodFunction<ROOT::Math::IGradientFunctionMultiDim>-;
#endif
// typedef's
#pragma link C++ typedef ROOT::Math::FitMethodFunction;
#pragma link C++ typedef ROOT::Math::FitMethodGradFunction;


#pragma link C++ class ROOT::Math::Factory+;

#pragma link C++ class ROOT::Math::GaussIntegrator+;
#pragma link C++ class ROOT::Math::GaussLegendreIntegrator+;
#pragma link C++ class ROOT::Math::RichardsonDerivator+;

#pragma link C++ class ROOT::Math::RootFinder+;
#pragma link C++ class ROOT::Math::IRootFinderMethod+;
#pragma link C++ class ROOT::Math::BrentRootFinder+;
#pragma link C++ class ROOT::Math::IMinimizer1D+;
#pragma link C++ class ROOT::Math::BrentMinimizer1D+;

#pragma link C++ class ROOT::Math::DistSampler+;
#pragma link C++ class ROOT::Math::DistSamplerOptions+;
#pragma link C++ class ROOT::Math::GoFTest+;
#pragma link C++ class std::vector<std::vector<double> >+;

#pragma link C++ class ROOT::Math::Delaunay2D+;


#pragma link C++ class ROOT::Math::TDataPoint<1,Float_t>+;
#pragma link C++ typedef ROOT::Math::TDataPoint1F;
#pragma link C++ class ROOT::Math::TDataPoint<1,Double_t>+;
#pragma link C++ typedef ROOT::Math::TDataPoint1F;
#pragma link C++ typedef ROOT::Math::TDataPoint1D;
#pragma link C++ class  ROOT::Math::TDataPointN<Double_t>+;
#pragma link C++ class  ROOT::Math::TDataPointN<Float_t>+;
//
//N.B. use old streamer (do not use +) for KDTree class because it will not work on Windows
// to work one would need to change the internal classes from private to public
#pragma link C++ class ROOT::Math::KDTree<ROOT::Math::TDataPoint1D>;



#include "LinkDef_Func.h"

#endif
