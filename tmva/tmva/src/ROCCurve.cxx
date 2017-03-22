// @(#)root/tmva $Id$
// Author: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer and Simon Pfreundschuh

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : ROCCurve                                                              *
 *                                                                                *
 * Description:                                                                   *
 *      This is class to compute ROC Integral (AUC)                               *
 *                                                                                *
 * Authors :                                                                      *
 *      Omar Zapata     <Omar.Zapata@cern.ch>    - UdeA/ITM Colombia              *
 *      Lorenzo Moneta  <Lorenzo.Moneta@cern.ch> - CERN, Switzerland              *
 *      Sergei Gleyzer  <Sergei.Gleyzer@cern.ch> - U of Florida & CERN            *
 *                                                                                *
 * Copyright (c) 2015:                                                            *
 *      CERN, Switzerland                                                         *
 *      UdeA/ITM, Colombia                                                        *
 *      U. of Florida, USA                                                        *
 **********************************************************************************/

/*! \class TMVA::ROCCurve
\ingroup TMVA

*/
#include "TMVA/Tools.h"
#include "TMVA/ROCCurve.h"
#include "TMVA/Config.h"
#include "TMVA/Version.h"
#include "TMVA/MsgLogger.h"
#include "TGraph.h"

#include<vector>
#include <cassert>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

TMVA::ROCCurve::ROCCurve(const std::vector<Float_t> & mva, const std::vector<Bool_t> & mvat) :
   fLogger ( new TMVA::MsgLogger("ROCCurve") ),fGraph(NULL)
{
   assert(mva.size() == mvat.size() );
   for(UInt_t i=0;i<mva.size();i++)
   {
      if(mvat[i] ) fMvaS.push_back(mva[i]);
      else fMvaB.push_back(mva[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::ROCCurve::~ROCCurve() {
   delete fLogger;
   if(fGraph) delete fGraph;
}

////////////////////////////////////////////////////////////////////////////////
/// Auxilary functions (Specificity and Sensitivity)

Float_t TMVA::ROCCurve::ComputeSen(Float_t threshold)
{
   Float_t true_positives = 0.0;
   Float_t false_negatives = 0.0;
   for (auto signal_res : fMvaS) {
       if (signal_res > threshold)
           ++true_positives;
       else
           ++false_negatives;
   }
   return true_positives + false_negatives == 0.0 ? 0.0 : true_positives / (true_positives + false_negatives);
}

Float_t TMVA::ROCCurve::ComputeSpe(Float_t threshold)
{
   Float_t true_negatives = 0.0;
   Float_t false_positives = 0.0;
   for (auto background_res : fMvaB) {
       if (background_res > threshold)
           ++false_positives;
       else
           ++true_negatives;
   }
   return true_negatives + false_positives == 0.0 ? 0.0 : true_negatives / (true_negatives + false_positives);
}

////////////////////////////////////////////////////////////////////////////////
/// ROC Integral (AUC)

Double_t TMVA::ROCCurve::GetROCIntegral()
{
   Float_t integral=0;
   int ndivisions = 40;
   fEpsilonSig.push_back(0);
   fEpsilonBgk.push_back(0);

   for (Float_t i=-1.0;i<1.0;i+=(1.0/ndivisions)) {
       fEpsilonSig.push_back(1 - ComputeSen(i));
       fEpsilonBgk.push_back(ComputeSpe(i));
   }

   fEpsilonSig.push_back(1.0);
   fEpsilonBgk.push_back(1.0);
   for (UInt_t i=0;i<fEpsilonSig.size()-1;i++) {
      integral += 0.5*(fEpsilonSig[i+1]-fEpsilonSig[i])*(fEpsilonBgk[i]+fEpsilonBgk[i+1]);
   }
   return integral;
}


////////////////////////////////////////////////////////////////////////////////

TGraph* TMVA::ROCCurve::GetROCCurve(const UInt_t points)
{
   const UInt_t ndivisions = points - 1;
   fEpsilonSig.resize(points);
   fEpsilonBgk.resize(points);
   // Fixed values.
   fEpsilonSig[0] = 0.0;
   fEpsilonSig[ndivisions] = 1.0;
   fEpsilonBgk[0] = 1.0;
   fEpsilonBgk[ndivisions] = 0.0;

   for (UInt_t i = 1; i < ndivisions; i++) {
      Float_t threshold = -1.0 + i * 2.0 / (Float_t) ndivisions;
      fEpsilonSig[ndivisions - i] = ComputeSen(threshold);
      fEpsilonBgk[ndivisions - i] = ComputeSpe(threshold);
   }
   if(!fGraph)    fGraph=new TGraph(fEpsilonSig.size(),&fEpsilonSig[0],&fEpsilonBgk[0]);
   return fGraph;
}
