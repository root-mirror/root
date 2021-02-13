// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RegressionVariance                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description: Calculate the separation criteria used in regression              *
 *                                                                                *
 *          There are two things: the Separation Index, and the Separation Gain   *
 *          Separation Index:                                                     *
 *          Measure of the "Variance" of a sample.                                *
 *                                                                                *
 *          Separation Gain:                                                      *
 *          the measure of how the quality of separation of the sample increases  *
 *          by splitting the sample e.g. into a "left-node" and a "right-node"    *
 *          (N * Index_parent) - (N_left * Index_left) - (N_right * Index_right)  *
 *          this is then the quality criteria which is optimized for when trying  *
 *          to increase the information in the system (making the best selection  *
 *                                                                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      Heidelberg U., Germany                                                    *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/
#include "TMath.h"
#include "TMVA/RegressionVariance.h"

ClassImp(TMVA::RegressionVariance);

/*! \class TMVA::RegressionVariance
\ingroup TMVA
Calculate the "SeparationGain" for Regression analysis
separation criteria used in various training algorithms

There are two things: the Separation Index, and the Separation Gain
Separation Index:
Measure of the "Variance" of a sample.

Separation Gain:
the measure of how the quality of separation of the sample increases
by splitting the sample e.g. into a "left-node" and a "right-node"
(N * Index_parent) - (N_left * Index_left) - (N_right * Index_right)
this is then the quality criteria which is optimized for when trying
to increase the information in the system (making the best selection
*/

////////////////////////////////////////////////////////////////////////////////
/// Separation Gain:
/// the measure of how the quality of separation of the sample increases
/// by splitting the sample e.g. into a "left-node" and a "right-node"
/// (N * Index_parent) - (N_left * Index_left) - (N_right * Index_right)
/// this is then the quality criteria which is optimized for when trying
/// to increase the information in the system
/// for the Regression: as the "Gain is maximised", the RMS (sqrt(variance))
/// which is used as a "separation" index should be as small as possible.
/// the "figure of merit" here has to be -(rms left+rms-right) or 1/rms...

Double_t TMVA::RegressionVariance::GetSeparationGain(const Double_t nLeft,
                                                     const Double_t targetLeft, const Double_t target2Left,
                                                     const Double_t nTot,
                                                     const Double_t targetTot, const Double_t target2Tot)
{

   if  ( nTot==nLeft || nLeft==0 ) return 0.;

   Double_t parentIndex = nTot * this->GetSeparationIndex(nTot,targetTot,target2Tot);
   Double_t leftIndex   = ( (nTot - nLeft) * this->GetSeparationIndex(nTot-nLeft,targetTot-targetLeft,target2Tot-target2Left) );
   Double_t rightIndex  =    nLeft * this->GetSeparationIndex(nLeft,targetLeft,target2Left);

   //  return 1/ (leftIndex + rightIndex);
   return (parentIndex - leftIndex - rightIndex)/(parentIndex);
}

////////////////////////////////////////////////////////////////////////////////
/// Separation Index:  a simple Variance

Double_t TMVA::RegressionVariance::GetSeparationIndex(const Double_t n,
                                                      const Double_t target, const Double_t target2)
{
   //   return TMath::Sqrt(( target2 - target*target/n) / n);
   return ( target2 - target*target/n) / n;

}



