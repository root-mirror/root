// @(#)root/gl:$Id$
// Author: Timur Pocheptsov  2009
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TKDEFGT
#define ROOT_TKDEFGT

#include <vector>

#include "RtypesCore.h"

//KDE - kenrel density estimator,
//based on fast Gauss transform.

class TGL5DDataSet;

class TKDEFGT {
private:
   //KDE-related stuff.
   std::vector<Double_t>    fXC;     //Centers.
   std::vector<Double_t>    fWeights;//Weights.
   std::vector<UInt_t>      fIndxc;  //Internal data.
   std::vector<Double_t>    fA_K;    //Polynomial coefficient (pd x K)
   std::vector<UInt_t>      fIndx;   //Internal data.
   std::vector<UInt_t>      fXhead;  //Internal data.
   std::vector<UInt_t>      fXboxsz; //Internal data.
   std::vector<Double_t>    fDistC;  //Internal data.
   std::vector<Double_t>    fC_K;    //Internal data.
   std::vector<UInt_t>      fCinds;  //Internal data.

   mutable std::vector<UInt_t>   fHeads; //Internal data.
   mutable std::vector<Double_t> fDx;    //Internal data.
   mutable std::vector<Double_t> fProds; //Internal data.

   UInt_t                  fDim;        //Number of dimensions.
   UInt_t                  fP;          //Order of trancation.
   UInt_t                  fK;          //Number of centers.
   Double_t                fSigma;      //Noise Standard deviation of the kernel (default sigma = 1)
   UInt_t                  fPD;         //nchoosek(fP + fDim - 1, fDim);
   Bool_t                  fModelValid; //Check, if coefficients are ok.

public:
   TKDEFGT();

   virtual ~TKDEFGT();

   //Generic version.
   //"sources" must be a vector of packed coordinates, if you have
   //dim == 3, vector will be [xyz|xyz|...].
   void BuildModel(const std::vector<Double_t> &sources, Double_t sigma = 1.,
                   UInt_t dim = 3, UInt_t p = 8, UInt_t k = 0);
   //Special version for data from TTree.
   void BuildModel(const TGL5DDataSet *sources, Double_t sigma = 1.,
                   UInt_t p = 8, UInt_t k = 0);
   //"targets" is a vector of packed coordinates, the same as above in BuildModel:
   //[xyz|xyz|xyz...] for dim == 3.
   void Predict(const std::vector<Double_t> &targets, std::vector<Double_t> &densities,
                Double_t e)const;

private:
   //Generic version.
   void Kcenter(const std::vector<double> &x);
   //Special version for data sources from TTree.
   void Kcenter(const TGL5DDataSet *sources);

   void Compute_C_k();
   //Generic version.
   void Compute_A_k(const std::vector<Double_t> &x);
   //Version for TTree.
   void Compute_A_k(const TGL5DDataSet *sources);

   TKDEFGT(const TKDEFGT &rhs);
   TKDEFGT &operator = (const TKDEFGT &rhs);
};

#endif
