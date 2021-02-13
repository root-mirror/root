// $Id$
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : LDA                                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Local LDA method used by MethodKNN to compute MVA value.                  *
 *      This is experimental code under development. This class computes          *
 *      parameters of signal and background PDFs using Gaussian approximation.    *
 *                                                                                *
 * Author:                                                                        *
 *      John Alison John.Alison@cern.ch - University of Pennsylvania, USA         *
 *                                                                                *
 * Copyright (c) 2007:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      University of Pennsylvania, USA                                           *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_LDA
#define ROOT_TMVA_LDA


// C/C++
#include <map>
#include <vector>

// ROOT
#include "RtypesCore.h"
#include "TMatrixFfwd.h"

typedef std::vector<std::vector<Float_t> >  LDAEvents;

namespace TMVA {

   class MsgLogger;

   class LDA {

   public:

      LDA(Float_t tolerence = 1.0e-5, Bool_t debug = false);
      ~LDA();

      // Signal probability with Gaussian approximation
      Float_t GetProb(const std::vector<Float_t>& x, Int_t k);

      // Log likelihood function with Gaussian approximation
      Float_t GetLogLikelihood(const std::vector<Float_t>& x, Int_t k);

      // Create LDA matrix using local events found by knn method
      void Initialize(const LDAEvents& inputSignal, const LDAEvents& inputBackground);

   private:

      // Probability value using Gaussian approximation
      Float_t FSub(const std::vector<Float_t>& x, Int_t k);

      MsgLogger& Log() const { return *fLogger; }

   private:

      // data members
      Float_t       fTolerence;                    // documentation!
      UInt_t        fNumParams;                    // documentation!
      std::map<Int_t, std::vector<Float_t> > fMu;  // documentation!
      TMatrixF*     fSigma;                        // documentation!
      TMatrixF*     fSigmaInverse;                 // documentation!
      std::map<Int_t, Float_t> fEventFraction;     // documentation!
      Bool_t        fDebug;                        // documentation!

      mutable MsgLogger *fLogger;                  // message logging service
   };
}
#endif
