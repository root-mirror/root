// @(#)root/tmva $Id$ 
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MinuitWrapper                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Wrapper around MINUIT                                                    *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <peter.speckmayer@cern.ch> - CERN, Switzerland           *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MinuitWrapper
#define ROOT_TMVA_MinuitWrapper

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MinuitWrapper                                                        //
//                                                                      //
// Wrapper around MINUIT                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMinuit.h"
#include "TMVA/IFitterTarget.h"
#include <vector>

class TMinuit;

namespace TMVA {

   class IFitterTarget;

   class MinuitWrapper : public TMinuit {

   public:

      MinuitWrapper( IFitterTarget& target, Int_t maxpar);
      virtual ~MinuitWrapper() {}

      Int_t Eval(Int_t, Double_t*, Double_t&, Double_t*, Int_t);
      void SetFitterTarget( IFitterTarget& target ) { fFitterTarget = target; }

      Int_t ExecuteCommand(const char *command, Double_t *args, Int_t nargs);
      void  Clear(Option_t * = 0);
      Int_t GetStats    (Double_t &amin, Double_t &edm, Double_t &errdef, Int_t &nvpar, Int_t &nparx);
      Int_t GetErrors   (Int_t ipar, Double_t &eplus, Double_t &eminus, Double_t &eparab, Double_t &globcc);
      Int_t SetParameter(Int_t ipar,const char *parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh);
      TObject *Clone(char const*) const;
      
   private:
      
      IFitterTarget&        fFitterTarget; // fitter Target
      std::vector<Double_t> fParameters;   // vector holding the current parameters 
      Int_t                 fNumPar;       // number of parameters
      
      ClassDef(MinuitWrapper,0); // Wrapper around TMinuit
   };

} // namespace TMVA

#endif


