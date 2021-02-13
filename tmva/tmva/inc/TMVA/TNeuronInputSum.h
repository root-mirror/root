// @(#)root/tmva $Id$
// Author: Matt Jachowski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::TNeuronInputSum                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      TNeuron input calculator -- calculates the weighted sum of inputs.        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Matt Jachowski  <jachowski@stanford.edu> - Stanford University, USA       *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/
 

#ifndef ROOT_TMVA_TNeuronInputSum
#define ROOT_TMVA_TNeuronInputSum

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TNeuronInputSum                                                      //
//                                                                      //
// TNeuron input calculator -- calculates the weighted sum of inputs    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/TNeuronInput.h"
#include "TMVA/TNeuron.h"

namespace TMVA {

   class TNeuronInputSum : public TNeuronInput {
    
   public:

      TNeuronInputSum() {}
      virtual ~TNeuronInputSum() {}

      // calculate input value for neuron
      Double_t GetInput( const TNeuron* neuron ) const {
         if (neuron->IsInputNeuron()) return 0;
         Double_t result = 0;
         Int_t npl = neuron->NumPreLinks();
         for (Int_t i=0; i < npl; i++) {
            result += neuron->PreLinkAt(i)->GetWeightedValue();
         }
         return result;
      }

      // name of class
      TString GetName() { return "Sum of weighted activations"; }

      ClassDef(TNeuronInputSum,0); // Calculates weighted sum of neuron inputs
   };

} // namespace TMVA

#endif
