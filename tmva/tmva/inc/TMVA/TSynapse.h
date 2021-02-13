// @(#)root/tmva $Id$
// Author: Matt Jachowski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::TSynapse                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Synapse class for use in derivatives of MethodANNBase                     *
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

#ifndef ROOT_TMVA_TSynapse
#define ROOT_TMVA_TSynapse

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSynapse                                                             //
//                                                                      //
// Synapse used by derivatives of MethodANNBase                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

namespace TMVA {

   class TNeuron;
   class MsgLogger;

   class TSynapse : public TObject {

   public:

      TSynapse();
      virtual ~TSynapse();

      // set the weight of the synapse
      void SetWeight(Double_t weight);

      // get the weight of the synapse
      Double_t GetWeight()                 { return fWeight;         }

      // set the learning rate
      void SetLearningRate(Double_t rate)  { fLearnRate = rate;      }

      // get the learning rate
      Double_t GetLearningRate()           { return fLearnRate;      }

      // decay the learning rate
      void DecayLearningRate(Double_t rate){ fLearnRate *= (1-rate); }

      // set the pre-neuron
      void SetPreNeuron(TNeuron* pre)      { fPreNeuron = pre;       }

      // set the post-neuron
      void SetPostNeuron(TNeuron* post)    { fPostNeuron = post;     }

      // get the weighted output of the pre-neuron
      Double_t GetWeightedValue();

      // get the weighted error field of the post-neuron
      Double_t GetWeightedDelta();

      // force the synapse to adjust its weight according to its error field
      void AdjustWeight();

      // calculate the error field of the synapse
      void CalculateDelta();

      // initialize the error field of the synapse to 0
      void InitDelta()           { fDelta = 0.0; fCount = 0; }

      void SetDEDw(Double_t DEDw)              { fDEDw = DEDw;           }
      Double_t GetDEDw()                       { return fDEDw;           }
      Double_t GetDelta()                      { return fDelta;          }

   private:

      Double_t fWeight;            // weight of the synapse
      Double_t fLearnRate;         // learning rate parameter
      Double_t fDelta;             // local error field
      Double_t fDEDw;              // sum of deltas
      Int_t    fCount;             // number of updates contributing to error field
      TNeuron* fPreNeuron;         // pointer to pre-neuron
      TNeuron* fPostNeuron;        // pointer to post-neuron

      MsgLogger& Log() const;

      ClassDef(TSynapse,0); // Synapse class used by MethodANNBase and derivatives
   };

} // namespace TMVA

#endif
