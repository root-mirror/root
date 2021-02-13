// @(#)root/mlp:$Id$
// Author: Christophe.Delaere@cern.ch   21/08/2002

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSynapse

This is a simple weighted bidirectional connection between
two neurons.
A network is built connecting two neurons by a synapse.
In addition to the value, the synapse can return the DeDw

*/

#include "TSynapse.h"
#include "TNeuron.h"

ClassImp(TSynapse);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TSynapse::TSynapse()
{
   fpre    = 0;
   fpost   = 0;
   fweight = 1;
   fDEDw   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor that connects two neurons

TSynapse::TSynapse(TNeuron * pre, TNeuron * post, Double_t w)
{
   fpre    = pre;
   fpost   = post;
   fweight = w;
   fDEDw   = 0;
   pre->AddPost(this);
   post->AddPre(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the pre-neuron

void TSynapse::SetPre(TNeuron * pre)
{
   if (fpre) {
      Error("SetPre","this synapse is already assigned to a pre-neuron.");
      return;
   }
   fpre = pre;
   pre->AddPost(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the post-neuron

void TSynapse::SetPost(TNeuron * post)
{
   if (fpost) {
      Error("SetPost","this synapse is already assigned to a post-neuron.");
      return;
   }
   fpost = post;
   post->AddPre(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the value: weighted input

Double_t TSynapse::GetValue() const
{
   if (fpre)
      return (fweight * fpre->GetValue());
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the derivative of the error wrt the synapse weight.

Double_t TSynapse::GetDeDw() const
{
   if (!(fpre && fpost))
      return 0;
   return (fpre->GetValue() * fpost->GetDeDw());
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the weight of the synapse.
/// This weight is the multiplying factor applied on the
/// output of a neuron in the linear combination given as input
/// of another neuron.

void TSynapse::SetWeight(Double_t w)
{
   fweight = w;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the derivative of the total error wrt the synapse weight

void TSynapse::SetDEDw(Double_t in)
{
   fDEDw = in;
}


