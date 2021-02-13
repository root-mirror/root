// @(#)root/roostats:$Id$
// Authors: Kevin Belasco        17/06/2009
// Authors: Kyle Cranmer         17/06/2009
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::MetropolisHastings
    \ingroup Roostats

This class uses the Metropolis-Hastings algorithm to construct a Markov Chain
of data points using Monte Carlo. In the main algorithm, new points in the
parameter space are proposed and then visited based on their relative
likelihoods.  This class can use any implementation of the ProposalFunction,
including non-symmetric proposal functions, to propose parameter points and
still maintain detailed balance when constructing the chain.



The "Likelihood" function that is sampled when deciding what steps to take in
the chain has been given a very generic implementation.  The user can create
any RooAbsReal based on the parameters and pass it to a MetropolisHastings
object with the method SetFunction(RooAbsReal&).  Be sure to tell
MetropolisHastings whether your RooAbsReal is on a (+/-) regular or log scale,
so that it knows what logic to use when sampling your RooAbsReal.  For example,
a common use is to sample from a -log(Likelihood) distribution (NLL), for which
the appropriate configuration calls are SetType(MetropolisHastings::kLog);
SetSign(MetropolisHastings::kNegative);
If you're using a traditional likelihood function:
SetType(MetropolisHastings::kRegular);  SetSign(MetropolisHastings::kPositive);
You must set these type and sign flags or MetropolisHastings will not construct
a MarkovChain.

Also note that in ConstructChain(), the values of the variables are randomized
uniformly over their intervals before construction of the MarkovChain begins.

*/

#include "RooStats/MetropolisHastings.h"

#include "RooStats/MarkovChain.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/ProposalFunction.h"

#include "Rtypes.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "TMath.h"
#include "TFile.h"

ClassImp(RooStats::MetropolisHastings);

using namespace RooFit;
using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

MetropolisHastings::MetropolisHastings()
{
   // default constructor
   fFunction = NULL;
   fPropFunc = NULL;
   fNumIters = 0;
   fNumBurnInSteps = 0;
   fSign = kSignUnset;
   fType = kTypeUnset;
}

////////////////////////////////////////////////////////////////////////////////

MetropolisHastings::MetropolisHastings(RooAbsReal& function, const RooArgSet& paramsOfInterest,
      ProposalFunction& proposalFunction, Int_t numIters)
{
   fFunction = &function;
   SetParameters(paramsOfInterest);
   SetProposalFunction(proposalFunction);
   fNumIters = numIters;
   fNumBurnInSteps = 0;
   fSign = kSignUnset;
   fType = kTypeUnset;
}

////////////////////////////////////////////////////////////////////////////////

MarkovChain* MetropolisHastings::ConstructChain()
{
   if (fParameters.getSize() == 0 || !fPropFunc || !fFunction) {
      coutE(Eval) << "Critical members unintialized: parameters, proposal " <<
                     " function, or (log) likelihood function" << endl;
         return NULL;
   }
   if (fSign == kSignUnset || fType == kTypeUnset) {
      coutE(Eval) << "Please set type and sign of your function using "
         << "MetropolisHastings::SetType() and MetropolisHastings::SetSign()" <<
         endl;
      return NULL;
   }

   if (fChainParams.getSize() == 0) fChainParams.add(fParameters);

   RooArgSet x;
   RooArgSet xPrime;
   x.addClone(fParameters);
   RandomizeCollection(x);
   xPrime.addClone(fParameters);
   RandomizeCollection(xPrime);

   MarkovChain* chain = new MarkovChain();
   // only the POI will be added to the chain
   chain->SetParameters(fChainParams);

   Int_t weight = 0;
   Double_t xL = 0.0, xPrimeL = 0.0, a = 0.0;

   // ibucur: i think the user should have the possibility to display all the message
   //    levels should they want to; maybe a setPrintLevel would be appropriate
   //    (maybe for the other classes that use this approach as well)?
   RooFit::MsgLevel oldMsgLevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

   // We will need to check if log-likelihood evaluation left an error status.
   // Now using faster eval error logging with CountErrors.
   if (fType == kLog) {
     RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
     //N.B: need to clear the count in case of previous errors !
     // the clear needs also to be done after calling setEvalErrorLoggingMode
     RooAbsReal::clearEvalErrorLog();
   }

   bool hadEvalError = true;

   Int_t i = 0;
   // get a good starting point for x
   // for fType == kLog, this means that fFunction->getVal() did not cause
   // an eval error
   // for fType == kRegular this means fFunction->getVal() != 0
   //
   // kbelasco: i < 1000 is sort of arbitrary, but way higher than the number of
   // steps we should have to take for any reasonable (log) likelihood function
   while (i < 1000 && hadEvalError) {
      RandomizeCollection(x);
      RooStats::SetParameters(&x, &fParameters);
      xL = fFunction->getVal();

      if (fType == kLog) {
         if (RooAbsReal::numEvalErrors() > 0) {
            RooAbsReal::clearEvalErrorLog();
            hadEvalError = true;
         } else
            hadEvalError = false;
      } else if (fType == kRegular) {
         if (xL == 0.0)
            hadEvalError = true;
         else
            hadEvalError = false;
      } else
         // for now the only 2 types are kLog and kRegular (won't get here)
         hadEvalError = false;
      ++i;
   }

   if(hadEvalError) {
      coutE(Eval) << "Problem finding a good starting point in " <<
                     "MetropolisHastings::ConstructChain() " << endl;
   }


   ooccoutP((TObject *)0, Generation) << "Metropolis-Hastings progress: ";

   // do main loop
   for (i = 0; i < fNumIters; i++) {
      // reset error handling flag
      hadEvalError = false;

      // print a dot every 1% of the chain construction
      if (i % (fNumIters / 100) == 0) ooccoutP((TObject*)0, Generation) << ".";

      fPropFunc->Propose(xPrime, x);

      RooStats::SetParameters(&xPrime, &fParameters);
      xPrimeL = fFunction->getVal();

      // check if log-likelihood for xprime had an error status
      if (fFunction->numEvalErrors() > 0 && fType == kLog) {
         xPrimeL = RooNumber::infinity();
         fFunction->clearEvalErrorLog();
         hadEvalError = true;
      }

      // why evaluate the last point again, can't we cache it?
      // kbelasco: commenting out lines below to add/test caching support
      //RooStats::SetParameters(&x, &fParameters);
      //xL = fFunction->getVal();

      if (fType == kLog) {
         if (fSign == kPositive)
            a = xL - xPrimeL;
         else
            a = xPrimeL - xL;
      }
      else
         a = xPrimeL / xL;
      //a = xL / xPrimeL;

      if (!hadEvalError && !fPropFunc->IsSymmetric(xPrime, x)) {
         Double_t xPrimePD = fPropFunc->GetProposalDensity(xPrime, x);
         Double_t xPD      = fPropFunc->GetProposalDensity(x, xPrime);
         if (fType == kRegular)
            a *= xPD / xPrimePD;
         else
            a += TMath::Log(xPrimePD) - TMath::Log(xPD);
      }

      if (!hadEvalError && ShouldTakeStep(a)) {
         // go to the proposed point xPrime

         // add the current point with the current weight
         if (weight != 0.0)
            chain->Add(x, CalcNLL(xL), (Double_t)weight);

         // reset the weight and go to xPrime
         weight = 1;
         RooStats::SetParameters(&xPrime, &x);
         xL = xPrimeL;
      } else {
         // stay at the current point
         weight++;
      }
   }

   // make sure to add the last point
   if (weight != 0.0)
      chain->Add(x, CalcNLL(xL), (Double_t)weight);
   ooccoutP((TObject *)0, Generation) << endl;

   RooMsgService::instance().setGlobalKillBelow(oldMsgLevel);

   Int_t numAccepted = chain->Size();
   coutI(Eval) << "Proposal acceptance rate: " <<
                   numAccepted/(Float_t)fNumIters * 100 << "%" << endl;
   coutI(Eval) << "Number of steps in chain: " << numAccepted << endl;

   //TFile chainDataFile("chainData.root", "recreate");
   //chain->GetDataSet()->Write();
   //chainDataFile.Close();

   return chain;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t MetropolisHastings::ShouldTakeStep(Double_t a)
{
   if ((fType == kLog && a <= 0.0) || (fType == kRegular && a >= 1.0)) {
      // The proposed point has a higher likelihood than the
      // current point, so we should go there
      return kTRUE;
   }
   else {
      // generate numbers on a log distribution to decide
      // whether to go to xPrime or stay at x
      //Double_t rand = fGen.Uniform(1.0);
      Double_t rand = RooRandom::uniform();
      if (fType == kLog) {
         rand = TMath::Log(rand);
         // kbelasco: should this be changed to just (-rand > a) for logical
         // consistency with below test when fType == kRegular?
         if (-1.0 * rand >= a)
            // we chose to go to the new proposed point
            // even though it has a lower likelihood than the current one
            return kTRUE;
      } else {
         // fType must be kRegular
         // kbelasco: ensure that we never visit a point where PDF == 0
         //if (rand <= a)
         if (rand < a)
            // we chose to go to the new proposed point
            // even though it has a lower likelihood than the current one
            return kTRUE;
      }
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////

Double_t MetropolisHastings::CalcNLL(Double_t xL)
{
   if (fType == kLog) {
      if (fSign == kNegative)
         return xL;
      else
         return -xL;
   } else {
      if (fSign == kPositive)
         return -1.0 * TMath::Log(xL);
      else
         return -1.0 * TMath::Log(-xL);
   }
}
