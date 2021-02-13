// @(#)root/roostats:$Id$
// Author: Kyle Cranmer   28/07/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef RooStats_NumberCountingUtils
#define RooStats_NumberCountingUtils

//_________________________________________________
/*
BEGIN_HTML
<h2>NumberCountingUtils</h2>
<p>
These are  RooStats standalone utilities
that calculate the p-value or Z value (eg. significance in
1-sided Gaussian standard deviations) for a number counting experiment.
This is a hypothesis test between background only and signal-plus-background.
The background estimate has uncertainty derived from an auxiliary or sideband
measurement.
</p>
<p>
This is based on code and comments from Bob Cousins 
and on the following papers:
<p>
<ul>
<li>Evaluation of three methods for calculating statistical significance when incorporating a
systematic uncertainty into a test of the background-only hypothesis for a Poisson process<br />
Authors: Robert D. Cousins, James T. Linnemann, Jordan Tucker<br />
http://arxiv.org/abs/physics/0702156<br />
NIM  A 595 (2008) 480--501</li>

<li>
Statistical Challenges for Searches for New Physics at the LHC<br />
Authors: Kyle Cranmer<br />
http://arxiv.org/abs/physics/0511028
</li>
<li>
 Measures of Significance in HEP and Astrophysics<br />
 Authors: J. T. Linnemann<br />
 http://arxiv.org/abs/physics/0312059
</li>
</ul>
<p>
The problem is treated in a fully frequentist fashion by 
interpreting the relative background uncertainty as
being due to an auxiliary or sideband observation 
that is also Poisson distributed with only background.
Finally, one considers the test as a ratio of Poisson means
where an interval is well known based on the conditioning on the total
number of events and the binomial distribution.
</p>

<p>
In short, this is an exact frequentist solution to the problem of
a main measurement x distributed as a Poisson around s+b and a sideband or 
auxiliary measurement y distributed as a Poisson around tau*b.  Eg. 
</p>
END_HTML
BEGIN_LATEX
L(x,y|s,b,#tau) = Pois(x|s+b) Pois(y|#tau b)
END_LATEX
BEGIN_HTML
<pre>
Naming conventions:
Exp = Expected
Obs = Observed
P   = p-value
Z   = Z-value or significance in sigma (one-sided convention)
</pre>
END_HTML
*/
//

#include "RtypesCore.h"


namespace RooStats{

   namespace  NumberCountingUtils {

  
  // Expected P-value for s=0 in a ratio of Poisson means.  
  // Here the background and its uncertainty are provided directly and 
  // assumed to be from the double Poisson counting setup described in the 
  // BinomialWithTau functions.  
  // Normally one would know tau directly, but here it is determiend from
  // the background uncertainty.  This is not strictly correct, but a useful 
  // approximation.
     Double_t BinomialExpZ(Double_t sExp, Double_t bExp, Double_t fractionalBUncertainty);

  // See BinomialWithTauExpP
     Double_t BinomialWithTauExpZ(Double_t sExp, Double_t bExp, Double_t tau);   

  // See BinomialObsP
     Double_t BinomialObsZ(Double_t nObs, Double_t bExp, Double_t fractionalBUncertainty);

  // See BinomialWithTauObsP
     Double_t BinomialWithTauObsZ(Double_t nObs, Double_t bExp, Double_t tau);
     
  // See BinomialExpP
     Double_t BinomialExpP(Double_t sExp, Double_t bExp, Double_t fractionalBUncertainty);

  // Expected P-value for s=0 in a ratio of Poisson means.  
  // Based on two expectations, a main measurement that might have signal
  // and an auxiliarly measurement for the background that is signal free.
  // The expected background in the auxiliary measurement is a factor
  // tau larger than in the main measurement.
     Double_t BinomialWithTauExpP(Double_t sExp, Double_t bExp, Double_t tau);

  // P-value for s=0 in a ratio of Poisson means.  
  // Here the background and its uncertainty are provided directly and 
  // assumed to be from the double Poisson counting setup.  
  // Normally one would know tau directly, but here it is determiend from
  // the background uncertainty.  This is not strictly correct, but a useful 
  // approximation.
     Double_t BinomialObsP(Double_t nObs, Double_t, Double_t fractionalBUncertainty);

  // P-value for s=0 in a ratio of Poisson means.  
  // Based on two observations, a main measurement that might have signal
  // and an auxiliarly measurement for the background that is signal free.
  // The expected background in the auxiliary measurement is a factor
  // tau larger than in the main measurement.
     Double_t BinomialWithTauObsP(Double_t nObs, Double_t bExp, Double_t tau);
      

   }
}

#endif
