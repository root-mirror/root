#ifndef ROOT_TEfficiency_cxx
#define ROOT_TEfficiency_cxx

//standard header
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <cassert>

//ROOT headers
#include "Math/DistFuncMathCore.h"
#include "TBinomialEfficiencyFitter.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TError.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/WrappedFunction.h"

//custom headers
#include "TEfficiency.h"

// file with extra class for FC method
#include "TEfficiencyHelper.h"

//default values
const Double_t kDefBetaAlpha = 1;
const Double_t kDefBetaBeta = 1;
const Double_t kDefConfLevel = 0.682689492137; // 1 sigma
const TEfficiency::EStatOption kDefStatOpt = TEfficiency::kFCP;
const Double_t kDefWeight = 1;

ClassImp(TEfficiency);

////////////////////////////////////////////////////////////////////////////////
/** \class TEfficiency
    \ingroup Hist
    \brief Class to handle efficiency histograms

## I. Overview
This class handles the calculation of efficiencies and their uncertainties. It
provides several statistical methods for calculating frequentist and Bayesian
confidence intervals as well as a function for combining several efficiencies.

Efficiencies have a lot of applications and meanings but in principle, they can
be described by the fraction of good/passed events k out of sample containing
N events. One is usually interested in the dependency of the efficiency on other
(binned) variables. The number of passed and total events is therefore stored
internally in two histograms (TEfficiency::fTotalHistogram and TEfficiency::fPassedHistogram).
Then the efficiency, as well as its upper and lower error, can be calculated for each bin
individually.

As the efficiency can be regarded as a parameter of a binomial distribution, the
number of passed and total events must always be integer numbers. Therefore a
filling with weights is not possible. However, you can assign a global weight to each
TEfficiency object (TEfficiency::SetWeight).
It is necessary to create one TEfficiency object
for each weight if you investigate a process involving different weights. This
procedure needs more effort but enables you to re-use the filled object in cases
where you want to change one or more weights. This would not be possible if all
events with different weights were filled in the same histogram.

## II. Creating a TEfficiency object
If you start a new analysis, it is highly recommended to use the TEfficiency class
from the beginning. You can then use one of the constructors for fixed or
variable bin size and your desired dimension. These constructors append the
created TEfficiency object to the current directory. So it will be written
automatically to a file during the next TFile::Write command.

Example: create a two-dimensional TEfficiency object with
- name = "eff"
- title = "my efficiency"
- axis titles: x, y and LaTeX-formatted epsilon as a label for Z axis
- 10 bins with constant bin width (= 1) along X axis starting at 0 (lower edge
  from the first bin) up to 10 (upper edge of last bin)
- 20 bins with constant bin width (= 0.5) along Y axis starting at -5 (lower
  edge from the first bin) up to 5 (upper edge of last bin)

        TEfficiency* pEff = new TEfficiency("eff","my efficiency;x;y;#epsilon",10,0,10,20,-5,5);

If you already have two histograms filled with the number of passed and total
events, you will use the constructor TEfficiency(const TH1& passed,const TH1& total)
to construct the TEfficiency object. The histograms "passed" and "total" have
to fulfill the conditions mentioned in TEfficiency::CheckConsistency, otherwise the construction will fail.
As the histograms already exist, the new TEfficiency is by default **not** attached
to the current directory to avoid duplication of data. If you want to store the
new object anyway, you can either write it directly by calling TObject::Write or attach it to a directory using TEfficiency::SetDirectory.
This also applies to TEfficiency objects created by the copy constructor TEfficiency::TEfficiency(const TEfficiency& rEff).


### Example 1

~~~~~~~~~~~~~~~{.cpp}
TEfficiency* pEff = 0;
TFile* pFile = new TFile("myfile.root","recreate");

//h_pass and h_total are valid and consistent histograms
if(TEfficiency::CheckConsistency(h_pass,h_total))
{
  pEff = new TEfficiency(h_pass,h_total);
  // this will write the TEfficiency object to "myfile.root"
  // AND pEff will be attached to the current directory
  pEff->Write();
}
~~~~~~~~~~~~~~~

### Example 2

~~~~~~~~~~~~~~~{.cpp}
TEfficiency* pEff = 0;
TFile* pFile = new TFile("myfile.root","recreate");

//h_pass and h_total are valid and consistent histograms
if(TEfficiency::CheckConsistency(h_pass,h_total))
{
  pEff = new TEfficiency(h_pass,h_total);
  //this will attach the TEfficiency object to the current directory
  pEff->SetDirectory(gDirectory);
  //now all objects in gDirectory will be written to "myfile.root"
  pFile->Write();
}
~~~~~~~~~~~~~~~

In case you already have two filled histograms and you only want to
plot them as a graph, you should rather use TGraphAsymmErrors::TGraphAsymmErrors(const TH1* pass,const TH1* total,Option_t* opt)
to create a graph object.

## III. Filling with events
You can fill the TEfficiency object by calling the TEfficiency::Fill(Bool_t bPassed,Double_t x,Double_t y,Double_t z) method.
The "bPassed" boolean flag indicates whether the current event is good
(both histograms are filled) or not (only TEfficiency::fTotalHistogram is filled).
The x, y and z variables determine the bin which is filled. For lower dimensions, the z- or even the y-value may be omitted.

Begin_Macro(source)
{
   //canvas only needed for this documentation
   TCanvas* c1 = new TCanvas("example","",600,400);
   c1->SetFillStyle(1001);
   c1->SetFillColor(kWhite);

   //create one-dimensional TEfficiency object with fixed bin size
   TEfficiency* pEff = new TEfficiency("eff","my efficiency;x;#epsilon",20,0,10);
   TRandom3 rand3;

   bool bPassed;
   double x;
   for(int i=0; i<10000; ++i)
   {
      //simulate events with variable under investigation
      x = rand3.Uniform(10);
      //check selection: bPassed = DoesEventPassSelection(x)
      bPassed = rand3.Rndm() < TMath::Gaus(x,5,4);
      pEff->Fill(bPassed,x);
   }

   pEff->Draw("AP");

   //only for this documentation
   return c1;
}
End_Macro

You can also set the number of passed or total events for a bin directly by
using the TEfficiency::SetPassedEvents or TEfficiency::SetTotalEvents method.

## IV. Statistic options
The calculation of the estimated efficiency depends on the chosen statistic
option. Let k denotes the number of passed events and N the number of total
events.

###Frequentist methods
The expectation value of the number of passed events is given by the true
efficiency times the total number of events. One can estimate the efficiency
by replacing the expected number of passed events by the observed number of
passed events.

\f[
      k = \epsilon \times N    \Rightarrow    \hat{\varepsilon} = \frac{k}{N}
\f]

### Bayesian methods
In Bayesian statistics a likelihood-function (how probable is it to get the
observed data assuming a true efficiency) and a prior probability (what is the
probability that a certain true efficiency is actually realised) are used to
determine a posterior probability by using Bayes theorem. At the moment, only
beta distributions (have 2 free parameters) are supported as prior
probabilities.

\f{eqnarray*}{
 P(\epsilon | k ; N) &=& \frac{1}{norm} \times P(k | \epsilon ; N) \times Prior(\epsilon) \\
 P(k | \epsilon ; N) &=& Binomial(N,k) \times \epsilon^{k} \times (1 - \epsilon)^{N - k} ...\  binomial\ distribution \\
 Prior(\epsilon) &=& \frac{1}{B(\alpha,\beta)} \times \epsilon ^{\alpha - 1} \times (1 - \epsilon)^{\beta - 1} \equiv Beta(\epsilon; \alpha,\beta) \\
 \Rightarrow P(\epsilon | k ; N) &=& \frac{1}{norm'} \times \epsilon^{k + \alpha - 1} \times (1 - \epsilon)^{N - k + \beta - 1} \equiv Beta(\epsilon; k + \alpha, N - k + \beta)
\f}

By default the expectation value of this posterior distribution is used as an estimator for the efficiency:

\f[
      \hat{\varepsilon} = \frac{k + \alpha}{N + \alpha + \beta}
\f]

Optionally the mode can also be used as a value for the estimated efficiency. This can be done by calling
SetBit(kPosteriorMode) or TEfficiency::SetPosteriorMode. In this case, the estimated efficiency is:

\f[
       \hat{\varepsilon} = \frac{k + \alpha -1}{N + \alpha + \beta - 2}
\f]

In the case of a uniform prior distribution, B(x,1,1), the posterior mode is k/n, equivalent to the frequentist
estimate (the maximum likelihood value).

The statistic options also specify which confidence interval is used for calculating
the uncertainties of the efficiency. The following properties define the error
calculation:
- **fConfLevel:** desired confidence level: 0 < fConfLevel < 1 (TEfficiency::GetConfidenceLevel / TEfficiency::SetConfidenceLevel)
- **fStatisticOption** defines which method is used to calculate the boundaries of the confidence interval (TEfficiency::SetStatisticOption)
- **fBeta_alpha, fBeta_beta:** parameters for the prior distribution which is only used in the bayesian case (TEfficiency::GetBetaAlpha / TEfficiency::GetBetaBeta / TEfficiency::SetBetaAlpha / TEfficiency::SetBetaBeta)
- **kIsBayesian:** flag whether bayesian statistics are used or not (TEfficiency::UsesBayesianStat)
- **kShortestInterval:** flag whether shortest interval (instead of central one) are used in case of Bayesian statistics  (TEfficiency::UsesShortestInterval). Normally shortest interval should be used in combination with the mode (see TEfficiency::UsesPosteriorMode)
- **fWeight:** global weight for this TEfficiency object which is used during combining or merging with other TEfficiency objects(TEfficiency::GetWeight / TEfficiency::SetWeight)

In the following table, the implemented confidence intervals are listed
with their corresponding statistic option. For more details on the calculation,
please have a look at the mentioned functions.


| name             | statistic option | function            | kIsBayesian | parameters |
|------------------|------------------|---------------------|-------------|------------|
| Clopper-Pearson | kFCP          | TEfficiency::ClopperPearson |false |total events, passed events, confidence level |
| normal approximation | kFNormal | TEfficiency::Normal | false | total events, passed events, confidence level |
| Wilson | kFWilson | TEfficiency::Wilson | false | total events, passed events, confidence level |
| Agresti-Coull | kFAC | TEfficiency::AgrestiCoull | false | total events, passed events. confidence level |
| Feldman-Cousins | kFFC | TEfficiency::FeldmanCousins | false | total events, passed events, confidence level |
| Mid-P Lancaster | kMidP | TEfficiency::MidPInterval | false | total events, passed events, confidence level |
| Jeffrey | kBJeffrey | TEfficiency::Bayesian | true | total events, passed events, confidence level, fBeta_alpha = 0.5, fBeta_beta = 0.5 |
| Uniform prior | kBUniform |TEfficiency::Bayesian | true |total events, passed events, confidence level, fBeta_alpha = 1, fBeta_beta = 1 |
| custom prior | kBBayesian |TEfficiency::Bayesian | true |total events, passed events, confidence level, fBeta_alpha, fBeta_beta |

The following example demonstrates the effect of different statistic options and
confidence levels.

Begin_Macro(source)
{
   //canvas only needed for the documentation
   TCanvas* c1 = new TCanvas("c1","",600,400);
   c1->Divide(2);
   c1->SetFillStyle(1001);
   c1->SetFillColor(kWhite);

   //create one-dimensional TEfficiency object with fixed bin size
   TEfficiency* pEff = new TEfficiency("eff","different confidence levels;x;#epsilon",20,0,10);
   TRandom3 rand3;

   bool bPassed;
   double x;
   for(int i=0; i<1000; ++i)
   {
      //simulate events with variable under investigation
      x = rand3.Uniform(10);
      //check selection: bPassed = DoesEventPassSelection(x)
      bPassed = rand3.Rndm() < TMath::Gaus(x,5,4);
      pEff->Fill(bPassed,x);
   }

   //set style attributes
   pEff->SetFillStyle(3004);
   pEff->SetFillColor(kRed);

   //copy current TEfficiency object and set new confidence level
   TEfficiency* pCopy = new TEfficiency(*pEff);
   pCopy->SetConfidenceLevel(0.90);

   //set style attributes
   pCopy->SetFillStyle(3005);
   pCopy->SetFillColor(kBlue);

   c1->cd(1);

   //add legend
   TLegend* leg1 = new TLegend(0.3,0.1,0.7,0.5);
   leg1->AddEntry(pEff,"95%","F");
   leg1->AddEntry(pCopy,"68.3%","F");

   pEff->Draw("A4");
   pCopy->Draw("same4");
   leg1->Draw("same");

   //use same confidence level but different statistic methods
   TEfficiency* pEff2 = new TEfficiency(*pEff);
   TEfficiency* pCopy2 = new TEfficiency(*pEff);

   pEff2->SetStatisticOption(TEfficiency::kFNormal);
   pCopy2->SetStatisticOption(TEfficiency::kFAC);

   pEff2->SetTitle("different statistic options;x;#epsilon");

   //set style attributes
   pCopy2->SetFillStyle(3005);
   pCopy2->SetFillColor(kBlue);

   c1->cd(2);

   //add legend
   TLegend* leg2 = new TLegend(0.3,0.1,0.7,0.5);
   leg2->AddEntry(pEff2,"kFNormal","F");
   leg2->AddEntry(pCopy2,"kFAC","F");

   pEff2->Draw("a4");
   pCopy2->Draw("same4");
   leg2->Draw("same");

   //only for this documentation
   c1->cd(0);
   return c1;
}
End_Macro

The prior probability of the efficiency in Bayesian statistics can be given
in terms of a beta distribution. The beta distribution has two positive shape
parameters. The resulting priors for different combinations of these shape
parameters are shown in the plot below.

Begin_Macro(source)
{
      //canvas only needed for the documentation
      TCanvas* c1 = new TCanvas("c1","",600,400);
      c1->SetFillStyle(1001);
      c1->SetFillColor(kWhite);

      //create different beta distributions
      TF1* f1 = new TF1("f1","TMath::BetaDist(x,1,1)",0,1);
      f1->SetLineColor(kBlue);
      TF1* f2 = new TF1("f2","TMath::BetaDist(x,0.5,0.5)",0,1);
      f2->SetLineColor(kRed);
      TF1* f3 = new TF1("f3","TMath::BetaDist(x,1,5)",0,1);
      f3->SetLineColor(kGreen+3);
      f3->SetTitle("Beta distributions as priors;#epsilon;P(#epsilon)");
      TF1* f4 = new TF1("f4","TMath::BetaDist(x,4,3)",0,1);
      f4->SetLineColor(kViolet);

      //add legend
      TLegend* leg = new TLegend(0.25,0.5,0.85,0.89);
      leg->SetFillColor(kWhite);
      leg->SetFillStyle(1001);
      leg->AddEntry(f1,"a=1, b=1","L");
      leg->AddEntry(f2,"a=0.5, b=0.5","L");
      leg->AddEntry(f3,"a=1, b=5","L");
      leg->AddEntry(f4,"a=4, b=3","L");

      f3->Draw();
      f1->Draw("same");
      f2->Draw("Same");
      f4->Draw("same");
      leg->Draw("same");

      //only for this documentation
      return c1;
}
End_Macro


## IV.1 Coverage probabilities for different methods
The following pictures illustrate the actual coverage probability for the
different values of the true efficiency and the total number of events when a
confidence level of 95% is desired.

\image html normal95.gif  "Normal Approximation"


\image html wilson95.gif  "Wilson"


\image html ac95.gif  "Agresti Coull"


\image html cp95.gif  "Clopper Pearson"


\image html uni95.gif  "Bayesian with Uniform Prior"


\image html jeffrey95.gif  "Bayesian with Jeffrey Prior"

The average (over all possible true efficiencies) coverage probability for
different number of total events is shown in the next picture.
\image html av_cov.png "Average Coverage"

## V. Merging and combining TEfficiency objects
In many applications, the efficiency should be calculated for an inhomogeneous
sample in the sense that it contains events with different weights. In order
to be able to determine the correct overall efficiency, it is necessary to
use for each subsample (= all events with the same weight) a different
TEfficiency object. After finishing your analysis you can then construct the
overall efficiency with its uncertainty.

This procedure has the advantage that you can change the weight of one
subsample easily without rerunning the whole analysis. On the other hand, more
effort is needed to handle several TEfficiency objects instead of one
histogram. In the case of many different or even continuously distributed
weights, this approach becomes cumbersome. One possibility to overcome this
problem is the usage of binned weights.

### Example
In particle physics weights arises from the fact that you want to
normalise your results to a certain reference value. A very common formula for
calculating weights is

\f{eqnarray*}{
   w &=& \frac{\sigma L}{N_{gen} \epsilon_{trig}} \\
     &-& \sigma ...\ cross\ section \\
     &-& L ...\ luminosity \\
     &-& N_{gen}\ ... number\ of\ generated\ events \\
     &-& \epsilon_{trig}\ ...\ (known)\ trigger\ efficiency \\
\f}

The reason for different weights can therefore be:
- different processes
- other integrated luminosity
- varying trigger efficiency
- different sample sizes
- ...
- or even combination of them

Depending on the actual meaning of different weights in your case, you
should either merge or combine them to get the overall efficiency.

### V.1 When should I use merging?
If the weights are artificial and do not represent real alternative hypotheses,
you should merge the different TEfficiency objects. That means especially for
the Bayesian case that the prior probability should be the same for all merged
TEfficiency objects. The merging can be done by invoking one of the following
operations:
- eff1.Add(eff2)
- eff1 += eff2
- eff1 = eff1 + eff2

The result of the merging is stored in the TEfficiency object which is marked
bold above. The contents of the internal histograms of both TEfficiency
objects are added and a new weight is assigned. The statistic options are not
changed.

\f[
   \frac{1}{w_{new}} = \frac{1}{w_{1}} + \frac{1}{w_{2}}
\f]

### Example:
If you use two samples with different numbers of generated events for the same
process and you want to normalise both to the same integrated luminosity and
trigger efficiency, the different weights then arise just from the fact that
you have different numbers of events. The TEfficiency objects should be merged
because the samples do not represent true alternatives. You expect the same
result as if you would have a big sample with all events in it.

\f[
   w_{1} = \frac{\sigma L}{\epsilon N_{1}}, w_{2} = \frac{\sigma L}{\epsilon N_{2}} \Rightarrow w_{new} = \frac{\sigma L}{\epsilon (N_{1} + N_{2})} = \frac{1}{\frac{1}{w_{1}} + \frac{1}{w_{2}}}
\f]

### V.2 When should I use combining?
You should combine TEfficiency objects whenever the weights represent
alternatives processes for the efficiency. As the combination of two TEfficiency
objects is not always consistent with the representation by two internal
histograms, the result is not stored in a TEfficiency object but a TGraphAsymmErrors
is returned which shows the estimated combined efficiency and its uncertainty
for each bin.
At the moment the combination method TEfficiency::Combine only supports a combination of 1-dimensional
efficiencies in a Bayesian approach.


For calculating the combined efficiency and its uncertainty for each bin only Bayesian statistics
is used. No frequentists methods are presently supported for computing the combined efficiency and
its confidence interval.
In the case of the Bayesian statistics, a combined posterior is constructed taking into account the
weight of each TEfficiency object. The same prior is used for all the TEfficiency objects.

\f{eqnarray*}{
  P_{comb}(\epsilon | {w_{i}}, {k_{i}} , {N_{i}}) = \frac{1}{norm} \prod_{i}{L(k_{i} | N_{i}, \epsilon)}^{w_{i}} \Pi( \epsilon )\\
L(k_{i} | N_{i}, \epsilon)\ is\ the\ likelihood\ function\ for\ the\ sample\ i\ (a\ Binomial\ distribution)\\
\Pi( \epsilon)\ is\ the\ prior,\ a\ beta\ distribution\ B(\epsilon, \alpha, \beta).\\
The\ resulting\ combined\ posterior\ is \\
P_{comb}(\epsilon |{w_{i}}; {k_{i}}; {N_{i}}) = B(\epsilon, \sum_{i}{ w_{i} k_{i}} + \alpha, \sum_{i}{ w_{i}(n_{i}-k_{i})}+\beta) \\
\hat{\varepsilon} = \int_{0}^{1} \epsilon \times P_{comb}(\epsilon | {k_{i}} , {N_{i}}) d\epsilon \\
confidence\ level = 1 - \alpha \\
\frac{\alpha}{2} = \int_{0}^{\epsilon_{low}} P_{comb}(\epsilon | {k_{i}} , {N_{i}}) d\epsilon ...\ defines\ lower\ boundary \\
1- \frac{\alpha}{2} = \int_{0}^{\epsilon_{up}} P_{comb}(\epsilon | {k_{i}} , {N_{i}}) d\epsilon ...\ defines\ upper\ boundary
\f}


###Example:
If you use cuts to select electrons which can originate from two different
processes, you can determine the selection efficiency for each process. The
overall selection efficiency is then the combined efficiency. The weights to be used in the
combination should be the probability that an
electron comes from the corresponding process.

\f[
p_{1} = \frac{\sigma_{1}}{\sigma_{1} + \sigma_{2}} = \frac{N_{1}w_{1}}{N_{1}w_{1} + N_{2}w_{2}}\\
p_{2} = \frac{\sigma_{2}}{\sigma_{1} + \sigma_{2}} = \frac{N_{2}w_{2}}{N_{1}w_{1} + N_{2}w_{2}}
\f]

## VI. Further operations

### VI.Information about the internal histograms
The methods TEfficiency::GetPassedHistogram and TEfficiency::GetTotalHistogram
return a constant pointer to the internal histograms. They can be used to
obtain information about the internal histograms (e.g., the binning, number of passed / total events in a bin, mean values...).
One can obtain a clone of the internal histograms by calling TEfficiency::GetCopyPassedHisto or TEfficiency::GetCopyTotalHisto.
The returned histograms are completely independent from the current
TEfficiency object. By default, they are not attached to a directory to
avoid the duplication of data and the user is responsible for deleting them.


~~~~~~~~~~~~~~~{.cpp}
//open a root file which contains a TEfficiency object
TFile* pFile = new TFile("myfile.root","update");

//get TEfficiency object with name "my_eff"
TEfficiency* pEff = (TEfficiency*)pFile->Get("my_eff");

//get clone of total histogram
TH1* clone = pEff->GetCopyTotalHisto();

//change clone...
//save changes of clone directly
clone->Write();
//or append it to the current directory and write the file
//clone->SetDirectory(gDirectory);
//pFile->Write();

//delete histogram object
delete clone;
clone = 0;
~~~~~~~~~~~~~~~

It is also possible to set the internal total or passed histogram by using the
methods TEfficiency::SetPassedHistogram or TEfficiency::SetTotalHistogram.

In order to ensure the validity of the TEfficiency object, the consistency of the
new histogram and the stored histogram is checked. It might be
impossible sometimes to change the histograms in a consistent way. Therefore one can force
the replacement by passing the "f" option. Then the user has to ensure that the
other internal histogram is replaced as well and that the TEfficiency object is
in a valid state.

### VI.2 Fitting
The efficiency can be fitted using the TEfficiency::Fit function which internally uses
the TBinomialEfficiencyFitter::Fit method.
As this method is using a maximum-likelihood-fit, it is necessary to initialise
the given fit function with reasonable start values.
The resulting fit function is attached to the list of associated functions and
will be drawn automatically during the next TEfficiency::Draw command.
The list of associated function can be modified by using the pointer returned
by TEfficiency::GetListOfFunctions.

Begin_Macro(source)
{
   //canvas only needed for this documentation
   TCanvas* c1 = new TCanvas("example","",600,400);
   c1->SetFillStyle(1001);
   c1->SetFillColor(kWhite);

   //create one-dimensional TEfficiency object with fixed bin size
   TEfficiency* pEff = new TEfficiency("eff","my efficiency;x;#epsilon",20,0,10);
   TRandom3 rand3;

   bool bPassed;
   double x;
   for(int i=0; i<10000; ++i)
   {
      //simulate events with variable under investigation
      x = rand3.Uniform(10);
      //check selection: bPassed = DoesEventPassSelection(x)
      bPassed = rand3.Rndm() < TMath::Gaus(x,5,4);
      pEff->Fill(bPassed,x);
   }

   //create a function for fitting and do the fit
   TF1* f1 = new TF1("f1","gaus",0,10);
   f1->SetParameters(1,5,2);
   pEff->Fit(f1);

   //create a threshold function
   TF1* f2 = new TF1("thres","0.8",0,10);
   f2->SetLineColor(kRed);
   //add it to the list of functions
   //use add first because the parameters of the last function will be displayed
   pEff->GetListOfFunctions()->AddFirst(f2);

   pEff->Draw("AP");

   //only for this documentation
   return c1;
}
End_Macro

### VI.3 Draw a TEfficiency object
A TEfficiency object can be drawn by calling the usual TEfficiency::Draw method.
At the moment drawing is only supported for 1- and 2-dimensional TEfficiency objects.
In the 1-dimensional case, you can use the same options as for the TGraphAsymmErrors::Draw
method. For 2-dimensional TEfficiency objects, you can pass the same options as
for a TH2::Draw object.

********************************************************************************/
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
///default constructor
///
///should not be used explicitly

TEfficiency::TEfficiency():
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fBoundary(0),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fPassedHistogram(0),
fTotalHistogram(0),
fWeight(kDefWeight)
{
   SetStatisticOption(kDefStatOpt);

   // create 2 dummy histograms
   fPassedHistogram = new TH1F("h_passed","passed",10,0,10);
   fTotalHistogram = new TH1F("h_total","total",10,0,10);
}

////////////////////////////////////////////////////////////////////////////////
///constructor using two existing histograms as input
///
///Input: passed - contains the events fulfilling some criteria
///       total  - contains all investigated events
///
///Notes: - both histograms have to fulfill the conditions of CheckConsistency
///       - dimension of the resulting efficiency object depends
///         on the dimension of the given histograms
///       - Clones of both histograms are stored internally
///       - The function SetName(total.GetName() + "_clone") is called to set
///         the names of the new object and the internal histograms..
///       - The created TEfficiency object is NOT appended to a directory. It
///         will not be written to disk during the next TFile::Write() command
///         in order to prevent duplication of data. If you want to save this
///         TEfficiency object anyway, you can either append it to a
///         directory by calling SetDirectory(TDirectory*) or write it
///         explicitly to disk by calling Write().

TEfficiency::TEfficiency(const TH1& passed,const TH1& total):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   //check consistency of histograms
   if(CheckConsistency(passed,total)) {
      // do not add cloned histograms to gDirectory
      TDirectory::TContext ctx(nullptr);
      fTotalHistogram = (TH1*)total.Clone();
      fPassedHistogram = (TH1*)passed.Clone();

      TString newName = total.GetName();
      newName += TString("_clone");
      SetName(newName);

      // are the histograms filled with weights?
      if(CheckWeights(passed,total))
      {
         Info("TEfficiency","given histograms are filled with weights");
         SetUseWeightedEvents();
      }
   }
   else {
      Error("TEfficiency(const TH1&,const TH1&)","histograms are not consistent -> results are useless");
      Warning("TEfficiency(const TH1&,const TH1&)","using two empty TH1D('h1','h1',10,0,10)");

      // do not add new created histograms to gDirectory
      TDirectory::TContext ctx(nullptr);
      fTotalHistogram = new TH1D("h1_total","h1 (total)",10,0,10);
      fPassedHistogram = new TH1D("h1_passed","h1 (passed)",10,0,10);
   }

   SetBit(kPosteriorMode,false);
   SetBit(kShortestInterval,false);

   SetStatisticOption(kDefStatOpt);
   SetDirectory(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 1-dimensional TEfficiency object with variable bin size.
///
/// Constructor creates two new and empty histograms with a given binning
///
/// Input:
///
///   - `name`: the common part of the name for both histograms (no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel".
///   - `nbins`: number of bins on the x-axis
///   - `xbins`: array of length (nbins + 1) with low-edges for each bin
///      xbins[nbinsx] ... lower edge for overflow bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbins,
                         const Double_t* xbins):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH1D("total","total",nbins,xbins);
   fPassedHistogram = new TH1D("passed","passed",nbins,xbins);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 1-dimensional TEfficiency object with fixed bins size.
///
/// Constructor creates two new and empty histograms with a fixed binning.
///
/// Input:
///
///   - `name`: the common part of the name for both histograms(no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel".
///   - `nbinsx`: number of bins on the x-axis
///   - `xlow`: lower edge of first bin
///   - `xup`: upper edge of last bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbinsx,
                         Double_t xlow,Double_t xup):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH1D("total","total",nbinsx,xlow,xup);
   fPassedHistogram = new TH1D("passed","passed",nbinsx,xlow,xup);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 2-dimensional TEfficiency object with fixed bin size.
///
/// Constructor creates two new and empty histograms with a fixed binning.
///
/// Input:
///
///   - `name`: the common part of the name for both histograms(no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel;zlabel".
///   - `nbinsx`: number of bins on the x-axis
///   - `xlow`: lower edge of first x-bin
///   - `xup`: upper edge of last x-bin
///   - `nbinsy`: number of bins on the y-axis
///   - `ylow`: lower edge of first y-bin
///   - `yup`: upper edge of last y-bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbinsx,
                         Double_t xlow,Double_t xup,Int_t nbinsy,
                         Double_t ylow,Double_t yup):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH2D("total","total",nbinsx,xlow,xup,nbinsy,ylow,yup);
   fPassedHistogram = new TH2D("passed","passed",nbinsx,xlow,xup,nbinsy,ylow,yup);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 2-dimensional TEfficiency object with variable bin size.
///
/// Constructor creates two new and empty histograms with a given binning.
///
/// Input:
///
///   - `name`: the common part of the name for both histograms(no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel;zlabel".
///   - `nbinsx`: number of bins on the x-axis
///   - `xbins`: array of length (nbins + 1) with low-edges for each bin
///      xbins[nbinsx] ... lower edge for overflow x-bin
///   - `nbinsy`: number of bins on the y-axis
///   - `ybins`: array of length (nbins + 1) with low-edges for each bin
///      ybins[nbinsy] ... lower edge for overflow y-bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbinsx,
                         const Double_t* xbins,Int_t nbinsy,
                         const Double_t* ybins):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH2D("total","total",nbinsx,xbins,nbinsy,ybins);
   fPassedHistogram = new TH2D("passed","passed",nbinsx,xbins,nbinsy,ybins);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 3-dimensional TEfficiency object with fixed bin size.
///
/// Constructor creates two new and empty histograms with a fixed binning.
///
/// Input:
///
///   - `name`: the common part of the name for both histograms(no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel;zlabel".
///   - `nbinsx`: number of bins on the x-axis
///   - `xlow`: lower edge of first x-bin
///   - `xup`: upper edge of last x-bin
///   - `nbinsy`: number of bins on the y-axis
///   - `ylow`: lower edge of first y-bin
///   - `yup`: upper edge of last y-bin
///   - `nbinsz`: number of bins on the z-axis
///   - `zlow`: lower edge of first z-bin
///   - `zup`: upper edge of last z-bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbinsx,
                         Double_t xlow,Double_t xup,Int_t nbinsy,
                         Double_t ylow,Double_t yup,Int_t nbinsz,
                         Double_t zlow,Double_t zup):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH3D("total","total",nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup);
   fPassedHistogram = new TH3D("passed","passed",nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Create 3-dimensional TEfficiency object with variable bin size.
///
/// Constructor creates two new and empty histograms with a given binning.
///
/// Input:
///
///   - `name`: the common part of the name for both histograms(no blanks)
///      fTotalHistogram has name: name + "_total"
///      fPassedHistogram has name: name + "_passed"
///   - `title`: the common part of the title for both histogram
///      fTotalHistogram has title: title + " (total)"
///      fPassedHistogram has title: title + " (passed)"
///      It is possible to label the axis by passing a title with
///      the following format: "title;xlabel;ylabel;zlabel".
///   - `nbinsx`: number of bins on the x-axis
///   - `xbins`: array of length (nbins + 1) with low-edges for each bin
///      xbins[nbinsx] ... lower edge for overflow x-bin
///   - `nbinsy`: number of bins on the y-axis
///   - `ybins`: array of length (nbins + 1) with low-edges for each bin
///      xbins[nbinsx] ... lower edge for overflow y-bin
///   - `nbinsz`: number of bins on the z-axis
///   - `zbins`: array of length (nbins + 1) with low-edges for each bin
///      xbins[nbinsx] ... lower edge for overflow z-bin

TEfficiency::TEfficiency(const char* name,const char* title,Int_t nbinsx,
                         const Double_t* xbins,Int_t nbinsy,
                         const Double_t* ybins,Int_t nbinsz,
                         const Double_t* zbins):
fBeta_alpha(kDefBetaAlpha),
fBeta_beta(kDefBetaBeta),
fConfLevel(kDefConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(kDefWeight)
{
   // do not add new created histograms to gDirectory
   TDirectory::TContext ctx(nullptr);
   fTotalHistogram = new TH3D("total","total",nbinsx,xbins,nbinsy,ybins,nbinsz,zbins);
   fPassedHistogram = new TH3D("passed","passed",nbinsx,xbins,nbinsy,ybins,nbinsz,zbins);

   Build(name,title);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
///
///The list of associated objects (e.g. fitted functions) is not copied.
///
///Note:
///
///   - SetName(rEff.GetName() + "_copy") is called to set the names of the
///        object and the histograms.
///   - The titles are set by calling SetTitle("[copy] " + rEff.GetTitle()).
///   - The copied TEfficiency object is NOT appended to a directory. It
///      will not be written to disk during the next TFile::Write() command
///      in order to prevent duplication of data. If you want to save this
///      TEfficiency object anyway, you can either append it to a directory
///      by calling SetDirectory(TDirectory*) or write it explicitly to disk
///      by calling Write().

TEfficiency::TEfficiency(const TEfficiency& rEff):
TNamed(),
TAttLine(),
TAttFill(),
TAttMarker(),
fBeta_alpha(rEff.fBeta_alpha),
fBeta_beta(rEff.fBeta_beta),
fBeta_bin_params(rEff.fBeta_bin_params),
fConfLevel(rEff.fConfLevel),
fDirectory(0),
fFunctions(0),
fPaintGraph(0),
fPaintHisto(0),
fWeight(rEff.fWeight)
{
   // copy TObject bits
   ((TObject&)rEff).Copy(*this);

   // do not add cloned histograms to gDirectory
   {
      TDirectory::TContext ctx(nullptr);
      fTotalHistogram = (TH1*)((rEff.fTotalHistogram)->Clone());
      fPassedHistogram = (TH1*)((rEff.fPassedHistogram)->Clone());
   }

   TString name = rEff.GetName();
   name += "_copy";
   SetName(name);
   TString title = "[copy] ";
   title += rEff.GetTitle();
   SetTitle(title);

   SetStatisticOption(rEff.GetStatisticOption());

   SetDirectory(0);

   //copy style
   rEff.TAttLine::Copy(*this);
   rEff.TAttFill::Copy(*this);
   rEff.TAttMarker::Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
///default destructor

TEfficiency::~TEfficiency()
{
   //delete all function in fFunctions
   // use same logic as in TH1 destructor
   // (see TH1::~TH1 code in TH1.cxx)
   if(fFunctions) {
      fFunctions->SetBit(kInvalidObject);
      TObject* obj = 0;
      while ((obj  = fFunctions->First())) {
         while(fFunctions->Remove(obj)) { }
         if (!obj->TestBit(kNotDeleted)) {
            break;
         }
         delete obj;
         obj = 0;
      }
      delete fFunctions;
      fFunctions = 0;
   }

   if(fDirectory)
      fDirectory->Remove(this);

   delete fTotalHistogram;
   delete fPassedHistogram;
   delete fPaintGraph;
   delete fPaintHisto;
}

////////////////////////////////////////////////////////////////////////////////
/**
    Calculates the boundaries for the frequentist Agresti-Coull interval

    \param total number of total events
    \param passed 0 <= number of passed events <= total
    \param level  confidence level
    \param bUpper true  - upper boundary is returned
    false - lower boundary is returned


   \f{eqnarray*}{
     \alpha &=& 1 - \frac{level}{2} \\
     \kappa &=& \Phi^{-1}(1 - \alpha,1)\ ... normal\ quantile\ function\\
     mode &=& \frac{passed + \frac{\kappa^{2}}{2}}{total + \kappa^{2}}\\
     \Delta &=& \kappa * \sqrt{\frac{mode * (1 - mode)}{total + \kappa^{2}}}\\
     return &=&  max(0,mode - \Delta)\ or\ min(1,mode + \Delta)
   \f}

*/

Double_t TEfficiency::AgrestiCoull(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   Double_t alpha = (1.0 - level)/2;
   Double_t kappa = ROOT::Math::normal_quantile(1 - alpha,1);

   Double_t mode = (passed + 0.5 * kappa * kappa) / (total + kappa * kappa);
   Double_t delta = kappa * std::sqrt(mode * (1 - mode) / (total + kappa * kappa));

   if(bUpper)
      return ((mode + delta) > 1) ? 1.0 : (mode + delta);
   else
      return ((mode - delta) < 0) ? 0.0 : (mode - delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the boundaries for the frequentist Feldman-Cousins interval
///
/// \param total number of total events
/// \param passed 0 <= number of passed events <= total
/// \param level confidence level
/// \param bUpper: true  - upper boundary is returned
///                false - lower boundary is returned

Double_t TEfficiency::FeldmanCousins(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   Double_t lower = 0;
   Double_t upper = 1;
   if (!FeldmanCousinsInterval(total,passed,level, lower, upper)) {
      ::Error("FeldmanCousins","Error running FC method - return 0 or 1");
   }
   return (bUpper) ? upper : lower;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the interval boundaries using the frequentist methods of Feldman-Cousins
///
/// \param[in] total number of total events
/// \param[in] passed 0 <= number of passed events <= total
/// \param[in] level  confidence level
/// \param[out] lower lower boundary returned on exit
/// \param[out] upper lower boundary returned on exit
/// \return a flag with the status of the calculation
///
/// Calculation:
///
/// The Feldman-Cousins is a frequentist method where the interval is estimated using a Neyman construction where the ordering
/// is based on the likelihood ratio:
/// \f[
///   LR =  \frac{Binomial(k | N, \epsilon)}{Binomial(k | N, \hat{\epsilon} ) }
/// \f]
/// See G. J. Feldman and R. D. Cousins, Phys. Rev. D57 (1998) 3873
/// and   R. D. Cousins, K. E. Hymes, J. Tucker, Nuclear Instruments and Methods in Physics Research A 612 (2010) 388
///
/// Implemented using classes developed by Jordan Tucker and Luca Lista
/// See File hist/hist/src/TEfficiencyHelper.h

Bool_t TEfficiency::FeldmanCousinsInterval(Double_t total,Double_t passed,Double_t level,Double_t & lower, Double_t & upper)
{
   FeldmanCousinsBinomialInterval fc;
   double alpha = 1.-level;
   fc.Init(alpha);
   fc.Calculate(passed, total);
   lower = fc.Lower();
   upper = fc.Upper();
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the boundaries using the  mid-P binomial
/// interval (Lancaster method)  from B. Cousing and J. Tucker.
/// See http://arxiv.org/abs/0905.3831 for a description and references for the method
///
/// Modify equal_tailed to get the kind of interval you want.
/// Can also be converted to interval on ratio of poisson means X/Y by the substitutions
/// ~~~ {.cpp}
///  X = passed
///  total = X + Y
///  lower_poisson = lower/(1 - lower)
///  upper_poisson = upper/(1 - upper)
/// ~~~

Double_t TEfficiency::MidPInterval(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   const double alpha = 1. - level;
   const bool equal_tailed = true;  // change if you don;t want equal tailed interval
   const double alpha_min = equal_tailed ? alpha/2 : alpha;
   const double tol = 1e-9; // tolerance
   double pmin = 0;
   double pmax = 0;
   double p = 0;

   pmin = 0; pmax = 1;


   // treat special case for 0<passed<1
   // do a linear interpolation of the upper limit values
   if ( passed > 0 && passed < 1) {
      double p0 =  MidPInterval(total,0.0,level,bUpper);
      double p1 =  MidPInterval(total,1.0,level,bUpper);
      p = (p1 - p0) * passed + p0;
      return p;
   }

   while (std::abs(pmax - pmin) > tol) {
      p = (pmin + pmax)/2;
      //double v = 0.5 * ROOT::Math::binomial_pdf(int(passed), p, int(total));
      // make it work for non integer using the binomial - beta relationship
      double v = 0.5 * ROOT::Math::beta_pdf(p, passed+1., total-passed+1)/(total+1);
      //if (passed > 0) v += ROOT::Math::binomial_cdf(int(passed - 1), p, int(total));
      // compute the binomial cdf at passed -1
      if ( (passed-1) >= 0) v += ROOT::Math::beta_cdf_c(p, passed, total-passed+1);

      double vmin =  (bUpper) ? alpha_min : 1.- alpha_min;
      if (v > vmin)
         pmin = p;
      else
         pmax = p;
   }

   return p;
}


////////////////////////////////////////////////////////////////////////////////
/**
Calculates the boundaries for a Bayesian confidence interval (shortest or central interval depending on the option)

\param[in] total number of total events
\param[in] passed 0 <= number of passed events <= total
\param[in] level  confidence level
\param[in] alpha  shape parameter > 0 for the prior distribution (fBeta_alpha)
\param[in] beta  shape parameter > 0 for the prior distribution (fBeta_beta)
\param[in] bUpper
           - true  - upper boundary is returned
           - false - lower boundary is returned
\param[in] bShortest ??

Note: In the case central confidence interval is calculated.
     when passed = 0 (or passed = total) the lower (or upper)
     interval values will be larger than 0 (or smaller than 1).

Calculation:

The posterior probability in bayesian statistics is given by:
\f[
   P(\varepsilon |k,N) \propto L(\varepsilon|k,N) \times Prior(\varepsilon)
\f]
As an efficiency can be interpreted as probability of a positive outcome of
a Bernoullli trial the likelihood function is given by the binomial
distribution:
\f[
  L(\varepsilon|k,N) = Binomial(N,k) \varepsilon ^{k} (1 - \varepsilon)^{N-k}
\f]
At the moment only beta distributions are supported as prior probabilities
of the efficiency (\f$ B(\alpha,\beta)\f$ is the beta function):
\f[
  Prior(\varepsilon) = \frac{1}{B(\alpha,\beta)} \varepsilon ^{\alpha - 1} (1 - \varepsilon)^{\beta - 1}
\f]
The posterior probability is therefore again given by a beta distribution:
\f[
  P(\varepsilon |k,N) \propto \varepsilon ^{k + \alpha - 1} (1 - \varepsilon)^{N - k + \beta - 1}
\f]
In case of central intervals
the lower boundary for the equal-tailed confidence interval is given by the
inverse cumulative (= quantile) function for the quantile \f$ \frac{1 - level}{2} \f$.
The upper boundary for the equal-tailed confidence interval is given by the
inverse cumulative (= quantile) function for the quantile \f$ \frac{1 + level}{2} \f$.
Hence it is the solution \f$ \varepsilon \f$ of the following equation:
\f[
  I_{\varepsilon}(k + \alpha,N - k + \beta) = \frac{1}{norm} \int_{0}^{\varepsilon} dt t^{k + \alpha - 1} (1 - t)^{N - k + \beta - 1} =  \frac{1 \pm level}{2}
\f]
In the case of shortest interval the minimum interval around the mode is found by minimizing the length of all intervals width the
given probability content. See TEfficiency::BetaShortestInterval
*/

Double_t TEfficiency::Bayesian(Double_t total,Double_t passed,Double_t level,Double_t alpha,Double_t beta,Bool_t bUpper, Bool_t bShortest)
{
   Double_t a = double(passed)+alpha;
   Double_t b = double(total-passed)+beta;

   if (bShortest) {
      double lower = 0;
      double upper = 1;
      BetaShortestInterval(level,a,b,lower,upper);
      return (bUpper) ? upper : lower;
   }
   else
      return BetaCentralInterval(level, a, b, bUpper);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the boundaries for a central confidence interval for a Beta distribution
///
/// \param[in] level  confidence level
/// \param[in] a  parameter > 0 for the beta distribution (for a posterior is passed + prior_alpha
/// \param[in] b  parameter > 0 for the beta distribution (for a posterior is (total-passed) + prior_beta
/// \param[in] bUpper true  - upper boundary is returned
///                   false - lower boundary is returned

Double_t TEfficiency::BetaCentralInterval(Double_t level,Double_t a,Double_t b,Bool_t bUpper)
{
   if(bUpper) {
      if((a > 0) && (b > 0))
         return ROOT::Math::beta_quantile((1+level)/2,a,b);
      else {
         gROOT->Error("TEfficiency::BayesianCentral","Invalid input parameters - return 1");
         return 1;
      }
   }
   else {
      if((a > 0) && (b > 0))
         return ROOT::Math::beta_quantile((1-level)/2,a,b);
      else {
         gROOT->Error("TEfficiency::BayesianCentral","Invalid input parameters - return 0");
         return 0;
      }
   }
}

struct Beta_interval_length {
   Beta_interval_length(Double_t level,Double_t alpha,Double_t beta ) :
   fCL(level), fAlpha(alpha), fBeta(beta)
   {}

   Double_t LowerMax() {
      // max allowed value of lower given the interval size
      return ROOT::Math::beta_quantile_c(fCL, fAlpha,fBeta);
   }

   Double_t operator() (double lower) const {
      // return length of interval
      Double_t plow = ROOT::Math::beta_cdf(lower, fAlpha, fBeta);
      Double_t pup = plow + fCL;
      double upper = ROOT::Math::beta_quantile(pup, fAlpha,fBeta);
      return upper-lower;
   }
   Double_t fCL; // interval size (confidence level)
   Double_t fAlpha; // beta distribution alpha parameter
   Double_t fBeta; // beta distribution beta parameter

};

////////////////////////////////////////////////////////////////////////////////
/// Calculates the boundaries for a shortest confidence interval for a Beta  distribution
///
/// \param[in] level  confidence level
/// \param[in] a  parameter > 0 for the beta distribution (for a posterior is passed + prior_alpha
/// \param[in] b  parameter > 0 for the beta distribution (for a posterior is (total-passed) + prior_beta
/// \param[out] upper upper boundary is returned
/// \param[out] lower lower boundary is returned
///
/// The lower/upper boundary are then obtained by finding the shortest interval of the beta distribution
/// contained the desired probability level.
/// The length of all possible intervals is minimized in order to find the shortest one

Bool_t TEfficiency::BetaShortestInterval(Double_t level,Double_t a,Double_t b, Double_t & lower, Double_t & upper)
{
   if (a <= 0 || b <= 0) {
      lower = 0; upper = 1;
      gROOT->Error("TEfficiency::BayesianShortest","Invalid input parameters - return [0,1]");
      return kFALSE;
   }

   // treat here special cases when mode == 0 or 1
   double mode = BetaMode(a,b);
   if (mode == 0.0) {
      lower = 0;
      upper = ROOT::Math::beta_quantile(level, a, b);
      return kTRUE;
   }
   if (mode == 1.0) {
      lower = ROOT::Math::beta_quantile_c(level, a, b);
      upper = 1.0;
      return kTRUE;
   }
   // special case when the shortest interval is undefined  return the central interval
   // can happen for a posterior when passed=total=0
   //
   if ( a==b && a<=1.0) {
      lower = BetaCentralInterval(level,a,b,kFALSE);
      upper = BetaCentralInterval(level,a,b,kTRUE);
      return kTRUE;
   }

   // for the other case perform a minimization
   // make a function of the length of the posterior interval as a function of lower bound
   Beta_interval_length intervalLength(level,a,b);
   // minimize the interval length
   ROOT::Math::WrappedFunction<const Beta_interval_length &> func(intervalLength);
   ROOT::Math::BrentMinimizer1D minim;
   minim.SetFunction(func, 0, intervalLength.LowerMax() );
   minim.SetNpx(2); // no need to bracket with many iterations. Just do few times to estimate some better points
   bool ret = minim.Minimize(100, 1.E-10,1.E-10);
   if (!ret) {
      gROOT->Error("TEfficiency::BayesianShortes","Error finding the shortest interval");
      return kFALSE;
   }
   lower = minim.XMinimum();
   upper = lower + minim.FValMinimum();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the mean (average) of the beta distribution
///
/// \param[in] a  parameter > 0 for the beta distribution (for a posterior is passed + prior_alpha
/// \param[in] b  parameter > 0 for the beta distribution (for a posterior is (total-passed) + prior_beta
///

Double_t TEfficiency::BetaMean(Double_t a,Double_t b)
{
   if (a <= 0 || b <= 0 ) {
      gROOT->Error("TEfficiency::BayesianMean","Invalid input parameters - return 0");
      return 0;
   }

   Double_t mean =  a / (a + b);
   return mean;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the mode of the beta distribution
///
/// \param[in] a  parameter > 0 for the beta distribution (for a posterior is passed + prior_alpha
/// \param[in] b  parameter > 0 for the beta distribution (for a posterior is (total-passed) + prior_beta
///
/// note the mode is defined for a Beta(a,b) only if (a,b)>1 (a = passed+alpha; b = total-passed+beta)
/// return then the following in case (a,b) < 1:
/// - if (a==b) return 0.5 (it is really undefined)
/// - if (a < b) return 0;
/// - if (a > b) return 1;

Double_t TEfficiency::BetaMode(Double_t a,Double_t b)
{
   if (a <= 0 || b <= 0 ) {
      gROOT->Error("TEfficiency::BayesianMode","Invalid input parameters - return 0");
      return 0;
   }
   if ( a <= 1 || b <= 1) {
      if ( a < b) return 0;
      if ( a > b) return 1;
      if (a == b) return 0.5; // cannot do otherwise
   }

   // since a and b are > 1 here denominator cannot be 0 or < 0
   Double_t mode =  (a - 1.0) / (a + b -2.0);
   return mode;
}
////////////////////////////////////////////////////////////////////////////////
/// Building standard data structure of a TEfficiency object
///
/// Notes:
/// - calls: SetName(name), SetTitle(title)
/// - set the statistic option to the default (kFCP)
/// - appends this object to the current directory SetDirectory(gDirectory)

void TEfficiency::Build(const char* name,const char* title)
{
   SetName(name);
   SetTitle(title);

   SetStatisticOption(kDefStatOpt);
   SetDirectory(gDirectory);

   SetBit(kPosteriorMode,false);
   SetBit(kShortestInterval,false);
   SetBit(kUseWeights,false);

   //set normalisation factors to 0, otherwise the += may not work properly
   fPassedHistogram->SetNormFactor(0);
   fTotalHistogram->SetNormFactor(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks binning for each axis
///
/// It is assumed that the passed histograms have the same dimension.

Bool_t TEfficiency::CheckBinning(const TH1& pass,const TH1& total)
{

   const TAxis* ax1 = 0;
   const TAxis* ax2 = 0;

   //check binning along axis
   for(Int_t j = 0; j < pass.GetDimension(); ++j) {
      switch(j) {
         case 0:
            ax1 = pass.GetXaxis();
            ax2 = total.GetXaxis();
            break;
         case 1:
            ax1 = pass.GetYaxis();
            ax2 = total.GetYaxis();
            break;
         case 2:
            ax1 = pass.GetZaxis();
            ax2 = total.GetZaxis();
            break;
      }

      if(ax1->GetNbins() != ax2->GetNbins()) {
         gROOT->Info("TEfficiency::CheckBinning","Histograms are not consistent: they have different number of bins");
         return false;
      }

      for(Int_t i = 1; i <= ax1->GetNbins() + 1; ++i)
         if(!TMath::AreEqualRel(ax1->GetBinLowEdge(i), ax2->GetBinLowEdge(i), 1.E-15)) {
            gROOT->Info("TEfficiency::CheckBinning","Histograms are not consistent: they have different bin edges");
            return false;
         }


   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks the consistence of the given histograms
///
/// The histograms are considered as consistent if:
/// - both have the same dimension
/// - both have the same binning
/// - pass.GetBinContent(i) <= total.GetBinContent(i) for each bin i
///

Bool_t TEfficiency::CheckConsistency(const TH1& pass,const TH1& total, Option_t*)
{
   if(pass.GetDimension() != total.GetDimension()) {
      gROOT->Error("TEfficiency::CheckConsistency","passed TEfficiency objects have different dimensions");
      return false;
   }

   if(!CheckBinning(pass,total)) {
      gROOT->Error("TEfficiency::CheckConsistency","passed TEfficiency objects have different binning");
      return false;
   }

   if(!CheckEntries(pass,total)) {
      gROOT->Error("TEfficiency::CheckConsistency","passed TEfficiency objects do not have consistent bin contents");
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks whether bin contents are compatible with binomial statistics
///
/// The following inequality has to be valid for each bin i:
/// total.GetBinContent(i) >= pass.GetBinContent(i)
///
///
///
/// Note:
///
///   - It is assumed that both histograms have the same dimension and binning.

Bool_t TEfficiency::CheckEntries(const TH1& pass,const TH1& total, Option_t*)
{

   //check: pass <= total
   Int_t nbinsx, nbinsy, nbinsz, nbins;

   nbinsx = pass.GetNbinsX();
   nbinsy = pass.GetNbinsY();
   nbinsz = pass.GetNbinsZ();

   switch(pass.GetDimension()) {
      case 1: nbins = nbinsx + 2; break;
      case 2: nbins = (nbinsx + 2) * (nbinsy + 2); break;
      case 3: nbins = (nbinsx + 2) * (nbinsy + 2) * (nbinsz + 2); break;
      default: nbins = 0;
   }

   for(Int_t i = 0; i < nbins; ++i) {
      if(pass.GetBinContent(i) > total.GetBinContent(i)) {
         gROOT->Info("TEfficiency::CheckEntries","Histograms are not consistent: passed bin content > total bin content");
         return false;
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if both histogram are weighted. If they are weighted a true is returned
///
Bool_t TEfficiency::CheckWeights(const TH1& pass,const TH1& total)
{
   if (pass.GetSumw2N() == 0 && total.GetSumw2N() == 0) return false;

   // check also that the total sum of weight and weight squares are consistent
   Double_t statpass[TH1::kNstat];
   Double_t stattotal[TH1::kNstat];

   pass.GetStats(statpass);
   total.GetStats(stattotal);

   double tolerance = (total.IsA() == TH1F::Class() ) ? 1.E-5 : 1.E-12;

   //require: sum of weights == sum of weights^2
   if(!TMath::AreEqualRel(statpass[0],statpass[1],tolerance) ||
      !TMath::AreEqualRel(stattotal[0],stattotal[1],tolerance) ) {
      return true;
   }

   // histograms are not weighted
   return false;

}


////////////////////////////////////////////////////////////////////////////////
/// Create the graph used be painted (for dim=1 TEfficiency)
/// The return object is managed by the caller

TGraphAsymmErrors * TEfficiency::CreateGraph(Option_t * opt) const
{
   if (GetDimension() != 1) {
      Error("CreatePaintingGraph","Call this function only for dimension == 1");
      return 0;
   }


   Int_t npoints = fTotalHistogram->GetNbinsX();
   TGraphAsymmErrors * graph = new TGraphAsymmErrors(npoints);
   graph->SetName("eff_graph");
   FillGraph(graph,opt);

   return graph;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the graph to be painted with information from TEfficiency
/// Internal method called by TEfficiency::Paint or TEfficiency::CreateGraph

void TEfficiency::FillGraph(TGraphAsymmErrors * graph, Option_t * opt) const
{
   TString option = opt;
   option.ToLower();

   Bool_t plot0Bins = false;
   if (option.Contains("e0") ) plot0Bins = true;

   Double_t x,y,xlow,xup,ylow,yup;
   //point i corresponds to bin i+1 in histogram
   // point j is point graph index
   // LM: cannot use TGraph::SetPoint because it deletes the underlying
   // histogram  each time (see TGraph::SetPoint)
   // so use it only when extra points are added to the graph
   Int_t j = 0;
   double * px = graph->GetX();
   double * py = graph->GetY();
   double * exl = graph->GetEXlow();
   double * exh = graph->GetEXhigh();
   double * eyl = graph->GetEYlow();
   double * eyh = graph->GetEYhigh();
   Int_t npoints = fTotalHistogram->GetNbinsX();
   for (Int_t i = 0; i < npoints; ++i) {
      if (!plot0Bins && fTotalHistogram->GetBinContent(i+1) == 0 )    continue;
      x = fTotalHistogram->GetBinCenter(i+1);
      y = GetEfficiency(i+1);
      xlow = fTotalHistogram->GetBinCenter(i+1) - fTotalHistogram->GetBinLowEdge(i+1);
      xup = fTotalHistogram->GetBinWidth(i+1) - xlow;
      ylow = GetEfficiencyErrorLow(i+1);
      yup = GetEfficiencyErrorUp(i+1);
      // in the case the graph already existed and extra points have been added
      if (j >= graph->GetN() ) {
         graph->SetPoint(j,x,y);
         graph->SetPointError(j,xlow,xup,ylow,yup);
      }
      else {
         px[j] = x;
         py[j] = y;
         exl[j] = xlow;
         exh[j] = xup;
         eyl[j] = ylow;
         eyh[j] = yup;
      }
      j++;
   }

   // tell the graph the effective number of points
   graph->Set(j);
   //refresh title before painting if changed
   TString oldTitle = graph->GetTitle();
   TString newTitle = GetTitle();
   if (oldTitle != newTitle ) {
      graph->SetTitle(newTitle);
   }

   // set the axis labels
   TString xlabel = fTotalHistogram->GetXaxis()->GetTitle();
   TString ylabel = fTotalHistogram->GetYaxis()->GetTitle();
   if (xlabel) graph->GetXaxis()->SetTitle(xlabel);
   if (ylabel) graph->GetYaxis()->SetTitle(ylabel);

   //copying style information
   TAttLine::Copy(*graph);
   TAttFill::Copy(*graph);
   TAttMarker::Copy(*graph);

   // copy axis labels if existing. Assume are there in the total histogram
   if (fTotalHistogram->GetXaxis()->GetLabels() != nullptr) {
      for (int ibin = 1; ibin <= fTotalHistogram->GetXaxis()->GetNbins(); ++ibin) {
         // we need to find the right bin for the Histogram representing the xaxis of the graph
         int grbin = graph->GetXaxis()->FindBin(fTotalHistogram->GetXaxis()->GetBinCenter(ibin));
         graph->GetXaxis()->SetBinLabel(grbin, fTotalHistogram->GetXaxis()->GetBinLabel(ibin));
      }
   }
   // this method forces the graph to compute correctly the axis
   // according to the given points
   graph->GetHistogram();

}

////////////////////////////////////////////////////////////////////////////////
/// Create the histogram used to be painted (for dim=2 TEfficiency)
/// The return object is managed by the caller

TH2 * TEfficiency::CreateHistogram(Option_t *) const
{
   if (GetDimension() != 2) {
      Error("CreatePaintingistogram","Call this function only for dimension == 2");
      return 0;
   }

   Int_t nbinsx = fTotalHistogram->GetNbinsX();
   Int_t nbinsy = fTotalHistogram->GetNbinsY();
   TAxis * xaxis = fTotalHistogram->GetXaxis();
   TAxis * yaxis = fTotalHistogram->GetYaxis();
   TH2 * hist = 0;

   if (xaxis->IsVariableBinSize() && yaxis->IsVariableBinSize() )
      hist = new TH2F("eff_histo",GetTitle(),nbinsx,xaxis->GetXbins()->GetArray(),
                      nbinsy,yaxis->GetXbins()->GetArray());
   else if (xaxis->IsVariableBinSize() && ! yaxis->IsVariableBinSize() )
      hist = new TH2F("eff_histo",GetTitle(),nbinsx,xaxis->GetXbins()->GetArray(),
                      nbinsy,yaxis->GetXmin(), yaxis->GetXmax());
   else if (!xaxis->IsVariableBinSize() &&  yaxis->IsVariableBinSize() )
      hist = new TH2F("eff_histo",GetTitle(),nbinsx,xaxis->GetXmin(), xaxis->GetXmax(),
                      nbinsy,yaxis->GetXbins()->GetArray());
   else
      hist = new TH2F("eff_histo",GetTitle(),nbinsx,xaxis->GetXmin(), xaxis->GetXmax(),
                      nbinsy,yaxis->GetXmin(), yaxis->GetXmax());


   hist->SetDirectory(0);

   FillHistogram(hist);

   return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the 2d histogram to be painted with information from TEfficiency 2D
/// Internal method called by TEfficiency::Paint or TEfficiency::CreatePaintingGraph

void TEfficiency::FillHistogram(TH2 * hist ) const
{
   //refresh title before each painting
   hist->SetTitle(GetTitle());

   // set the axis labels
   TString xlabel = fTotalHistogram->GetXaxis()->GetTitle();
   TString ylabel = fTotalHistogram->GetYaxis()->GetTitle();
   TString zlabel = fTotalHistogram->GetZaxis()->GetTitle();
   if (xlabel) hist->GetXaxis()->SetTitle(xlabel);
   if (ylabel) hist->GetYaxis()->SetTitle(ylabel);
   if (zlabel) hist->GetZaxis()->SetTitle(zlabel);

   Int_t bin;
   Int_t nbinsx = hist->GetNbinsX();
   Int_t nbinsy = hist->GetNbinsY();
   for(Int_t i = 0; i < nbinsx + 2; ++i) {
      for(Int_t j = 0; j < nbinsy + 2; ++j) {
         bin = GetGlobalBin(i,j);
         hist->SetBinContent(bin,GetEfficiency(bin));
      }
   }

   // copy axis labels if existing. Assume are there in the total histogram
   if (fTotalHistogram->GetXaxis()->GetLabels() != nullptr) {
      for (int ibinx = 1; ibinx <= fTotalHistogram->GetXaxis()->GetNbins(); ++ibinx)
         hist->GetXaxis()->SetBinLabel(ibinx, fTotalHistogram->GetXaxis()->GetBinLabel(ibinx));
   }
   if (fTotalHistogram->GetYaxis()->GetLabels() != nullptr) {
      for (int ibiny = 1; ibiny <= fTotalHistogram->GetYaxis()->GetNbins(); ++ibiny)
         hist->GetYaxis()->SetBinLabel(ibiny, fTotalHistogram->GetYaxis()->GetBinLabel(ibiny));
   }

   //copying style information
   TAttLine::Copy(*hist);
   TAttFill::Copy(*hist);
   TAttMarker::Copy(*hist);
   hist->SetStats(0);

   return;

}
////////////////////////////////////////////////////////////////////////////////
/**
Calculates the boundaries for the frequentist Clopper-Pearson interval

This interval is recommended by the PDG.

\param[in] total number of total events
\param[in] passed 0 <= number of passed events <= total
\param[in] level confidence level
\param[in] bUpper true  - upper boundary is returned
                 ;false - lower boundary is returned

Calculation:

The lower boundary of the Clopper-Pearson interval is the "exact" inversion
of the test:
   \f{eqnarray*}{
      P(x \geq passed; total) &=& \frac{1 - level}{2}\\
      P(x \geq passed; total) &=& 1 - P(x \leq passed - 1; total)\\
      &=& 1 - \frac{1}{norm} * \int_{0}^{1 - \varepsilon} t^{total - passed} (1 - t)^{passed - 1} dt\\
      &=& 1 - \frac{1}{norm} * \int_{\varepsilon}^{1} t^{passed - 1} (1 - t)^{total - passed} dt\\
      &=& \frac{1}{norm} * \int_{0}^{\varepsilon} t^{passed - 1} (1 - t)^{total - passed} dt\\
      &=& I_{\varepsilon}(passed,total - passed + 1)
    \f}
The lower boundary is therefore given by the \f$ \frac{1 - level}{2}\f$ quantile
of the beta distribution.

The upper boundary of the Clopper-Pearson interval is the "exact" inversion
of the test:
   \f{eqnarray*}{
      P(x \leq passed; total) &=& \frac{1 - level}{2}\\
      P(x \leq passed; total) &=& \frac{1}{norm} * \int_{0}^{1 - \varepsilon} t^{total - passed - 1} (1 - t)^{passed} dt\\
      &=& \frac{1}{norm} * \int_{\varepsilon}^{1} t^{passed} (1 - t)^{total - passed - 1} dt\\
      &=& 1 - \frac{1}{norm} * \int_{0}^{\varepsilon} t^{passed} (1 - t)^{total - passed - 1} dt\\
      \Rightarrow 1 - \frac{1 - level}{2} &=& \frac{1}{norm} * \int_{0}^{\varepsilon} t^{passed} (1 - t)^{total - passed -1} dt\\
      \frac{1 + level}{2} &=& I_{\varepsilon}(passed + 1,total - passed)
    \f}
The upper boundary is therefore given by the \f$\frac{1 + level}{2}\f$ quantile
of the beta distribution.

Note: The connection between the binomial distribution and the regularized
     incomplete beta function \f$ I_{\varepsilon}(\alpha,\beta)\f$ has been used.
*/

Double_t TEfficiency::ClopperPearson(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   Double_t alpha = (1.0 - level) / 2;
   if(bUpper)
      return ((passed == total) ? 1.0 : ROOT::Math::beta_quantile(1 - alpha,passed + 1,total-passed));
   else
      return ((passed == 0) ? 0.0 : ROOT::Math::beta_quantile(alpha,passed,total-passed+1.0));
}
////////////////////////////////////////////////////////////////////////////////
/**
    Calculates the combined efficiency and its uncertainties

    This method does a bayesian combination of the given samples.

    \param[in] up  contains the upper limit of the confidence interval afterwards
    \param[in] low  contains the lower limit of the confidence interval afterwards
    \param[in] n    number of samples which are combined
    \param[in] pass array of length n containing the number of passed events
    \param[in] total array of length n containing the corresponding numbers of total events
    \param[in] alpha  shape parameters for the beta distribution as prior
    \param[in] beta   shape parameters for the beta distribution as prior
    \param[in] level  desired confidence level
    \param[in] w weights for each sample; if not given, all samples get the weight 1
              The weights do not need to be normalized, since they are internally renormalized
              to the number of effective entries.
    \param[in] opt
      -  mode : The mode is returned instead of the mean of the posterior as best value
                When using the mode the shortest interval is also computed instead of the central one
      -  shortest: compute shortest interval (done by default if mode option is set)
      -  central: compute central interval (done by default if mode option is NOT set)

    Calculation:

    The combined posterior distributions is calculated from the Bayes theorem assuming a common prior Beta distribution.
        It is easy to proof that the combined posterior is then:
 \f{eqnarray*}{
      P_{comb}(\epsilon |{w_{i}}; {k_{i}}; {N_{i}}) &=& B(\epsilon, \sum_{i}{ w_{i} k_{i}} + \alpha, \sum_{i}{ w_{i}(n_{i}-k_{i})}+\beta)\\
      w_{i} &=& weight\ for\ each\ sample\ renormalized\ to\ the\ effective\ entries\\
      w^{'}_{i} &=&  w_{i} \frac{ \sum_{i} {w_{i} } } { \sum_{i} {w_{i}^{2} } }
    \f}

    The estimated efficiency is the mode (or the mean) of the obtained posterior distribution

    The boundaries of the confidence interval for a confidence level (1 - a)
    are given by the a/2 and 1-a/2 quantiles of the resulting cumulative
    distribution.

    Example (uniform prior distribution):

Begin_Macro(source)
{
     TCanvas* c1 = new TCanvas("c1","",600,800);
     c1->Divide(1,2);
     c1->SetFillStyle(1001);
     c1->SetFillColor(kWhite);

     TF1* p1 = new TF1("p1","TMath::BetaDist(x,19,9)",0,1);
     TF1* p2 = new TF1("p2","TMath::BetaDist(x,4,8)",0,1);
     TF1* comb = new TF1("comb2","TMath::BetaDist(x,[0],[1])",0,1);
     double nrm = 1./(0.6*0.6+0.4*0.4); // weight normalization
     double a = 0.6*18.0 + 0.4*3.0 + 1.0;  // new alpha parameter of combined beta dist.
     double b = 0.6*10+0.4*7+1.0;  // new beta parameter of combined beta dist.
     comb->SetParameters(nrm*a ,nrm *b );
     TF1* const1 = new TF1("const1","0.05",0,1);
     TF1* const2 = new TF1("const2","0.95",0,1);

     p1->SetLineColor(kRed);
     p1->SetTitle("combined posteriors;#epsilon;P(#epsilon|k,N)");
     p2->SetLineColor(kBlue);
     comb->SetLineColor(kGreen+2);

     TLegend* leg1 = new TLegend(0.12,0.65,0.5,0.85);
     leg1->AddEntry(p1,"k1 = 18, N1 = 26","l");
     leg1->AddEntry(p2,"k2 = 3, N2 = 10","l");
     leg1->AddEntry(comb,"combined: p1 = 0.6, p2=0.4","l");

     c1->cd(1);
     comb->Draw();
     p1->Draw("same");
     p2->Draw("same");
     leg1->Draw("same");
     c1->cd(2);
     const1->SetLineWidth(1);
     const2->SetLineWidth(1);
     TGraph* gr = (TGraph*)comb->DrawIntegral();
     gr->SetTitle("cumulative function of combined posterior with boundaries for cl = 95%;#epsilon;CDF");
     const1->Draw("same");
     const2->Draw("same");

     c1->cd(0);
     return c1;
}
End_Macro

**/
////////////////////////////////////////////////////////////////////
Double_t TEfficiency::Combine(Double_t& up,Double_t& low,Int_t n,
                              const Int_t* pass,const Int_t* total,
                              Double_t alpha, Double_t beta,
                              Double_t level,const Double_t* w,Option_t* opt)
{
   TString option(opt);
   option.ToLower();

   //LM:  new formula for combination
   // works only if alpha beta are the same always
   // the weights are normalized to w(i) -> N_eff w(i)/ Sum w(i)
   // i.e. w(i) -> Sum (w(i) / Sum (w(i)^2) * w(i)
   // norm = Sum (w(i) / Sum (w(i)^2)
   double ntot = 0;
   double ktot = 0;
   double sumw = 0;
   double sumw2 = 0;
   for (int i = 0; i < n ; ++i) {
      if(pass[i] > total[i]) {
         ::Error("TEfficiency::Combine","total events = %i < passed events %i",total[i],pass[i]);
         ::Info("TEfficiency::Combine","stop combining");
         return -1;
      }

      ntot += w[i] * total[i];
      ktot += w[i] * pass[i];
      sumw += w[i];
      sumw2 += w[i]*w[i];
      //mean += w[i] * (pass[i] + alpha[i])/(total[i] + alpha[i] + beta[i]);
   }
   double norm = sumw/sumw2;
   ntot *= norm;
   ktot *= norm;
   if(ktot > ntot) {
      ::Error("TEfficiency::Combine","total  = %f < passed  %f",ntot,ktot);
      ::Info("TEfficiency::Combine","stop combining");
      return -1;
   }

   double a = ktot + alpha;
   double b = ntot - ktot + beta;

   double mean = a/(a+b);
   double mode = BetaMode(a,b);


   Bool_t shortestInterval = option.Contains("sh") || ( option.Contains("mode") && !option.Contains("cent") );

   if (shortestInterval)
      BetaShortestInterval(level, a, b, low, up);
   else {
      low = BetaCentralInterval(level, a, b, false);
      up = BetaCentralInterval(level, a, b, true);
   }

   if (option.Contains("mode")) return mode;
   return mean;

}
////////////////////////////////////////////////////////////////////////////////
/// Combines a list of 1-dimensional TEfficiency objects
///
/// A TGraphAsymmErrors object is returned which contains the estimated
/// efficiency and its uncertainty for each bin.
/// If the combination fails, a zero pointer is returned.
///
/// At the moment the combining is only implemented for bayesian statistics.
///
/// \param[in] pList list containing TEfficiency objects which should be combined
///            only one-dimensional efficiencies are taken into account
/// \param[in] option
/// - s     : strict combining; only TEfficiency objects with the same beta
///           prior and the flag kIsBayesian == true are combined
///           If not specified the prior parameter of the first TEfficiency object is used
/// - v     : verbose mode; print information about combining
/// - cl=x  : set confidence level (0 < cl < 1). If not specified, the
///           confidence level of the first TEfficiency object is used.
/// - mode    Use mode of combined posterior as estimated value for the efficiency
/// - shortest: compute shortest interval (done by default if mode option is set)
/// - central: compute central interval (done by default if mode option is NOT set)
/// \param[in] n number of weights (has to be the number of one-dimensional
///            TEfficiency objects in pList)
///            If no weights are passed, the internal weights GetWeight() of
///            the given TEfficiency objects are used.
/// \param[in] w array of length n with weights for each TEfficiency object in
///           pList (w[0] correspond to pList->First ... w[n-1] -> pList->Last)
///           The weights do not have to be normalised.
///
/// For each bin the calculation is done by the Combine(double&, double& ...) method.

TGraphAsymmErrors* TEfficiency::Combine(TCollection* pList,Option_t* option,
                                        Int_t n,const Double_t* w)
{
   TString opt = option;
   opt.ToLower();

   //parameter of prior distribution, confidence level and normalisation factor
   Double_t alpha = -1;
   Double_t beta = -1;
   Double_t level = 0;

   //flags for combining
   Bool_t bStrict = false;
   Bool_t bOutput = false;
   Bool_t bWeights = false;
   //list of all information needed to weight and combine efficiencies
   std::vector<TH1*> vTotal;    vTotal.reserve(n);
   std::vector<TH1*> vPassed;   vPassed.reserve(n);
   std::vector<Double_t> vWeights;  vWeights.reserve(n);
   //    std::vector<Double_t> vAlpha;
   //    std::vector<Double_t> vBeta;

   if(opt.Contains("s")) {
      opt.ReplaceAll("s","");
      bStrict = true;
   }

   if(opt.Contains("v")) {
      opt.ReplaceAll("v","");
      bOutput = true;
   }

   if(opt.Contains("cl=")) {
      Ssiz_t pos = opt.Index("cl=") + 3;
      level = atof( opt(pos,opt.Length() ).Data() );
      if((level <= 0) || (level >= 1))
         level = 0;
      opt.ReplaceAll("cl=","");
   }

   //are weights explicitly given
   if(n && w) {
      bWeights = true;
      for(Int_t k = 0; k < n; ++k) {
         if(w[k] > 0)
            vWeights.push_back(w[k]);
         else {
            gROOT->Error("TEfficiency::Combine","invalid custom weight found w = %.2lf",w[k]);
            gROOT->Info("TEfficiency::Combine","stop combining");
            return 0;
         }
      }
   }

   TIter next(pList);
   TObject* obj = 0;
   TEfficiency* pEff = 0;
   while((obj = next())) {
      pEff = dynamic_cast<TEfficiency*>(obj);
      //is object a TEfficiency object?
      if(pEff) {
         if(pEff->GetDimension() > 1)
            continue;
         if(!level) level = pEff->GetConfidenceLevel();

         if(alpha<1) alpha = pEff->GetBetaAlpha();
         if(beta<1) beta = pEff->GetBetaBeta();

         //if strict combining, check priors, confidence level and statistic
         if(bStrict) {
            if(alpha != pEff->GetBetaAlpha())
               continue;
            if(beta != pEff->GetBetaBeta())
               continue;
            if(!pEff->UsesBayesianStat())
               continue;
         }

         vTotal.push_back(pEff->fTotalHistogram);
         vPassed.push_back(pEff->fPassedHistogram);

         //no weights given -> use weights of TEfficiency objects
         if(!bWeights)
            vWeights.push_back(pEff->fWeight);

         //strict combining -> using global prior
         //   if(bStrict) {
         //      vAlpha.push_back(alpha);
         //      vBeta.push_back(beta);
         //   }
         //   else {
         //      vAlpha.push_back(pEff->GetBetaAlpha());
         //      vBeta.push_back(pEff->GetBetaBeta());
         //   }
      }
   }

   //no TEfficiency objects found
   if(vTotal.empty()) {
      gROOT->Error("TEfficiency::Combine","no TEfficiency objects in given list");
      gROOT->Info("TEfficiency::Combine","stop combining");
      return 0;
   }

   //invalid number of custom weights
   if(bWeights && (n != (Int_t)vTotal.size())) {
      gROOT->Error("TEfficiency::Combine","number of weights n=%i differs from number of TEfficiency objects k=%i which should be combined",n,(Int_t)vTotal.size());
      gROOT->Info("TEfficiency::Combine","stop combining");
      return 0;
   }

   Int_t nbins_max = vTotal.at(0)->GetNbinsX();
   //check binning of all histograms
   for(UInt_t i=0; i<vTotal.size(); ++i) {
      if (!TEfficiency::CheckBinning(*vTotal.at(0),*vTotal.at(i)) )
         gROOT->Warning("TEfficiency::Combine","histograms have not the same binning -> results may be useless");
      if(vTotal.at(i)->GetNbinsX() < nbins_max) nbins_max = vTotal.at(i)->GetNbinsX();
   }

   //display information about combining
   if(bOutput) {
      gROOT->Info("TEfficiency::Combine","combining %i TEfficiency objects",(Int_t)vTotal.size());
      if(bWeights)
         gROOT->Info("TEfficiency::Combine","using custom weights");
      if(bStrict) {
         gROOT->Info("TEfficiency::Combine","using the following prior probability for the efficiency: P(e) ~ Beta(e,%.3lf,%.3lf)",alpha,beta);
      }
      else
         gROOT->Info("TEfficiency::Combine","using individual priors of each TEfficiency object");
      gROOT->Info("TEfficiency::Combine","confidence level = %.2lf",level);
   }

   //create TGraphAsymmErrors with efficiency
   std::vector<Double_t> x(nbins_max);
   std::vector<Double_t> xlow(nbins_max);
   std::vector<Double_t> xhigh(nbins_max);
   std::vector<Double_t> eff(nbins_max);
   std::vector<Double_t> efflow(nbins_max);
   std::vector<Double_t> effhigh(nbins_max);

   //parameters for combining:
   //number of objects
   Int_t num = vTotal.size();
   std::vector<Int_t> pass(num);
   std::vector<Int_t> total(num);

   //loop over all bins
   Double_t low = 0;
   Double_t up = 0;
   for(Int_t i=1; i <= nbins_max; ++i) {
      //the binning of the x-axis is taken from the first total histogram
      x[i-1] = vTotal.at(0)->GetBinCenter(i);
      xlow[i-1] = x[i-1] - vTotal.at(0)->GetBinLowEdge(i);
      xhigh[i-1] = vTotal.at(0)->GetBinWidth(i) - xlow[i-1];

      for(Int_t j = 0; j < num; ++j) {
         pass[j] = (Int_t)(vPassed.at(j)->GetBinContent(i) + 0.5);
         total[j] = (Int_t)(vTotal.at(j)->GetBinContent(i) + 0.5);
      }

      //fill efficiency and errors
      eff[i-1] = Combine(up,low,num,&pass[0],&total[0],alpha,beta,level,&vWeights[0],opt.Data());
      //did an error occurred ?
      if(eff[i-1] == -1) {
         gROOT->Error("TEfficiency::Combine","error occurred during combining");
         gROOT->Info("TEfficiency::Combine","stop combining");
         return 0;
      }
      efflow[i-1]= eff[i-1] - low;
      effhigh[i-1]= up - eff[i-1];
   }//loop over all bins

   TGraphAsymmErrors* gr = new TGraphAsymmErrors(nbins_max,&x[0],&eff[0],&xlow[0],&xhigh[0],&efflow[0],&effhigh[0]);

   return gr;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a graph.
///
///  Compute the closest distance of approach from point px,py to this line.
///  The distance is computed in pixels units.
///
/// Forward the call to the painted graph

Int_t TEfficiency::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (fPaintGraph) return fPaintGraph->DistancetoPrimitive(px,py);
   if (fPaintHisto) return fPaintHisto->DistancetoPrimitive(px,py);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Draws the current TEfficiency object
///
/// \param[in] opt
///  - 1-dimensional case: same options as TGraphAsymmErrors::Draw()
///     but as default "AP" is used
///  - 2-dimensional case: same options as TH2::Draw()
///  - 3-dimensional case: not yet supported
///
/// Specific TEfficiency drawing options:
/// - E0 - plot bins where the total number of passed events is zero
///      (the error interval will be [0,1] )

void TEfficiency::Draw(Option_t* opt)
{
   //check options
   TString option = opt;
   option.ToLower();
   // use by default "AP"
   if (option.IsNull() ) option = "ap";

   if(gPad && !option.Contains("same"))
      gPad->Clear();
   else {
      // add always "a" if not present
      if (!option.Contains("a") ) option += "a";
   }

   // add always p to the option
   if (!option.Contains("p") ) option += "p";


   AppendPad(option.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function is called when the drawn class is clicked with the locator
///  If Left button clicked on one of the line end points, this point
///     follows the cursor until button is released.
///
///  if Middle button clicked, the line is moved parallel to itself
///     until the button is released.
/// Forward the call to the underlying graph

void TEfficiency::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fPaintGraph) fPaintGraph->ExecuteEvent(event,px,py);
   else if (fPaintHisto) fPaintHisto->ExecuteEvent(event,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// This function is used for filling the two histograms.
///
/// \param[in] bPassed flag whether the current event passed the selection
///                 - true: both histograms are filled
///                 - false: only the total histogram is filled
/// \param[in] x x-value
/// \param[in] y y-value (use default=0 for 1-D efficiencies)
/// \param[in] z z-value (use default=0 for 2-D or 1-D efficiencies)

void TEfficiency::Fill(Bool_t bPassed,Double_t x,Double_t y,Double_t z)
{
   switch(GetDimension()) {
      case 1:
         fTotalHistogram->Fill(x);
         if(bPassed)
            fPassedHistogram->Fill(x);
         break;
      case 2:
         ((TH2*)(fTotalHistogram))->Fill(x,y);
         if(bPassed)
            ((TH2*)(fPassedHistogram))->Fill(x,y);
         break;
      case 3:
         ((TH3*)(fTotalHistogram))->Fill(x,y,z);
         if(bPassed)
            ((TH3*)(fPassedHistogram))->Fill(x,y,z);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
///This function is used for filling the two histograms with a weight.
///
/// \param[in] bPassed flag whether the current event passed the selection
///                 - true: both histograms are filled
///                 - false: only the total histogram is filled
/// \param[in] weight weight for the event
/// \param[in] x x-value
/// \param[in] y y-value (use default=0 for 1-D efficiencies)
/// \param[in] z z-value (use default=0 for 2-D or 1-D efficiencies)
///
/// Note: - this function will call SetUseWeightedEvents if it was not called by the user before

void TEfficiency::FillWeighted(Bool_t bPassed,Double_t weight,Double_t x,Double_t y,Double_t z)
{
   if(!TestBit(kUseWeights))
   {
      // Info("FillWeighted","call SetUseWeightedEvents() manually to ensure correct storage of sum of weights squared");
      SetUseWeightedEvents();
   }

   switch(GetDimension()) {
      case 1:
         fTotalHistogram->Fill(x,weight);
         if(bPassed)
            fPassedHistogram->Fill(x,weight);
         break;
      case 2:
         ((TH2*)(fTotalHistogram))->Fill(x,y,weight);
         if(bPassed)
            ((TH2*)(fPassedHistogram))->Fill(x,y,weight);
         break;
      case 3:
         ((TH3*)(fTotalHistogram))->Fill(x,y,z,weight);
         if(bPassed)
            ((TH3*)(fPassedHistogram))->Fill(x,y,z,weight);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the global bin number containing the given values
///
/// Note:
///
///   - values which belong to dimensions higher than the current dimension
///     of the TEfficiency object are ignored (i.e. for 1-dimensional
///     efficiencies only the x-value is considered)

Int_t TEfficiency::FindFixBin(Double_t x,Double_t y,Double_t z) const
{
   Int_t nx = fTotalHistogram->GetXaxis()->FindFixBin(x);
   Int_t ny = 0;
   Int_t nz = 0;

   switch(GetDimension()) {
      case 3: nz = fTotalHistogram->GetZaxis()->FindFixBin(z);
      case 2: ny = fTotalHistogram->GetYaxis()->FindFixBin(y);break;
   }

   return GetGlobalBin(nx,ny,nz);
}

////////////////////////////////////////////////////////////////////////////////
/// Fits the efficiency using the TBinomialEfficiencyFitter class
///
/// The resulting fit function is added to the list of associated functions.
///
/// Options:
/// - "+": previous fitted functions in the list are kept, by default
///   all functions in the list are deleted
/// - for more fitting options see TBinomialEfficiencyFitter::Fit

TFitResultPtr TEfficiency::Fit(TF1* f1,Option_t* opt)
{
   TString option = opt;
   option.ToLower();

   //replace existing functions in list with same name
   Bool_t bDeleteOld = true;
   if(option.Contains("+")) {
      option.ReplaceAll("+","");
      bDeleteOld = false;
   }

   TBinomialEfficiencyFitter Fitter(fPassedHistogram,fTotalHistogram);

   TFitResultPtr result = Fitter.Fit(f1,option.Data());

   //create copy which is appended to the list
   TF1* pFunc = new TF1(*f1);

   if(bDeleteOld) {
      TIter next(fFunctions);
      TObject* obj = 0;
      while((obj = next())) {
         if(obj->InheritsFrom(TF1::Class())) {
            fFunctions->Remove(obj);
            delete obj;
         }
      }
   }

   // create list if necessary
   if(!fFunctions)
      fFunctions = new TList();

   fFunctions->Add(pFunc);

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a cloned version of fPassedHistogram
///
/// Notes:
///       - The histogram is filled with unit weights. You might want to scale
///         it with the global weight GetWeight().
///       - The returned object is owned by the user who has to care about the
///         deletion of the new TH1 object.
///       - This histogram is by default NOT attached to the current directory
///         to avoid duplication of data. If you want to store it automatically
///         during the next TFile::Write() command, you have to attach it to
///         the corresponding directory.
///
/// ~~~~~~~{.cpp}
///  TFile* pFile = new TFile("passed.root","update");
///  TEfficiency* pEff = (TEfficiency*)gDirectory->Get("my_eff");
///  TH1* copy = pEff->GetCopyPassedHisto();
///  copy->SetDirectory(gDirectory);
///  pFile->Write();
/// ~~~~~~~

TH1* TEfficiency::GetCopyPassedHisto() const
{
   // do not add cloned histogram to gDirectory
   TDirectory::TContext ctx(nullptr);
   TH1* tmp = (TH1*)(fPassedHistogram->Clone());

   return tmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a cloned version of fTotalHistogram
///
/// Notes:
///       - The histogram is filled with unit weights. You might want to scale
///         it with the global weight GetWeight().
///       - The returned object is owned by the user who has to care about the
///         deletion of the new TH1 object.
///       - This histogram is by default NOT attached to the current directory
///         to avoid duplication of data. If you want to store it automatically
///         during the next TFile::Write() command, you have to attach it to
///         the corresponding directory.
///
/// ~~~~~~~{.cpp}
///  TFile* pFile = new TFile("total.root","update");
///  TEfficiency* pEff = (TEfficiency*)gDirectory->Get("my_eff");
///  TH1* copy = pEff->GetCopyTotalHisto();
///  copy->SetDirectory(gDirectory);
///  pFile->Write();
/// ~~~~~~~

TH1* TEfficiency::GetCopyTotalHisto() const
{
   // do not add cloned histogram to gDirectory
   TDirectory::TContext ctx(nullptr);
   TH1* tmp = (TH1*)(fTotalHistogram->Clone());

   return tmp;
}

////////////////////////////////////////////////////////////////////////////////
///returns the dimension of the current TEfficiency object

Int_t TEfficiency::GetDimension() const
{
   return fTotalHistogram->GetDimension();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the efficiency in the given global bin
///
/// Note:
///      - The estimated efficiency depends on the chosen statistic option:
///        for frequentist ones:
///        \f$ \hat{\varepsilon} = \frac{passed}{total} \f$
///        for bayesian ones the expectation value of the resulting posterior
///        distribution is returned:
///        \f$ \hat{\varepsilon} = \frac{passed + \alpha}{total + \alpha + \beta} \f$
///        If the bit kPosteriorMode is set (or the method TEfficiency::UsePosteriorMode() has been called ) the
///        mode (most probable value) of the posterior is returned:
///        \f$ \hat{\varepsilon} = \frac{passed + \alpha -1}{total + \alpha + \beta -2} \f$
///       - If the denominator is equal to 0, an efficiency of 0 is returned.
///       - When \f$ passed + \alpha < 1 \f$ or \f$ total - passed + \beta < 1 \f$ the above
///        formula for the mode is not valid. In these cases values the estimated efficiency is 0 or 1.

Double_t TEfficiency::GetEfficiency(Int_t bin) const
{
   Double_t total = fTotalHistogram->GetBinContent(bin);
   Double_t passed = fPassedHistogram->GetBinContent(bin);

   if(TestBit(kIsBayesian)) {

      // parameters for the beta prior distribution
      Double_t alpha = TestBit(kUseBinPrior) ? GetBetaAlpha(bin) : GetBetaAlpha();
      Double_t beta  = TestBit(kUseBinPrior) ? GetBetaBeta(bin)  : GetBetaBeta();

      Double_t aa,bb;
      if(TestBit(kUseWeights))
      {
         Double_t tw =  fTotalHistogram->GetBinContent(bin);
         Double_t tw2 = fTotalHistogram->GetSumw2()->At(bin);
         Double_t pw =  fPassedHistogram->GetBinContent(bin);

         if (tw2 <= 0 ) return pw/tw;

         // tw/tw2 renormalize the weights
         double norm = tw/tw2;
         aa =  pw * norm + alpha;
         bb =  (tw - pw) * norm + beta;
      }
      else
      {
         aa = passed + alpha;
         bb = total - passed + beta;
      }

      if (!TestBit(kPosteriorMode) )
         return BetaMean(aa,bb);
      else
         return BetaMode(aa,bb);

   }
   else
      return (total)? ((Double_t)passed)/total : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the lower error on the efficiency in the given global bin
///
/// The result depends on the current confidence level fConfLevel and the
/// chosen statistic option fStatisticOption. See SetStatisticOption(Int_t) for
/// more details.
///
/// Note: If the histograms are filled with weights, only bayesian methods and the
///       normal approximation are supported.

Double_t TEfficiency::GetEfficiencyErrorLow(Int_t bin) const
{
   Double_t total = fTotalHistogram->GetBinContent(bin);
   Double_t passed = fPassedHistogram->GetBinContent(bin);

   Double_t eff = GetEfficiency(bin);

   // check whether weights have been used
   if(TestBit(kUseWeights))
   {
      Double_t tw =  fTotalHistogram->GetBinContent(bin);
      Double_t tw2 = fTotalHistogram->GetSumw2()->At(bin);
      Double_t pw =  fPassedHistogram->GetBinContent(bin);
      Double_t pw2 = fPassedHistogram->GetSumw2()->At(bin);

      if(TestBit(kIsBayesian))
      {
         Double_t alpha = TestBit(kUseBinPrior) ? GetBetaAlpha(bin) : GetBetaAlpha();
         Double_t beta  = TestBit(kUseBinPrior) ? GetBetaBeta(bin)  : GetBetaBeta();

         if (tw2 <= 0) return 0;

         // tw/tw2 renormalize the weights
         Double_t norm = tw/tw2;
         Double_t aa =  pw * norm + alpha;
         Double_t bb =  (tw - pw) * norm + beta;
         Double_t low = 0;
         Double_t upper = 1;
         if(TestBit(kShortestInterval)) {
            TEfficiency::BetaShortestInterval(fConfLevel,aa,bb,low,upper);
         }
         else {
            low = TEfficiency::BetaCentralInterval(fConfLevel,aa,bb,false);
         }

         return eff - low;
      }
      else
      {
         if(fStatisticOption != kFNormal)
         {
            Warning("GetEfficiencyErrorLow","frequentist confidence intervals for weights are only supported by the normal approximation");
            Info("GetEfficiencyErrorLow","setting statistic option to kFNormal");
            const_cast<TEfficiency*>(this)->SetStatisticOption(kFNormal);
         }

         Double_t variance = ( pw2 * (1. - 2 * eff) + tw2 * eff *eff ) / ( tw * tw) ;
         Double_t sigma = sqrt(variance);

         Double_t prob = 0.5 * (1.- fConfLevel);
         Double_t delta = ROOT::Math::normal_quantile_c(prob, sigma);

         // avoid to return errors which makes eff-err < 0
         return (eff - delta < 0) ? eff : delta;
      }
   }
   else
   {
      if(TestBit(kIsBayesian))
      {
         // parameters for the beta prior distribution
         Double_t alpha = TestBit(kUseBinPrior) ? GetBetaAlpha(bin) : GetBetaAlpha();
         Double_t beta  = TestBit(kUseBinPrior) ? GetBetaBeta(bin)  : GetBetaBeta();
         return (eff - Bayesian(total,passed,fConfLevel,alpha,beta,false,TestBit(kShortestInterval)));
      }
      else
         return (eff - fBoundary(total,passed,fConfLevel,false));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the upper error on the efficiency in the given global bin
///
/// The result depends on the current confidence level fConfLevel and the
/// chosen statistic option fStatisticOption. See SetStatisticOption(Int_t) for
/// more details.
///
/// Note: If the histograms are filled with weights, only bayesian methods and the
///       normal approximation are supported.

Double_t TEfficiency::GetEfficiencyErrorUp(Int_t bin) const
{
   Double_t total = fTotalHistogram->GetBinContent(bin);
   Double_t passed = fPassedHistogram->GetBinContent(bin);

   Double_t eff = GetEfficiency(bin);

   // check whether weights have been used
   if(TestBit(kUseWeights))
   {
      Double_t tw =  fTotalHistogram->GetBinContent(bin);
      Double_t tw2 = fTotalHistogram->GetSumw2()->At(bin);
      Double_t pw =  fPassedHistogram->GetBinContent(bin);
      Double_t pw2 = fPassedHistogram->GetSumw2()->At(bin);

      if(TestBit(kIsBayesian))
      {
         Double_t alpha = TestBit(kUseBinPrior) ? GetBetaAlpha(bin) : GetBetaAlpha();
         Double_t beta  = TestBit(kUseBinPrior) ? GetBetaBeta(bin)  : GetBetaBeta();

         if (tw2 <= 0) return 0;

         // tw/tw2 renormalize the weights
         Double_t norm = tw/tw2;
         Double_t aa =  pw * norm + alpha;
         Double_t bb =  (tw - pw) * norm + beta;
         Double_t low = 0;
         Double_t upper = 1;
         if(TestBit(kShortestInterval)) {
            TEfficiency::BetaShortestInterval(fConfLevel,aa,bb,low,upper);
         }
         else {
            upper = TEfficiency::BetaCentralInterval(fConfLevel,aa,bb,true);
         }

         return upper - eff;
      }
      else
      {
         if(fStatisticOption != kFNormal)
         {
            Warning("GetEfficiencyErrorUp","frequentist confidence intervals for weights are only supported by the normal approximation");
            Info("GetEfficiencyErrorUp","setting statistic option to kFNormal");
            const_cast<TEfficiency*>(this)->SetStatisticOption(kFNormal);
         }

         Double_t variance = ( pw2 * (1. - 2 * eff) + tw2 * eff *eff ) / ( tw * tw) ;
         Double_t sigma = sqrt(variance);

         Double_t prob = 0.5 * (1.- fConfLevel);
         Double_t delta = ROOT::Math::normal_quantile_c(prob, sigma);

         return (eff + delta > 1) ? 1.-eff : delta;
      }
   }
   else
   {
      if(TestBit(kIsBayesian))
      {
         // parameters for the beta prior distribution
         Double_t alpha = TestBit(kUseBinPrior) ? GetBetaAlpha(bin) : GetBetaAlpha();
         Double_t beta  = TestBit(kUseBinPrior) ? GetBetaBeta(bin)  : GetBetaBeta();
         return (Bayesian(total,passed,fConfLevel,alpha,beta,true,TestBit(kShortestInterval)) - eff);
      }
      else
         return fBoundary(total,passed,fConfLevel,true) - eff;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the global bin number which can be used as argument for the
/// following functions:
///
/// - GetEfficiency(bin), GetEfficiencyErrorLow(bin), GetEfficiencyErrorUp(bin)
/// - SetPassedEvents(bin), SetTotalEvents(bin)
///
/// see TH1::GetBin() for conventions on numbering bins

Int_t TEfficiency::GetGlobalBin(Int_t binx,Int_t biny,Int_t binz) const
{
   return fTotalHistogram->GetBin(binx,biny,binz);
}

////////////////////////////////////////////////////////////////////////////////

TList* TEfficiency::GetListOfFunctions()
{
   return (fFunctions) ? fFunctions : fFunctions = new TList();
}

////////////////////////////////////////////////////////////////////////////////
/// Merges the TEfficiency objects in the given list to the given
/// TEfficiency object using the operator+=(TEfficiency&)
///
/// The merged result is stored in the current object. The statistic options and
/// the confidence level are taken from the current object.
///
/// This function should be used when all TEfficiency objects correspond to
/// the same process.
///
/// The new weight is set according to:
/// \f$  \frac{1}{w_{new}} = \sum_{i} \frac{1}{w_{i}} \f$

Long64_t TEfficiency::Merge(TCollection* pList)
{
   if(!pList->IsEmpty()) {
      TIter next(pList);
      TObject* obj = 0;
      TEfficiency* pEff = 0;
      while((obj = next())) {
         pEff = dynamic_cast<TEfficiency*>(obj);
         if(pEff) {
            *this += *pEff;
         }
      }
   }
   return (Long64_t)fTotalHistogram->GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/**
Returns the confidence limits for the efficiency supposing that the
efficiency follows a normal distribution with the rms below

\param[in] total number of total events
\param[in] passed 0 <= number of passed events <= total
\param[in] level  confidence level
\param[in] bUpper
                - true  - upper boundary is returned
                - false - lower boundary is returned

Calculation:

\f{eqnarray*}{
      \hat{\varepsilon} &=& \frac{passed}{total}\\
      \sigma_{\varepsilon} &=& \sqrt{\frac{\hat{\varepsilon} (1 - \hat{\varepsilon})}{total}}\\
      \varepsilon_{low} &=& \hat{\varepsilon} \pm \Phi^{-1}(\frac{level}{2},\sigma_{\varepsilon})
\f}
*/

Double_t TEfficiency::Normal(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   Double_t alpha = (1.0 - level)/2;
   if (total == 0) return (bUpper) ? 1 : 0;
   Double_t average = passed / total;
   Double_t sigma = std::sqrt(average * (1 - average) / total);
   Double_t delta = ROOT::Math::normal_quantile(1 - alpha,sigma);

   if(bUpper)
      return ((average + delta) > 1) ? 1.0 : (average + delta);
   else
      return ((average - delta) < 0) ? 0.0 : (average - delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Adds the histograms of another TEfficiency object to current histograms
///
/// The statistic options and the confidence level remain unchanged.
///
/// fTotalHistogram += rhs.fTotalHistogram;
/// fPassedHistogram += rhs.fPassedHistogram;
///
/// calculates a new weight:
/// current weight of this TEfficiency object = \f$ w_{1} \f$
/// weight of rhs = \f$ w_{2} \f$
/// \f$ w_{new} = \frac{w_{1} \times w_{2}}{w_{1} + w_{2}} \f$

TEfficiency& TEfficiency::operator+=(const TEfficiency& rhs)
{

   if (fTotalHistogram == 0 && fPassedHistogram == 0) {
      // efficiency is empty just copy it over
      *this = rhs;
      return *this;
   }
   else if (fTotalHistogram == 0 || fPassedHistogram == 0) {
      Fatal("operator+=","Adding to a non consistent TEfficiency object which has not a total or a passed histogram ");
      return *this;
   }

   if (rhs.fTotalHistogram == 0 && rhs.fPassedHistogram == 0 ) {
      Warning("operator+=","no operation: adding an empty object");
      return *this;
   }
   else  if (rhs.fTotalHistogram == 0  || rhs.fPassedHistogram == 0 ) {
      Fatal("operator+=","Adding a non consistent TEfficiency object which has not a total or a passed histogram ");
      return *this;
   }

   fTotalHistogram->ResetBit(TH1::kIsAverage);
   fPassedHistogram->ResetBit(TH1::kIsAverage);

   fTotalHistogram->Add(rhs.fTotalHistogram);
   fPassedHistogram->Add(rhs.fPassedHistogram);

   SetWeight((fWeight * rhs.GetWeight())/(fWeight + rhs.GetWeight()));

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator
///
/// The histograms, statistic option, confidence level, weight and paint styles
/// of rhs are copied to the this TEfficiency object.
///
/// Note: - The list of associated functions is not copied. After this
///         operation the list of associated functions is empty.

TEfficiency& TEfficiency::operator=(const TEfficiency& rhs)
{
   if(this != &rhs)
   {
      //statistic options
      SetStatisticOption(rhs.GetStatisticOption());
      SetConfidenceLevel(rhs.GetConfidenceLevel());
      SetBetaAlpha(rhs.GetBetaAlpha());
      SetBetaBeta(rhs.GetBetaBeta());
      SetWeight(rhs.GetWeight());

      //associated list of functions
      if(fFunctions)
         fFunctions->Delete();

      //copy histograms
      delete fTotalHistogram;
      delete fPassedHistogram;

      // do not add cloned histogram to gDirectory
      {
         TDirectory::TContext ctx(nullptr);
         fTotalHistogram = (TH1*)(rhs.fTotalHistogram->Clone());
         fPassedHistogram = (TH1*)(rhs.fPassedHistogram->Clone());
      }
      //delete temporary paint objects
      delete fPaintHisto;
      delete fPaintGraph;
      fPaintHisto = 0;
      fPaintGraph = 0;

      //copy style
      rhs.TAttLine::Copy(*this);
      rhs.TAttFill::Copy(*this);
      rhs.TAttMarker::Copy(*this);
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Paints this TEfficiency object
///
/// For details on the possible option see Draw(Option_t*)
///
/// Note for 1D classes
/// In 1D the TEfficiency uses a TGraphAsymmErrors for drawing
/// The TGraph is created only the first time Paint is used. The user can manipulate the
/// TGraph via the method TEfficiency::GetPaintedGraph()
/// The TGraph creates behing an histogram for the axis. The histogram is created also only the first time.
/// If the axis needs to be updated because in the meantime the class changed use this trick
/// which will trigger a re-calculation of the axis of the graph
/// TEfficiency::GetPaintedGraph()->Set(0)
///
/// Note that in order to access the painted graph via GetPaintedGraph() you need either to call Paint or better
/// gPad->Update();
///

void TEfficiency::Paint(const Option_t* opt)
{


   if(!gPad)
      return;


   //use TGraphAsymmErrors for painting
   if(GetDimension() == 1) {
      if(!fPaintGraph) {
         fPaintGraph = CreateGraph(opt);
      }
      else
         // update existing graph already created
         FillGraph(fPaintGraph, opt);

      //paint graph

      fPaintGraph->Paint(opt);

      //paint all associated functions
      if(fFunctions) {
         //paint box with fit parameters
         //the fit statistics will be painted if gStyle->SetOptFit(1) has been
         // called by the user
         TIter next(fFunctions);
         TObject* obj = 0;
         while((obj = next())) {
            if(obj->InheritsFrom(TF1::Class())) {
               fPaintGraph->PaintStats((TF1*)obj);
               ((TF1*)obj)->Paint("sameC");
            }
         }
      }

      return;
   }

   //use TH2 for painting
   if(GetDimension() == 2) {
      if(!fPaintHisto) {
         fPaintHisto = CreateHistogram();
      }
      else
         FillHistogram(fPaintHisto);

      //paint histogram
      fPaintHisto->Paint(opt);
      return;
   }
   Warning("Paint","Painting 3D efficiency is not implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Have histograms fixed bins along each axis?

void TEfficiency::SavePrimitive(std::ostream& out,Option_t* opt)
{
   Bool_t equi_bins = true;

   //indentation
   TString indent = "   ";
   //names for arrays containing the bin edges
   //static counter needed if more objects are saved
   static Int_t naxis = 0;
   TString sxaxis="xAxis",syaxis="yAxis",szaxis="zAxis";

   //note the missing break statements!
   switch(GetDimension()) {
      case 3:
         equi_bins = equi_bins && !fTotalHistogram->GetZaxis()->GetXbins()->fArray
         && !fTotalHistogram->GetZaxis()->GetXbins()->fN;
      case 2:
         equi_bins = equi_bins && !fTotalHistogram->GetYaxis()->GetXbins()->fArray
         && !fTotalHistogram->GetYaxis()->GetXbins()->fN;
      case 1:
         equi_bins = equi_bins && !fTotalHistogram->GetXaxis()->GetXbins()->fArray
         && !fTotalHistogram->GetXaxis()->GetXbins()->fN;
   }

   //create arrays containing the variable binning
   if(!equi_bins) {
      Int_t i;
      ++naxis;
      sxaxis += naxis;
      syaxis += naxis;
      szaxis += naxis;
      //x axis
      out << indent << "Double_t " << sxaxis << "["
      << fTotalHistogram->GetXaxis()->GetXbins()->fN << "] = {";
      for (i = 0; i < fTotalHistogram->GetXaxis()->GetXbins()->fN; ++i) {
         if (i != 0) out << ", ";
         out << fTotalHistogram->GetXaxis()->GetXbins()->fArray[i];
      }
      out << "}; " << std::endl;
      //y axis
      if(GetDimension() > 1) {
         out << indent << "Double_t " << syaxis << "["
         << fTotalHistogram->GetYaxis()->GetXbins()->fN << "] = {";
         for (i = 0; i < fTotalHistogram->GetYaxis()->GetXbins()->fN; ++i) {
            if (i != 0) out << ", ";
            out << fTotalHistogram->GetYaxis()->GetXbins()->fArray[i];
         }
         out << "}; " << std::endl;
      }
      //z axis
      if(GetDimension() > 2) {
         out << indent << "Double_t " << szaxis << "["
         << fTotalHistogram->GetZaxis()->GetXbins()->fN << "] = {";
         for (i = 0; i < fTotalHistogram->GetZaxis()->GetXbins()->fN; ++i) {
            if (i != 0) out << ", ";
            out << fTotalHistogram->GetZaxis()->GetXbins()->fArray[i];
         }
         out << "}; " << std::endl;
      }
   }//creating variable binning

   //TEfficiency pointer has efficiency name + counter
   static Int_t eff_count = 0;
   ++eff_count;
   TString eff_name = GetName();
   eff_name += eff_count;

   const char* name = eff_name.Data();

   //construct TEfficiency object
   const char quote = '"';
   out << indent << std::endl;
   out << indent << ClassName() << " * " << name << " = new " << ClassName()
   << "(" << quote << GetName() << quote << "," << quote
   << GetTitle() << quote <<",";
   //fixed bin size -> use n,min,max constructor
   if(equi_bins) {
      out << fTotalHistogram->GetXaxis()->GetNbins() << ","
      << fTotalHistogram->GetXaxis()->GetXmin() << ","
      << fTotalHistogram->GetXaxis()->GetXmax();
      if(GetDimension() > 1) {
         out << "," << fTotalHistogram->GetYaxis()->GetNbins() << ","
         << fTotalHistogram->GetYaxis()->GetXmin() << ","
         << fTotalHistogram->GetYaxis()->GetXmax();
      }
      if(GetDimension() > 2) {
         out << "," << fTotalHistogram->GetZaxis()->GetNbins() << ","
         << fTotalHistogram->GetZaxis()->GetXmin() << ","
         << fTotalHistogram->GetZaxis()->GetXmax();
      }
   }
   //variable bin size -> use n,*bins constructor
   else {
      out << fTotalHistogram->GetXaxis()->GetNbins() << "," << sxaxis;
      if(GetDimension() > 1)
         out << "," << fTotalHistogram->GetYaxis()->GetNbins() << ","
         << syaxis;
      if(GetDimension() > 2)
         out << "," << fTotalHistogram->GetZaxis()->GetNbins() << ","
         << szaxis;
   }
   out << ");" << std::endl;
   out << indent << std::endl;

   //set statistic options
   out << indent << name << "->SetConfidenceLevel(" << fConfLevel << ");"
   << std::endl;
   out << indent << name << "->SetBetaAlpha(" << fBeta_alpha << ");"
   << std::endl;
   out << indent << name << "->SetBetaBeta(" << fBeta_beta << ");" << std::endl;
   out << indent << name << "->SetWeight(" << fWeight << ");" << std::endl;
   out << indent << name << "->SetStatisticOption(" << fStatisticOption << ");"
   << std::endl;
   out << indent << name << "->SetPosteriorMode(" << TestBit(kPosteriorMode) << ");" << std::endl;
   out << indent << name << "->SetShortestInterval(" << TestBit(kShortestInterval) << ");" << std::endl;
   if(TestBit(kUseWeights))
      out << indent << name << "->SetUseWeightedEvents();" << std::endl;

   // save bin-by-bin prior parameters
   for(unsigned int i = 0; i < fBeta_bin_params.size(); ++i)
   {
      out << indent << name << "->SetBetaBinParameters(" << i << "," << fBeta_bin_params.at(i).first
      << "," << fBeta_bin_params.at(i).second << ");" << std::endl;
   }

   //set bin contents
   Int_t nbins = fTotalHistogram->GetNbinsX() + 2;
   if(GetDimension() > 1)
      nbins *= fTotalHistogram->GetNbinsY() + 2;
   if(GetDimension() > 2)
      nbins *= fTotalHistogram->GetNbinsZ() + 2;

   //important: set first total number than passed number
   for(Int_t i = 0; i < nbins; ++i) {
      out << indent << name <<"->SetTotalEvents(" << i << "," <<
      fTotalHistogram->GetBinContent(i) << ");" << std::endl;
      out << indent << name <<"->SetPassedEvents(" << i << "," <<
      fPassedHistogram->GetBinContent(i) << ");" << std::endl;
   }

   //save list of functions
   TIter next(fFunctions);
   TObject* obj = 0;
   while((obj = next())) {
      obj->SavePrimitive(out,"nodraw");
      if(obj->InheritsFrom(TF1::Class())) {
         out << indent << name << "->GetListOfFunctions()->Add("
         << obj->GetName() << ");" << std::endl;
      }
   }

   //set style
   SaveFillAttributes(out,name);
   SaveLineAttributes(out,name);
   SaveMarkerAttributes(out,name);

   //draw TEfficiency object
   TString option = opt;
   option.ToLower();
   if (!option.Contains("nodraw"))
      out<< indent << name<< "->Draw(" << quote << opt << quote << ");"
      << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the shape parameter &alpha;
///
/// The prior probability of the efficiency is given by the beta distribution:
/// \f[
///   f(\varepsilon;\alpha;\beta) = \frac{1}{B(\alpha,\beta)} \varepsilon^{\alpha-1} (1 - \varepsilon)^{\beta-1}
/// \f]
///
/// Note: - both shape parameters have to be positive (i.e. > 0)

void TEfficiency::SetBetaAlpha(Double_t alpha)
{
   if(alpha > 0)
      fBeta_alpha = alpha;
   else
      Warning("SetBetaAlpha(Double_t)","invalid shape parameter %.2lf",alpha);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the shape parameter &beta;
///
/// The prior probability of the efficiency is given by the beta distribution:
/// \f[
///   f(\varepsilon;\alpha,\beta) = \frac{1}{B(\alpha,\beta)} \varepsilon^{\alpha-1} (1 - \varepsilon)^{\beta-1}
/// \f]
///
/// Note: - both shape parameters have to be positive (i.e. > 0)

void TEfficiency::SetBetaBeta(Double_t beta)
{
   if(beta > 0)
      fBeta_beta = beta;
   else
      Warning("SetBetaBeta(Double_t)","invalid shape parameter %.2lf",beta);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets different  shape parameter &alpha; and &beta;
/// for the prior distribution for each bin. By default the global parameter are used if they are not set
/// for the specific bin
/// The prior probability of the efficiency is given by the beta distribution:
/// \f[
///   f(\varepsilon;\alpha;\beta) = \frac{1}{B(\alpha,\beta)} \varepsilon^{\alpha-1} (1 - \varepsilon)^{\beta-1}
/// \f]
///
/// Note:
/// - both shape parameters have to be positive (i.e. > 0)
/// - bin gives the global bin number (cf. GetGlobalBin)

void TEfficiency::SetBetaBinParameters(Int_t bin, Double_t alpha, Double_t beta)
{
   if (!fPassedHistogram || !fTotalHistogram) return;
   TH1 * h1 = fTotalHistogram;
   // doing this I get h1->fN which is available only for a TH1D
   UInt_t n = h1->GetBin(h1->GetNbinsX()+1, h1->GetNbinsY()+1, h1->GetNbinsZ()+1 ) + 1;

   // in case vector is not created do with default alpha, beta params
   if (fBeta_bin_params.size() != n )
      fBeta_bin_params = std::vector<std::pair<Double_t, Double_t> >(n, std::make_pair(fBeta_alpha, fBeta_beta) );

   // vector contains also values for under/overflows
   fBeta_bin_params[bin] = std::make_pair(alpha,beta);
   SetBit(kUseBinPrior,true);

}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, Double_t xmin, Double_t xmax)
{
   if (GetDimension() != 1) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xmin,xmax);
   fTotalHistogram->SetBins(nx,xmin,xmax);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, const Double_t *xBins)
{
   if (GetDimension() != 1) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xBins);
   fTotalHistogram->SetBins(nx,xBins);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
   if (GetDimension() != 2) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xmin,xmax,ny,ymin,ymax);
   fTotalHistogram->SetBins(nx,xmin,xmax,ny,ymin,ymax);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins)
{
   if (GetDimension() != 2) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xBins,ny,yBins);
   fTotalHistogram->SetBins(nx,xBins,ny,yBins);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax,
                            Int_t nz, Double_t zmin, Double_t zmax)
{
   if (GetDimension() != 3) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax);
   fTotalHistogram->SetBins (nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bins for the underlined passed and total histograms
/// If the class have been already filled the previous contents will be lost

Bool_t TEfficiency::SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins, Int_t nz,
                            const Double_t *zBins )
{
   if (GetDimension() != 3) {
      Error("SetBins","Using wrong SetBins function for a %d-d histogram",GetDimension());
      return kFALSE;
   }
   if (fTotalHistogram->GetEntries() != 0 ) {
      Warning("SetBins","Histogram entries will be lost after SetBins");
      fPassedHistogram->Reset();
      fTotalHistogram->Reset();
   }
   fPassedHistogram->SetBins(nx,xBins,ny,yBins,nz,zBins);
   fTotalHistogram->SetBins(nx,xBins,ny,yBins,nz,zBins);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the confidence level (0 < level < 1)
/// The default value is 1-sigma :~ 0.683

void TEfficiency::SetConfidenceLevel(Double_t level)
{
   if((level > 0) && (level < 1))
      fConfLevel = level;
   else
      Warning("SetConfidenceLevel(Double_t)","invalid confidence level %.2lf",level);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the directory holding this TEfficiency object
///
/// A reference to this TEfficiency object is removed from the current
/// directory (if it exists) and a new reference to this TEfficiency object is
/// added to the given directory.
///
/// Notes: - If the given directory is 0, the TEfficiency object does not
///         belong to any directory and will not be written to file during the
///         next TFile::Write() command.

void TEfficiency::SetDirectory(TDirectory* dir)
{
   if(fDirectory == dir)
      return;
   if(fDirectory)
      fDirectory->Remove(this);
   fDirectory = dir;
   if(fDirectory)
      fDirectory->Append(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the name
///
/// Note: The names of the internal histograms are set to "name + _total" and
///      "name + _passed" respectively.

void TEfficiency::SetName(const char* name)
{
   TNamed::SetName(name);

   //setting the names (appending the correct ending)
   TString name_total = name + TString("_total");
   TString name_passed = name + TString("_passed");
   fTotalHistogram->SetName(name_total);
   fPassedHistogram->SetName(name_passed);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the number of passed events in the given global bin
///
/// returns "true" if the number of passed events has been updated
/// otherwise "false" ist returned
///
/// Note: - requires: 0 <= events <= fTotalHistogram->GetBinContent(bin)

Bool_t TEfficiency::SetPassedEvents(Int_t bin,Int_t events)
{
   if(events <= fTotalHistogram->GetBinContent(bin)) {
      fPassedHistogram->SetBinContent(bin,events);
      return true;
   }
   else {
      Error("SetPassedEvents(Int_t,Int_t)","total number of events (%.1lf) in bin %i is less than given number of passed events %i",fTotalHistogram->GetBinContent(bin),bin,events);
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the histogram containing the passed events
///
/// The given histogram is cloned and stored internally as histogram containing
/// the passed events. The given histogram has to be consistent with the current
/// fTotalHistogram (see CheckConsistency(const TH1&,const TH1&)).
/// The method returns whether the fPassedHistogram has been replaced (true) or
/// not (false).
///
/// Note: The list of associated functions fFunctions is cleared.
///
/// Option:
/// - "f": force the replacement without checking the consistency
///               This can lead to inconsistent histograms and useless results
///               or unexpected behaviour. But sometimes it might be the only
///               way to change the histograms. If you use this option, you
///               should ensure that the fTotalHistogram is replaced by a
///               consistent one (with respect to rPassed) as well.

Bool_t TEfficiency::SetPassedHistogram(const TH1& rPassed,Option_t* opt)
{
   TString option = opt;
   option.ToLower();

   Bool_t bReplace = option.Contains("f");

   if(!bReplace)
      bReplace = CheckConsistency(rPassed,*fTotalHistogram);

   if(bReplace) {
      delete fPassedHistogram;
      // do not add cloned histogram to gDirectory
      {
         TDirectory::TContext ctx(nullptr);
         fPassedHistogram = (TH1*)(rPassed.Clone());
         fPassedHistogram->SetNormFactor(0);
      }

      if(fFunctions)
         fFunctions->Delete();

      //check whether both histograms are filled with weights
      bool useWeights = CheckWeights(rPassed,*fTotalHistogram);

      SetUseWeightedEvents(useWeights);

      return true;
   }
   else
      return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the statistic option which affects the calculation of the confidence interval
///
/// Options:
/// - kFCP (=0)(default): using the Clopper-Pearson interval (recommended by PDG)
///                       sets kIsBayesian = false
///                       see also ClopperPearson
/// - kFNormal   (=1)   : using the normal approximation
///                       sets kIsBayesian = false
///                       see also Normal
/// - kFWilson   (=2)   : using the Wilson interval
///                       sets kIsBayesian = false
///                       see also Wilson
/// - kFAC       (=3)   : using the Agresti-Coull interval
///                       sets kIsBayesian = false
///                       see also AgrestiCoull
/// - kFFC       (=4)   : using the Feldman-Cousins frequentist method
///                       sets kIsBayesian = false
///                       see also FeldmanCousins
/// - kBJeffrey  (=5)   : using the Jeffrey interval
///                       sets kIsBayesian = true, fBeta_alpha = 0.5 and fBeta_beta = 0.5
///                       see also Bayesian
/// - kBUniform  (=6)   : using a uniform prior
///                       sets kIsBayesian = true, fBeta_alpha = 1 and fBeta_beta = 1
///                       see also Bayesian
/// - kBBayesian (=7)   : using a custom prior defined by fBeta_alpha and fBeta_beta
///                      sets kIsBayesian = true
///                      see also Bayesian
/// - kMidP (=8)       : using the Lancaster Mid-P method
///                      sets kIsBayesian = false


void TEfficiency::SetStatisticOption(EStatOption option)
{
   fStatisticOption = option;

   switch(option)
   {
      case kFCP:
         fBoundary = &ClopperPearson;
         SetBit(kIsBayesian,false);
         break;
      case kFNormal:
         fBoundary = &Normal;
         SetBit(kIsBayesian,false);
         break;
      case kFWilson:
         fBoundary = &Wilson;
         SetBit(kIsBayesian,false);
         break;
      case kFAC:
         fBoundary = &AgrestiCoull;
         SetBit(kIsBayesian,false);
         break;
      case kFFC:
         fBoundary = &FeldmanCousins;
         SetBit(kIsBayesian,false);
         break;
      case kMidP:
         fBoundary = &MidPInterval;
         SetBit(kIsBayesian,false);
         break;
      case kBJeffrey:
         fBeta_alpha = 0.5;
         fBeta_beta = 0.5;
         SetBit(kIsBayesian,true);
         SetBit(kUseBinPrior,false);
         break;
      case kBUniform:
         fBeta_alpha = 1;
         fBeta_beta = 1;
         SetBit(kIsBayesian,true);
         SetBit(kUseBinPrior,false);
         break;
      case kBBayesian:
         SetBit(kIsBayesian,true);
         break;
      default:
         fStatisticOption = kFCP;
         fBoundary = &ClopperPearson;
         SetBit(kIsBayesian,false);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the title
///
/// Notes:
///       - The titles of the internal histograms are set to "title + (total)"
///         or "title + (passed)" respectively.
///       - It is possible to label the axis of the histograms as usual (see
///         TH1::SetTitle).
///
/// Example: Setting the title to "My Efficiency" and label the axis
///     pEff->SetTitle("My Efficiency;x label;eff");

void TEfficiency::SetTitle(const char* title)
{

   //setting the titles (looking for the first semicolon and insert the tokens there)
   TString title_passed = title;
   TString title_total = title;
   Ssiz_t pos = title_passed.First(";");
   if (pos != kNPOS) {
      title_passed.Insert(pos," (passed)");
      title_total.Insert(pos," (total)");
   }
   else {
      title_passed.Append(" (passed)");
      title_total.Append(" (total)");
   }
   fPassedHistogram->SetTitle(title_passed);
   fTotalHistogram->SetTitle(title_total);

   // strip (total) for the TEfficiency title
   // HIstogram SetTitle has already stripped the axis
   TString teffTitle = fTotalHistogram->GetTitle();
   teffTitle.ReplaceAll(" (total)","");
   TNamed::SetTitle(teffTitle);

}

////////////////////////////////////////////////////////////////////////////////
/// Sets the number of total events in the given global bin
///
/// returns "true" if the number of total events has been updated
/// otherwise "false" ist returned
///
/// Note: - requires: fPassedHistogram->GetBinContent(bin) <= events

Bool_t TEfficiency::SetTotalEvents(Int_t bin,Int_t events)
{
   if(events >= fPassedHistogram->GetBinContent(bin)) {
      fTotalHistogram->SetBinContent(bin,events);
      return true;
   }
   else {
      Error("SetTotalEvents(Int_t,Int_t)","passed number of events (%.1lf) in bin %i is bigger than given number of total events %i",fPassedHistogram->GetBinContent(bin),bin,events);
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the histogram containing all events
///
/// The given histogram is cloned and stored internally as histogram containing
/// all events. The given histogram has to be consistent with the current
/// fPassedHistogram (see CheckConsistency(const TH1&,const TH1&)).
/// The method returns whether the fTotalHistogram has been replaced (true) or
/// not (false).
///
/// Note: The list of associated functions fFunctions is cleared.
///
/// Option:
/// - "f": force the replacement without checking the consistency
///               This can lead to inconsistent histograms and useless results
///               or unexpected behaviour. But sometimes it might be the only
///               way to change the histograms. If you use this option, you
///               should ensure that the fPassedHistogram is replaced by a
///               consistent one (with respect to rTotal) as well.

Bool_t TEfficiency::SetTotalHistogram(const TH1& rTotal,Option_t* opt)
{
   TString option = opt;
   option.ToLower();

   Bool_t bReplace = option.Contains("f");

   if(!bReplace)
      bReplace = CheckConsistency(*fPassedHistogram,rTotal);

   if(bReplace) {
      delete fTotalHistogram;
      // do not add cloned histogram to gDirectory
      {
         TDirectory::TContext ctx(nullptr);
         fTotalHistogram = (TH1*)(rTotal.Clone());
      }
      fTotalHistogram->SetNormFactor(0);

      if(fFunctions)
         fFunctions->Delete();

      //check whether both histograms are filled with weights
      bool useWeights = CheckWeights(*fPassedHistogram,rTotal);
      SetUseWeightedEvents(useWeights);

      return true;
   }
   else
      return false;
}

////////////////////////////////////////////////////////////////////////////////

void TEfficiency::SetUseWeightedEvents(bool on)
{
   if (on && !TestBit(kUseWeights) )
       gROOT->Info("TEfficiency::SetUseWeightedEvents","Handle weighted events for computing efficiency");

   SetBit(kUseWeights,on);

   if (on && fTotalHistogram->GetSumw2N() != fTotalHistogram->GetNcells())
      fTotalHistogram->Sumw2();
   if (on && fPassedHistogram->GetSumw2N() != fTotalHistogram->GetNcells() )
      fPassedHistogram->Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the global weight for this TEfficiency object
///
/// Note: - weight has to be positive ( > 0)

void TEfficiency::SetWeight(Double_t weight)
{
   if(weight > 0)
      fWeight = weight;
   else
      Warning("SetWeight","invalid weight %.2lf",weight);
}

////////////////////////////////////////////////////////////////////////////////
/**
Calculates the boundaries for the frequentist Wilson interval

\param[in] total number of total events
\param[in] passed 0 <= number of passed events <= total
\param[in] level  confidence level
\param[in] bUpper
                - true  - upper boundary is returned
                - false - lower boundary is returned

Calculation:
\f{eqnarray*}{
      \alpha &=& 1 - \frac{level}{2}\\
      \kappa &=& \Phi^{-1}(1 - \alpha,1) ...\ normal\ quantile\ function\\
      mode &=& \frac{passed + \frac{\kappa^{2}}{2}}{total + \kappa^{2}}\\
      \Delta &=& \frac{\kappa}{total + \kappa^{2}} * \sqrt{passed (1 - \frac{passed}{total}) + \frac{\kappa^{2}}{4}}\\
      return &=& max(0,mode - \Delta)\ or\ min(1,mode + \Delta)
\f}

*/

Double_t TEfficiency::Wilson(Double_t total,Double_t passed,Double_t level,Bool_t bUpper)
{
   Double_t alpha = (1.0 - level)/2;
   if (total == 0) return (bUpper) ? 1 : 0;
   Double_t average = ((Double_t)passed) / total;
   Double_t kappa = ROOT::Math::normal_quantile(1 - alpha,1);

   Double_t mode = (passed + 0.5 * kappa * kappa) / (total + kappa * kappa);
   Double_t delta = kappa / (total + kappa*kappa) * std::sqrt(total * average
                                                              * (1 - average) + kappa * kappa / 4);
   if(bUpper)
      return ((mode + delta) > 1) ? 1.0 : (mode + delta);
   else
      return ((mode - delta) < 0) ? 0.0 : (mode - delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Addition operator
///
/// adds the corresponding histograms:
/// ~~~ {.cpp}
/// lhs.GetTotalHistogram() + rhs.GetTotalHistogram()
/// lhs.GetPassedHistogram() + rhs.GetPassedHistogram()
/// ~~~
/// the statistic option and the confidence level are taken from lhs

const TEfficiency operator+(const TEfficiency& lhs,const TEfficiency& rhs)
{
   TEfficiency tmp(lhs);
   tmp += rhs;
   return tmp;
}

#endif
