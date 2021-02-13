// @(#)root/splot:$Id$
// Author: Muriel Pivk, Anna Kreshuk    10/2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#include "TSPlot.h"
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TTreePlayer.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "TSelectorDraw.h"
#include "TBrowser.h"
#include "TClass.h"
#include "TMath.h"
#include "snprintf.h"

extern void Yields(Int_t &, Double_t *, Double_t &f, Double_t *x, Int_t iflag);

ClassImp(TSPlot);

/** \class TSPlot

A common method used in High Energy Physics to perform measurements is
the maximum Likelihood method, exploiting discriminating variables to
disentangle signal from background. The crucial point for such an
analysis to be reliable is to use an exhaustive list of sources of
events combined with an accurate description of all the Probability
Density Functions (PDF).

To assess the validity of the fit, a convincing quality check
is to explore further the data sample by examining the distributions of
control variables. A control variable can be obtained for instance by
removing one of the discriminating variables before performing again
the maximum Likelihood fit: this removed variable is a control
variable. The expected distribution of this control variable, for
signal, is to be compared to the one extracted, for signal, from the
data sample. In order to be able to do so, one must be able to unfold
from the distribution of the whole data sample.

The TSPlot method allows to reconstruct the distributions for
the control variable, independently for each of the various sources of
events, without making use of any <em>a priori</em> knowledge on <u>this</u>
variable. The aim is thus to use the knowledge available for the
discriminating variables to infer the behaviour of the individual
sources of events with respect to the control variable.

TSPlot is optimal if the control variable is uncorrelated with the discriminating variables.


A detail description of the formalism itself, called \f$\hbox{$_s$}{\cal P}lot\f$, is given
in [<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/node1.html#bib:sNIM">1</a>].

### The method


The \f$\hbox{$_s$}{\cal P}lot\f$ technique is developed in the above context of a
maximum Likelihood method making use of discriminating variables.

One considers a data sample in which are merged several species
of events. These species represent various signal components and
background components which all together account for the data sample.
The different terms of the log-Likelihood are:

  - \f$N\f$ : the total number of events in the data sample,
  - \f${\rm N}_{\rm s}\f$ : the number of species of events populating the data sample,
  - \f$N_i\f$ : the number of events expected on the average for the \f$i^{\rm th}\f$ species,
  - \f${\rm f}_i(y_e)\f$" : the value of the PDFs of the discriminating variables
    \f$y\f$" for the\f$i^{th}\f$ species and for event\f$e\f$",
  - \f$x\f$" : the set of control variables which, by definition, do not appear in
    the expression of the Likelihood function \f${\cal L}\f$.

The extended log-Likelihood reads:

 \f[
{\cal L}=\sum_{e=1}^{N}\ln \Big\{ \sum_{i=1}^{{\rm N}_{\rm s}}N_i{\rm f}_i(y_e) \Big\} -\sum_{i=1}^{{\rm N}_{\rm s}}N_i \tag{1}
\f]

From this expression, after maximization of \f${\cal L}\f$ with respect to the \f$N_i\f$ parameters,
a weight can be computed for every event and each species, in order to obtain later the true distribution
\f$\hbox{M}_i(x)\f$ of variable \f$x\f$. If \f${\rm n}\f$ is one of the
 \f${\rm N}_{\rm s}\f$ species present in the data sample, the weight for this species is defined by:


\f[
\fbox{$
{_s{\cal P}}_{\rm n}(y_e)={\sum_{j=1}^{{\rm N}_{\rm s}} \hbox{V}_{{\rm n}j}{\rm f}_j(y_e)\over\sum_{k=1}^{{\rm N}_{\rm s}}N_k{\rm f}_k(y_e) } $} , \tag{2}
\f]


where \f$\hbox{V}_{{\rm n}j}\f$

is the covariance matrix resulting from the Likelihood maximization.
This matrix can be used directly from the fit, but this is numerically
less accurate than the direct computation:


\f[
\hbox{ V}^{-1}_{{\rm n}j}~=~
{\partial^2(-{\cal L})\over\partial N_{\rm n}\partial N_j}~=~
\sum_{e=1}^N {{\rm f}_{\rm n}(y_e){\rm f}_j(y_e)\over(\sum_{k=1}^{{\rm N}_{\rm s}}N_k{\rm f}_k(y_e))^2} . \tag{3}
\f]


The distribution of the control variable \f$x\f$ obtained by histogramming the weighted
events reproduces, on average, the true distribution
\f${\hbox{ {M}}}_{\rm n}(x)\f$

The class TSPlot allows to reconstruct the true distribution
\f${\hbox{ {M}}}_{\rm n}(x)\f$

of a control variable \f$x\f$ for each of the \f${\rm N}_{\rm s}\f$ species from
the sole knowledge of the PDFs of the discriminating variables \f${\rm f}_i(y)\f$.
The plots obtained thanks to the TSPlot class are called \f$\hbox {$_s$}{\cal P}lots\f$.


### Some properties and checks


Beside reproducing the true distribution,\f$\hbox {$_s$}{\cal P}lots\f$ bear remarkable properties:


  - Each \f$x\f$ - distribution is properly normalized:

\f[
\sum_{e=1}^{N} {_s{\cal P}}_{\rm n}(y_e)~=~N_{\rm n} ~. \tag{4}
\f]


  - For any event:

\f[
\sum_{l=1}^{{\rm N}_{\rm s}} {_s{\cal P}}_l(y_e) ~=~1 ~. \tag{5}
\f]

    That is to say that, summing up the \f${\rm N}_{\rm s}\f$ \f$\hbox {$_s$}{\cal P}lots\f$,
    one recovers the data sample distribution in \f$x\f$, and summing up the number of events
    entering in a \f$\hbox{$_s$}{\cal P}lot\f$ for a given species, one recovers the yield of the
    species, as provided by the fit.
    The property <a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#eq:NormalizationOK">4</a> is implemented in the TSPlot class as a check.

  - the sum of the statistical uncertainties per bin

\f[
\sigma[N_{\rm n}\  _s\tilde{\rm M}_{\rm n}(x) {\delta x}]~=~\sqrt{\sum_{e \subset {\delta x}} ({_s{\cal P}}_{\rm n})^2} ~. \tag{6}
\f]

    reproduces the statistical uncertainty on the yield \f$N_{\rm n}\f$, as provided by the fit:
    \f$\sigma[N_{\rm n}]\equiv\sqrt{\hbox{ V}_{{\rm n}{\rm n}}}\f$ .
    Because of that and since the determination of the yields is optimal
    when obtained using a Likelihood fit, one can conclude that the \f$\hbox{$_s$}{\cal P}lot\f$
    technique is itself an optimal method to reconstruct distributions of control variables.


### Different steps followed by TSPlot


  1. A maximum Likelihood fit is performed to obtain the yields \f$N_i\f$
     of the various species.The fit relies on discriminating variables \f$y\f$
     uncorrelated with a control variable \f$x\f$:
     the later is therefore totally absent from the fit.

  2. The weights \f${_s{\cal P}}\f$ are calculated using Eq.
     (<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#eq:weightxnotiny">2</a>)
     where the covariance matrix is taken from Minuit.

  3. Histograms of \f$x\f$ are filled by weighting the events with \f${_s{\cal P}}\f$ .

  4. Error bars per bin are given by Eq. (<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#eq:ErrorPerBin">6</a>).


The \f$\hbox {$_s$}{\cal P}lots\f$ reproduce the true distributions of the species
in the control variable \f$x\f$, within the above defined statistical uncertainties.

### Illustrations


To illustrate the technique, one considers an example derived from the analysis where
\f$\hbox {$_s$}{\cal P}lots\f$
have been first used (charmless B decays). One is dealing with a data
sample in which two species are present: the first is termed signal and
the second background. A maximum Likelihood fit is performed to obtain
the two yields \f$N_1\f$ and \f$N_2\f$ . The fit relies on two discriminating
variables collectively denoted \f$y\f$ which are chosen within three possible
variables denoted \f${m_{\rm ES}}\f$ , \f$\Delta E\f$ and \f${\cal F}\f$.
The variable which is not incorporated in \f$y\f$ is used as the control variable
\f$x\f$ . The six distributions of the three variables are assumed to be the ones
depicted in Fig. <a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:pdfs">1</a>.


\image html splot_pdfmesNIM.png width=800

#### Figure 1:

Distributions of the three discriminating variables available to perform the Likelihood fit:
\f${m_{\rm ES}}\f$ , \f$\Delta E\f$ , \f${\cal F}\f$ .
Among the three variables, two are used to perform the fit while one is
kept out of the fit to serve the purpose of a control variable. The
three distributions on the top (resp. bottom) of the figure correspond
to the signal (resp. background). The unit of the vertical axis is
chosen such that it indicates the number of entries per bin, if one
slices the histograms in 25 bins.

A data sample being built through a Monte Carlo simulation based on the
distributions shown in Fig.
<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:pdfs">1</a>,
one obtains the three distributions of Fig.
<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:pdfstot">2</a>.
Whereas the distribution of \f$\Delta E\f$  clearly indicates the presence of the signal,
the distribution of \f${m_{\rm ES}}\f$ and \f${\cal F}\f$  are less obviously populated by signal.


\image html splot_genfiTOTNIM.png  width=800

#### Figure 2:

Distributions of the three discriminating variables for signal plus
background. The three distributions are the ones obtained from a data
sample obtained through a Monte Carlo simulation based on the
distributions shown in Fig.
<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:pdfs">1</a>.
The data sample consists of 500 signal events and 5000 background events.


Choosing \f$\Delta E\f$ and \f${\cal F}\f$ as discriminating variables to determine
\f$N_1\f$ and \f$N_2\f$ through a maximum Likelihood fit, one builds, for the control
variable \f${m_{\rm ES}}\f$  which is unknown to the fit, the two \f$\hbox {$_s$}{\cal P}lots\f$
for signal and background shown in
Fig. <a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:messPlots">3</a>.
One observes that the \f$\hbox{$_s$}{\cal P}lot\f$
for signal reproduces correctly the PDF even where the latter vanishes,
although the error bars remain sizeable. This results from the almost
complete cancellation between positive and negative weights: the sum of
weights is close to zero while the sum of weights squared is not. The
occurence of negative weights occurs through the appearance of the
covariance matrix, and its negative components, in the definition of
Eq. (<a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#eq:weightxnotiny">2</a>).


A word of caution is in order with respect to the error bars. Whereas
their sum in quadrature is identical to the statistical uncertainties
of the yields determined by the fit, and if, in addition, they are
asymptotically correct, the error bars should be handled with care for
low statistics and/or for too fine binning. This is because the error
bars do not incorporate two known properties of the PDFs: PDFs are
positive definite and can be non-zero in a given x-bin, even if in the
particular data sample at hand, no event is observed in this bin. The
latter limitation is not specific to \f$\hbox {$_s$}{\cal P}lots\f$ ,
rather it is always present when one is willing to infer the PDF at the
origin of an histogram, when, for some bins, the number of entries does
not guaranty the applicability of the Gaussian regime. In such
situations, a satisfactory practice is to attach allowed ranges to the
histogram to indicate the upper and lower limits of the PDF value which
are consistent with the actual observation, at a given confidence
level.


\image html splot_mass-bkg-sPlot.png  width=600

#### Figure 3:

The \f$\hbox {$_s$}{\cal P}lots\f$ (signal on top, background on bottom)
obtained for \f${m_{\rm ES}}\f$ are represented as dots with error bars.
They are obtained from a fit using only information from \f$\Delta E\f$ and
\f${\cal F}\f$

<p>
Choosing \f${m_{\rm ES}}\f$ and \f$\Delta E\f$ as discriminating variables to
determine \f$N_1\f$ and \f$N_2\f$ through a maximum Likelihood fit, one builds,
for the control variable \f${\cal F}\f$ which is unknown to the fit, the two
\f$\hbox {$_s$}{\cal P}lots\f$ for signal and background shown in
Fig. <a href="http://www.slac.stanford.edu/%7Epivk/sPlot/sPlot_ROOT/sPlot_ROOT.html#fig:FisPlots">4</a>.
In the \f$\hbox{$_s$}{\cal P}lot\f$ for signal one observes that error bars are
the largest in the \f$x\f$ regions where the background is the largest.


\image html splot_fisher-bkg-sPlot.png width=600

#### Figure 4:

The \f$\hbox {$_s$}{\cal P}lots\f$ (signal on top, background on bottom) obtained
for \f${\cal F}\f$ are represented as dots with error bars. They are obtained
from a fit using only information from \f${m_{\rm ES}}\f$ and \f$\Delta E\f$

The results above can be obtained by running the tutorial TestSPlot.C
*/


////////////////////////////////////////////////////////////////////////////////
/// default constructor (used by I/O only)

TSPlot::TSPlot() :
 fTree(0),
 fTreename(0),
 fVarexp(0),
 fSelection(0)
{
   fNx = 0;
   fNy=0;
   fNevents = 0;
   fNSpecies=0;
   fNumbersOfEvents=0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal TSPlot constructor
///  - nx :  number of control variables
///  - ny :  number of discriminating variables
///  - ne :  total number of events
///  - ns :  number of species
///  - tree: input data

TSPlot::TSPlot(Int_t nx, Int_t ny, Int_t ne, Int_t ns, TTree *tree) :
 fTreename(0),
 fVarexp(0),
 fSelection(0)

{
   fNx = nx;
   fNy=ny;
   fNevents = ne;
   fNSpecies=ns;

   fXvar.ResizeTo(fNevents, fNx);
   fYvar.ResizeTo(fNevents, fNy);
   fYpdf.ResizeTo(fNevents, fNSpecies*fNy);
   fSWeights.ResizeTo(fNevents, fNSpecies*(fNy+1));
   fTree = tree;
   fNumbersOfEvents = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TSPlot::~TSPlot()
{
   if (fNumbersOfEvents)
      delete [] fNumbersOfEvents;
   if (!fXvarHists.IsEmpty())
      fXvarHists.Delete();
   if (!fYvarHists.IsEmpty())
      fYvarHists.Delete();
   if (!fYpdfHists.IsEmpty())
      fYpdfHists.Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// To browse the histograms

void TSPlot::Browse(TBrowser *b)
{
   if (!fSWeightsHists.IsEmpty()) {
      TIter next(&fSWeightsHists);
      TH1D* h = 0;
      while ((h = (TH1D*)next()))
         b->Add(h,h->GetName());
   }

   if (!fYpdfHists.IsEmpty()) {
      TIter next(&fYpdfHists);
      TH1D* h = 0;
      while ((h = (TH1D*)next()))
         b->Add(h,h->GetName());
   }
   if (!fYvarHists.IsEmpty()) {
      TIter next(&fYvarHists);
      TH1D* h = 0;
      while ((h = (TH1D*)next()))
         b->Add(h,h->GetName());
   }
   if (!fXvarHists.IsEmpty()) {
      TIter next(&fXvarHists);
      TH1D* h = 0;
      while ((h = (TH1D*)next()))
         b->Add(h,h->GetName());
   }
   b->Add(&fSWeights, "sWeights");
}

////////////////////////////////////////////////////////////////////////////////
/// Set the initial number of events of each species - used
/// as initial estimates in minuit

void TSPlot::SetInitialNumbersOfSpecies(Int_t *numbers)
{
   if (!fNumbersOfEvents)
      fNumbersOfEvents = new Double_t[fNSpecies];
   for (Int_t i=0; i<fNSpecies; i++)
      fNumbersOfEvents[i]=numbers[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the sWeights
///
/// The option controls the print level
///  - "Q" - no print out
///  - "V" - prints the estimated #of events in species - default
///  - "VV" - as "V" + the minuit printing + sums of weights for control

void TSPlot::MakeSPlot(Option_t *option)
{

   if (!fNumbersOfEvents){
      Error("MakeSPlot","Initial numbers of events in species have not been set");
      return;
   }
   Int_t i, j, ispecies;

   TString opt = option;
   opt.ToUpper();
   opt.ReplaceAll("VV", "W");

   //make sure that global fitter is minuit
   char s[]="TFitter";
   if (TVirtualFitter::GetFitter()){
      Int_t strdiff=strcmp(TVirtualFitter::GetFitter()->IsA()->GetName(), s);
      if (strdiff!=0)
         delete TVirtualFitter::GetFitter();
   }


   TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 2);
   fPdfTot.ResizeTo(fNevents, fNSpecies);

   //now let's do it, excluding different yvars
   //for iplot = -1 none is excluded
   for (Int_t iplot=-1; iplot<fNy; iplot++){
      for (i=0; i<fNevents; i++){
         for (ispecies=0; ispecies<fNSpecies; ispecies++){
            fPdfTot(i, ispecies)=1;
            for (j=0; j<fNy; j++){
               if (j!=iplot)
                  fPdfTot(i, ispecies)*=fYpdf(i, ispecies*fNy+j);
            }
         }
      }
      minuit->Clear();
      minuit->SetFCN(Yields);
      Double_t arglist[10];
      //set the print level
      if (opt.Contains("Q")||opt.Contains("V")){
         arglist[0]=-1;
      }
      if (opt.Contains("W"))
         arglist[0]=0;
      minuit->ExecuteCommand("SET PRINT", arglist, 1);

      minuit->SetObjectFit(&fPdfTot); //a tricky way to get fPdfTot matrix to fcn
      for (ispecies=0; ispecies<fNSpecies; ispecies++)
         minuit->SetParameter(ispecies, "", fNumbersOfEvents[ispecies], 1, 0, 0);

      minuit->ExecuteCommand("MIGRAD", arglist, 0);
      for (ispecies=0; ispecies<fNSpecies; ispecies++){
         fNumbersOfEvents[ispecies]=minuit->GetParameter(ispecies);
         if (!opt.Contains("Q"))
            printf("estimated #of events in species %d = %f\n", ispecies, fNumbersOfEvents[ispecies]);
      }
      if (!opt.Contains("Q"))
         printf("\n");
      Double_t *covmat = minuit->GetCovarianceMatrix();
      SPlots(covmat, iplot);

      if (opt.Contains("W")){
         Double_t *sumweight = new Double_t[fNSpecies];
         for (i=0; i<fNSpecies; i++){
            sumweight[i]=0;
            for (j=0; j<fNevents; j++)
               sumweight[i]+=fSWeights(j, (iplot+1)*fNSpecies + i);
            printf("checking sum of weights[%d]=%f\n", i, sumweight[i]);
         }
         printf("\n");
         delete [] sumweight;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the sWeights from the covariance matrix

void TSPlot::SPlots(Double_t *covmat, Int_t i_excl)
{
   Int_t i, ispecies, k;
   Double_t numerator, denominator;
   for (i=0; i<fNevents; i++){
      denominator=0;
      for (ispecies=0; ispecies<fNSpecies; ispecies++)
         denominator+=fNumbersOfEvents[ispecies]*fPdfTot(i, ispecies);
      for (ispecies=0; ispecies<fNSpecies; ispecies++){
         numerator=0;
         for (k=0; k<fNSpecies; k++)
            numerator+=covmat[ispecies*fNSpecies+k]*fPdfTot(i, k);
         fSWeights(i, (i_excl+1)*fNSpecies + ispecies)=numerator/denominator;
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Returns the matrix of sweights

void TSPlot::GetSWeights(TMatrixD &weights)
{
   if (weights.GetNcols()!=fNSpecies*(fNy+1) || weights.GetNrows()!=fNevents)
      weights.ResizeTo(fNevents, fNSpecies*(fNy+1));
   weights = fSWeights;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the matrix of sweights. It is assumed that the array passed in the
/// argurment is big enough

void TSPlot::GetSWeights(Double_t *weights)
{
   for (Int_t i=0; i<fNevents; i++){
      for (Int_t j=0; j<fNSpecies; j++){
         weights[i*fNSpecies+j]=fSWeights(i, j);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the histograms of x variables (not weighted) with nbins

void TSPlot::FillXvarHists(Int_t nbins)
{
   Int_t i, j;

   if (!fXvarHists.IsEmpty()){
      if (((TH1D*)fXvarHists.First())->GetNbinsX()!=nbins)
         fXvarHists.Delete();
      else
         return;
   }

   //make the histograms
   char name[12];
   for (i=0; i<fNx; i++){
      snprintf(name,sizeof(name), "x%d", i);
      TH1D *h = new TH1D(name, name, nbins, fMinmax(0, i), fMinmax(1, i));
      for (j=0; j<fNevents; j++)
         h->Fill(fXvar(j, i));
      fXvarHists.Add(h);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Returns the array of histograms of x variables (not weighted).
/// If histograms have not already
/// been filled, they are filled with default binning 100.

TObjArray* TSPlot::GetXvarHists()
{
   Int_t nbins = 100;
   if (fXvarHists.IsEmpty())
      FillXvarHists(nbins);
   else if (((TH1D*)fXvarHists.First())->GetNbinsX()!=nbins)
      FillXvarHists(nbins);
   return &fXvarHists;
}

////////////////////////////////////////////////////////////////////////////////
///Returns the histogram of variable ixvar.
/// If histograms have not already
/// been filled, they are filled with default binning 100.

TH1D *TSPlot::GetXvarHist(Int_t ixvar)
{
   Int_t nbins = 100;
   if (fXvarHists.IsEmpty())
      FillXvarHists(nbins);
   else if (((TH1D*)fXvarHists.First())->GetNbinsX()!=nbins)
      FillXvarHists(nbins);

   return (TH1D*)fXvarHists.UncheckedAt(ixvar);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the histograms of y variables

void TSPlot::FillYvarHists(Int_t nbins)
{
   Int_t i, j;

   if (!fYvarHists.IsEmpty()){
      if (((TH1D*)fYvarHists.First())->GetNbinsX()!=nbins)
         fYvarHists.Delete();
      else
         return;
   }

   //make the histograms
   char name[12];
   for (i=0; i<fNy; i++){
      snprintf(name,sizeof(name), "y%d", i);
      TH1D *h=new TH1D(name, name, nbins, fMinmax(0, fNx+i), fMinmax(1, fNx+i));
      for (j=0; j<fNevents; j++)
         h->Fill(fYvar(j, i));
      fYvarHists.Add(h);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the array of histograms of y variables. If histograms have not already
/// been filled, they are filled with default binning 100.

TObjArray* TSPlot::GetYvarHists()
{
   Int_t nbins = 100;
   if (fYvarHists.IsEmpty())
      FillYvarHists(nbins);
   else if (((TH1D*)fYvarHists.First())->GetNbinsX()!=nbins)
      FillYvarHists(nbins);
   return &fYvarHists;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the histogram of variable iyvar.If histograms have not already
/// been filled, they are filled with default binning 100.

TH1D *TSPlot::GetYvarHist(Int_t iyvar)
{
   Int_t nbins = 100;
   if (fYvarHists.IsEmpty())
      FillYvarHists(nbins);
   else if (((TH1D*)fYvarHists.First())->GetNbinsX()!=nbins)
      FillYvarHists(nbins);
   return (TH1D*)fYvarHists.UncheckedAt(iyvar);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the histograms of pdf-s of y variables with binning nbins

void TSPlot::FillYpdfHists(Int_t nbins)
{
   Int_t i, j, ispecies;

   if (!fYpdfHists.IsEmpty()){
      if (((TH1D*)fYpdfHists.First())->GetNbinsX()!=nbins)
         fYpdfHists.Delete();
      else
         return;
   }

   char name[34];
   for (ispecies=0; ispecies<fNSpecies; ispecies++){
      for (i=0; i<fNy; i++){
         snprintf(name,sizeof(name), "pdf_species%d_y%d", ispecies, i);
         //TH1D *h = new TH1D(name, name, nbins, ypdfmin[ispecies*fNy+i], ypdfmax[ispecies*fNy+i]);
         TH1D *h = new TH1D(name, name, nbins, fMinmax(0, fNx+fNy+ispecies*fNy+i), fMinmax(1, fNx+fNy+ispecies*fNy+i));
         for (j=0; j<fNevents; j++)
            h->Fill(fYpdf(j, ispecies*fNy+i));
         fYpdfHists.Add(h);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the array of histograms of pdf's of y variables with binning nbins.
/// If histograms have not already
/// been filled, they are filled with default binning 100.

TObjArray* TSPlot::GetYpdfHists()
{
   Int_t nbins = 100;
   if (fYpdfHists.IsEmpty())
      FillYpdfHists(nbins);

   return &fYpdfHists;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the histogram of the pdf of variable iyvar for species #ispecies, binning nbins.
/// If histograms have not already
/// been filled, they are filled with default binning 100.

TH1D *TSPlot::GetYpdfHist(Int_t iyvar, Int_t ispecies)
{
   Int_t nbins = 100;
   if (fYpdfHists.IsEmpty())
      FillYpdfHists(nbins);

   return (TH1D*)fYpdfHists.UncheckedAt(fNy*ispecies+iyvar);
}

////////////////////////////////////////////////////////////////////////////////
/// The order of histograms in the array:
///
/// x0_species0, x0_species1,..., x1_species0, x1_species1,..., y0_species0, y0_species1,...
///
/// If the histograms have already been filled with a different binning, they are refilled
/// and all histograms are deleted

void TSPlot::FillSWeightsHists(Int_t nbins)
{
   if (fSWeights.GetNoElements()==0){
      Error("GetSWeightsHists", "SWeights were not computed");
      return;
   }

   if (!fSWeightsHists.IsEmpty()){
      if (((TH1D*)fSWeightsHists.First())->GetNbinsX()!=nbins)
         fSWeightsHists.Delete();
      else
         return;
   }

   char name[30];

   //Fill histograms of x-variables weighted with sWeights
   for (Int_t ivar=0; ivar<fNx; ivar++){
      for (Int_t ispecies=0; ispecies<fNSpecies; ispecies++){
            snprintf(name,30, "x%d_species%d", ivar, ispecies);
            TH1D *h = new TH1D(name, name, nbins, fMinmax(0, ivar), fMinmax(1, ivar));
            h->Sumw2();
            for (Int_t ievent=0; ievent<fNevents; ievent++)
               h->Fill(fXvar(ievent, ivar), fSWeights(ievent, ispecies));
            fSWeightsHists.AddLast(h);
         }
   }

   //Fill histograms of y-variables (excluded from the fit), weighted with sWeights
   for (Int_t iexcl=0; iexcl<fNy; iexcl++){
      for(Int_t ispecies=0; ispecies<fNSpecies; ispecies++){
            snprintf(name,30, "y%d_species%d", iexcl, ispecies);
            TH1D *h = new TH1D(name, name, nbins, fMinmax(0, fNx+iexcl), fMinmax(1, fNx+iexcl));
            h->Sumw2();
            for (Int_t ievent=0; ievent<fNevents; ievent++)
               h->Fill(fYvar(ievent, iexcl), fSWeights(ievent, fNSpecies*(iexcl+1)+ispecies));
            fSWeightsHists.AddLast(h);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an array of all histograms of variables, weighted with sWeights.
/// If histograms have not been already filled, they are filled with default binning 50
/// The order of histograms in the array:
///
/// x0_species0, x0_species1,..., x1_species0, x1_species1,..., y0_species0, y0_species1,...

TObjArray *TSPlot::GetSWeightsHists()
{
   Int_t nbins = 50; //default binning
   if (fSWeightsHists.IsEmpty())
      FillSWeightsHists(nbins);

   return &fSWeightsHists;
}

////////////////////////////////////////////////////////////////////////////////
/// The Fill...Hist() methods fill the histograms with the real limits on the variables
/// This method allows to refill the specified histogram with user-set boundaries min and max
///
///Parameters:
///
///  - type = 1 - histogram of x variable #nvar
///  - type = 2 - histogram of y variable #nvar
///  - type = 3 - histogram of y_pdf for y #nvar and species #nspecies
///  - type = 4 - histogram of x variable #nvar, species #nspecies, WITH sWeights
///  - type = 5 - histogram of y variable #nvar, species #nspecies, WITH sWeights

void TSPlot::RefillHist(Int_t type, Int_t nvar, Int_t nbins, Double_t min, Double_t max, Int_t nspecies)
{
   if (type<1 || type>5){
      Error("RefillHist", "type must lie between 1 and 5");
      return;
   }
   char name[20];
   Int_t j;
   TH1D *hremove;
   if (type==1){
      hremove = (TH1D*)fXvarHists.RemoveAt(nvar);
      delete hremove;
      snprintf(name,20,"x%d",nvar);
      TH1D *h = new TH1D(name, name, nbins, min, max);
      for (j=0; j<fNevents;j++)
         h->Fill(fXvar(j, nvar));
      fXvarHists.AddAt(h, nvar);
   }
   if (type==2){
      hremove = (TH1D*)fYvarHists.RemoveAt(nvar);
      delete hremove;
      snprintf(name,20, "y%d", nvar);
      TH1D *h = new TH1D(name, name, nbins, min, max);
      for (j=0; j<fNevents;j++)
         h->Fill(fYvar(j, nvar));
      fXvarHists.AddAt(h, nvar);
   }
   if (type==3){
      hremove = (TH1D*)fYpdfHists.RemoveAt(nspecies*fNy+nvar);
      delete hremove;
      snprintf(name,20, "pdf_species%d_y%d", nspecies, nvar);
      TH1D *h=new TH1D(name, name, nbins, min, max);
      for (j=0; j<fNevents; j++)
         h->Fill(fYpdf(j, nspecies*fNy+nvar));
      fYpdfHists.AddAt(h, nspecies*fNy+nvar);
   }
   if (type==4){
      hremove = (TH1D*)fSWeightsHists.RemoveAt(fNSpecies*nvar+nspecies);
      delete hremove;
      snprintf(name,20, "x%d_species%d", nvar, nspecies);
      TH1D *h = new TH1D(name, name, nbins, min, max);
      h->Sumw2();
      for (Int_t ievent=0; ievent<fNevents; ievent++)
         h->Fill(fXvar(ievent, nvar), fSWeights(ievent, nspecies));
      fSWeightsHists.AddAt(h, fNSpecies*nvar+nspecies);
   }
   if (type==5){
      hremove = (TH1D*)fSWeightsHists.RemoveAt(fNx*fNSpecies + fNSpecies*nvar+nspecies);
      delete hremove;
      snprintf(name,20, "y%d_species%d", nvar, nspecies);
      TH1D *h = new TH1D(name, name, nbins, min, max);
      h->Sumw2();
      for (Int_t ievent=0; ievent<fNevents; ievent++)
         h->Fill(fYvar(ievent, nvar), fSWeights(ievent, nspecies));
      fSWeightsHists.AddAt(h, fNx*fNSpecies + fNSpecies*nvar+nspecies);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Returns the histogram of a variable, weighted with sWeights.
///  - If histograms have not been already filled, they are filled with default binning 50
///  - If parameter ixvar!=-1, the histogram of x-variable ixvar is returned for species ispecies
///  - If parameter ixvar==-1, the histogram of y-variable iyexcl is returned for species ispecies
///  - If the histogram has already been filled and the binning is different from the parameter nbins
///    all histograms with old binning will be deleted and refilled.

TH1D *TSPlot::GetSWeightsHist(Int_t ixvar, Int_t ispecies,Int_t iyexcl)
{

   Int_t nbins = 50; //default binning
   if (fSWeightsHists.IsEmpty())
      FillSWeightsHists(nbins);

   if (ixvar==-1)
      return (TH1D*)fSWeightsHists.UncheckedAt(fNx*fNSpecies + fNSpecies*iyexcl+ispecies);
   else
      return (TH1D*)fSWeightsHists.UncheckedAt(fNSpecies*ixvar + ispecies);

}


////////////////////////////////////////////////////////////////////////////////
/// Set the input Tree

void TSPlot::SetTree(TTree *tree)
{
   fTree = tree;
}

////////////////////////////////////////////////////////////////////////////////
///Specifies the variables from the tree to be used for splot
///
///Variables fNx, fNy, fNSpecies and fNEvents should already be set!
///
///In the 1st parameter it is assumed that first fNx variables are x(control variables),
///then fNy y variables (discriminating variables),
///then fNy*fNSpecies ypdf variables (probability distribution functions of discriminating
///variables for different species). The order of pdfs should be: species0_y0, species0_y1,...
///species1_y0, species1_y1,...species[fNSpecies-1]_y0...
///The 2nd parameter allows to make a cut
///TTree::Draw method description contains more details on specifying expression and selection

void TSPlot::SetTreeSelection(const char* varexp, const char *selection, Long64_t firstentry)
{
   TTreeFormula **var;
   std::vector<TString> cnames;
   TList formulaList;
   TSelectorDraw *selector = (TSelectorDraw*)(((TTreePlayer*)fTree->GetPlayer())->GetSelector());

   Long64_t entry, entryNumber;
   Int_t i,nch;
   Int_t ncols;
   TObjArray *leaves = fTree->GetListOfLeaves();

   fTreename= new TString(fTree->GetName());
   if (varexp)
      fVarexp = new TString(varexp);
   if (selection)
      fSelection= new TString(selection);

   nch = varexp ? strlen(varexp) : 0;


//*-*- Compile selection expression if there is one
   TTreeFormula *select = nullptr;
   if (selection && strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select) return;
      if (!select->GetNdim()) { delete select; return; }
      formulaList.Add(select);
   }
//*-*- if varexp is empty, take first nx + ny + ny*nspecies columns by default

   if (nch == 0) {
      ncols = fNx + fNy + fNy*fNSpecies;
      for (i=0;i<ncols;i++) {
         cnames.push_back( leaves->At(i)->GetName() );
      }
//*-*- otherwise select only the specified columns
   } else {
      ncols = selector->SplitNames(varexp,cnames);
   }
   var = new TTreeFormula* [ncols];
   Double_t *xvars = new Double_t[ncols];

   fMinmax.ResizeTo(2, ncols);
   for (i=0; i<ncols; i++){
      fMinmax(0, i)=1e30;
      fMinmax(1, i)=-1e30;
   }

//*-*- Create the TreeFormula objects corresponding to each column
   for (i=0;i<ncols;i++) {
      var[i] = new TTreeFormula("Var1",cnames[i].Data(),fTree);
      formulaList.Add(var[i]);
   }

//*-*- Create a TreeFormulaManager to coordinate the formulas
   TTreeFormulaManager *manager = nullptr;
   if (formulaList.LastIndex()>=0) {
      manager = new TTreeFormulaManager;
      for(i=0;i<=formulaList.LastIndex();i++) {
         manager->Add((TTreeFormula*)formulaList.At(i));
      }
      manager->Sync();
      formulaList.Clear();
   }
//*-*- loop on all selected entries
   // fSelectedRows = 0;
   Int_t tnumber = -1;
   Long64_t selectedrows=0;
   for (entry=firstentry;entry<firstentry+fNevents;entry++) {
      entryNumber = fTree->GetEntryNumber(entry);
      if (entryNumber < 0) break;
      Long64_t localEntry = fTree->LoadTree(entryNumber);
      if (localEntry < 0) break;
      if (tnumber != fTree->GetTreeNumber()) {
         tnumber = fTree->GetTreeNumber();
         if (manager) manager->UpdateFormulaLeaves();
      }
      Int_t ndata = 1;
      if (manager && manager->GetMultiplicity()) {
         ndata = manager->GetNdata();
      }

      for(Int_t inst=0;inst<ndata;inst++) {
         Bool_t loaded = kFALSE;
         if (select) {
            if (select->EvalInstance(inst) == 0) {
               continue;
            }
         }

         if (inst==0) loaded = kTRUE;
         else if (!loaded) {
            // EvalInstance(0) always needs to be called so that
            // the proper branches are loaded.
            for (i=0;i<ncols;i++) {
               var[i]->EvalInstance(0);
            }
            loaded = kTRUE;
         }

         for (i=0;i<ncols;i++) {
            xvars[i] = var[i]->EvalInstance(inst);
         }

         // curentry = entry-firstentry;
         //printf("event#%d\n", curentry);
         //for (i=0; i<ncols; i++)
          //  printf("xvars[%d]=%f\n", i, xvars[i]);
         //selectedrows++;
         for (i=0; i<fNx; i++){
            fXvar(selectedrows, i) = xvars[i];
            if (fXvar(selectedrows, i) < fMinmax(0, i))
               fMinmax(0, i)=fXvar(selectedrows, i);
            if (fXvar(selectedrows, i) > fMinmax(1, i))
               fMinmax(1, i)=fXvar(selectedrows, i);
         }
         for (i=0; i<fNy; i++){
            fYvar(selectedrows, i) = xvars[i+fNx];
            //printf("y_in_loop(%d, %d)=%f, xvars[%d]=%f\n", selectedrows, i, fYvar(selectedrows, i), i+fNx, xvars[i+fNx]);
            if (fYvar(selectedrows, i) < fMinmax(0, i+fNx))
               fMinmax(0, i+fNx) = fYvar(selectedrows, i);
            if (fYvar(selectedrows, i) > fMinmax(1, i+fNx))
               fMinmax(1, i+fNx) = fYvar(selectedrows, i);
            for (Int_t j=0; j<fNSpecies; j++){
               fYpdf(selectedrows, j*fNy + i)=xvars[j*fNy + i+fNx+fNy];
               if (fYpdf(selectedrows, j*fNy+i) < fMinmax(0, j*fNy+i+fNx+fNy))
                  fMinmax(0, j*fNy+i+fNx+fNy) = fYpdf(selectedrows, j*fNy+i);
               if (fYpdf(selectedrows, j*fNy+i) > fMinmax(1, j*fNy+i+fNx+fNy))
                  fMinmax(1, j*fNy+i+fNx+fNy) = fYpdf(selectedrows, j*fNy+i);
            }
         }
      selectedrows++;
      }
   }
   fNevents=selectedrows;
  // for (i=0; i<fNevents; i++){
    //  printf("event#%d\n", i);
      //for (Int_t iy=0; iy<fNy; iy++)
        // printf("y[%d]=%f\n", iy, fYvar(i, iy));
      //for (Int_t ispecies=0; ispecies<fNSpecies; ispecies++){
      //   for (Int_t iy=0; iy<fNy; iy++)
        //    printf("ypdf[sp. %d, y %d]=%f\n", ispecies, iy, fYpdf(i, ispecies*fNy+iy));
     // }
   //}

   delete [] xvars;
   delete [] var;
}

////////////////////////////////////////////////////////////////////////////////
/// FCN-function for Minuit

void Yields(Int_t &, Double_t *, Double_t &f, Double_t *x, Int_t /*iflag*/)
{
   Double_t lik;
   Int_t i, ispecies;

   TVirtualFitter *fitter = TVirtualFitter::GetFitter();
   TMatrixD *pdftot = (TMatrixD*)fitter->GetObjectFit();
   Int_t nev = pdftot->GetNrows();
   Int_t nes = pdftot->GetNcols();
   f=0;
   for (i=0; i<nev; i++){
      lik=0;
      for (ispecies=0; ispecies<nes; ispecies++)
         lik+=x[ispecies]*(*pdftot)(i, ispecies);
      if (lik<0) lik=1;
      f+=TMath::Log(lik);
   }
   //extended likelihood, equivalent to chi2
   Double_t ntot=0;
   for (i=0; i<nes; i++)
      ntot += x[i];
   f = -2*(f-ntot);
}

