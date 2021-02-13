// @(#)root/unfold:$Id$
// Authors: Stefan Schmitt, Amnon Harel DESY and CERN, 11/08/11

/** \class TUnfoldDensity
\ingroup Unfold
An algorithm to unfold distributions from detector to truth level

TUnfoldDensity is used to decompose a measurement y into several sources x,
given the measurement uncertainties, background b and a matrix of migrations A.
The method can be applied to a large number of problems,
where the measured distribution y is a linear superposition
of several Monte Carlo shapes. Beyond such a simple template fit,
TUnfoldDensity has an adjustable regularisation term and also supports an
optional constraint on the total number of events.
Background sources can be specified, with a normalisation constant and
normalisation uncertainty. In addition, variants of the response
matrix may be specified, these are taken to determine systematic
uncertainties. Complex, multidimensional arrangements of signal and
background bins are managed with the help of the class TUnfoldBinning.

If you use this software, please consider the following citation

<b>S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]</b>

Detailed documentation and updates are available on
http://www.desy.de/~sschmitt

### Brief recipe to use TUnfoldSys:

  - Set up binning schemes for the truth and measured
distributions. The binning schemes may be coded in the XML language,
for reading use TUnfoldBinningXML.
  - A matrix (truth,reconstructed) is given as a two-dimensional histogram
    as argument to the constructor of TUnfold
  - A vector of measurements is given as one-dimensional histogram using
    the SetInput() method
  - Repeated calls to SubtractBackground() to specify background sources
  - Repeated calls to AddSysError() to specify systematic uncertainties
    - The unfolding is performed

    - either once with a fixed parameter tau, method DoUnfold(tau)
    - or multiple times in a scan to determine the best choice of tau,
      method ScanLCurve()
    - or multiple times in a scan to determine the best choice of tau,
      method ScanTau()

  - Unfolding results are retrieved using various GetXXX() methods

A detailed documentation of the various GetXXX() methods to control
systematic uncertainties is given with the method TUnfoldSys.

### Why to use complex binning schemes

in literature on unfolding, the "standard" test case is a
one-dimensional distribution without underflow or overflow bins.
The migration matrix is almost diagonal.

<b>This "standard" case is rarely realized for real problems.</b>

Often one has to deal with multi-dimensional distributions.
In addition, there are underflow and overflow bins
or other background bins, possibly determined with the help of auxiliary
measurements.

In TUnfoldDensity, such complex binning schemes are handled with the help
of the class TUnfoldBinning. For both the measurement and the truth
there is a tree structure. The tree nodes may correspond to single
bins (e.g. nuisance parameters) or may hold multi-dimensional distributions.

For example, the "measurement" tree could have two leaves, one for
the primary distribution and one for auxiliary measurements.
Similarly, the "truth" tree could have two leaves, one for the
signal and one for the background.
Each of the leaves may then have a multi-dimensional distribution.

The class TUnfoldBinning takes care to map all bins of the
"measurement" to a one-dimensional vector y.
Similarly, the "truth" bins are mapped to the vector x.

### How to choose the regularisation settings

In TUnfoldDensity, two methods are implemented to determine tau**2

  1. ScanLcurve()  locate the tau where the L-curve plot has a "kink"
     this function is implemented in the TUnfold class
  2. ScanTau() finds the solution such that some variable
     (e.g. global correlation coefficient) is minimized.
     This function is implemented in the TUnfoldDensity class

Each of the algorithms has its own advantages and disadvantages.
The algorithm (1) does not work if the input data are too similar to the
MC prediction. Typical no-go cases of the L-curve scan are:

  - the number of measurements is too small (e.g. ny=nx)
  - the input data have no statistical fluctuations
 [identical MC events are used to fill the matrix of migrations
 and the vector y for a "closure test"]

The algorithm (2) only works if the variable does have a real minimum
as a function of tau. If global correlations are minimized, the situation
is as follows:
The matrix of migration typically introduces negative correlations.
The area constraint introduces some positive correlation.
Regularisation on the "size" introduces no correlation.
Regularisation on 1st or 2nd derivatives adds positive correlations.

For these reasons, "size" regularisation does not work well with
the tau-scan: the higher tau, the smaller rho, but there is no minimum.
As a result, large values of tau (too strong regularisation) are found.
In contrast, the tau-scan is expected to work better with 1st or 2nd
derivative regularisation, because at some point the negative
correlations from migrations are approximately cancelled by the
positive correlations from the regularisation conditions.

whichever algorithm is used, the output has to be checked:

  1. The L-curve should have approximate L-shape
     and the final choice of tau should not be at the very edge of the
     scanned region
  2. The scan result should have a well-defined minimum and the
     final choice of tau should sit right in the minimum


--------------------------------------------------------------------------------
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.

<b>Version 17.6, with updated doxygen comments and bug-fixes in TUnfoldBinning</b>

#### History:
  - Version 17.5, bug fix in TUnfold also corrects GetEmatrixSysUncorr()
  - Version 17.4, in parallel to changes in TUnfoldBinning
  - Version 17.3, in parallel to changes in TUnfoldBinning
  - Version 17.2, with new options 'N' and 'c' for axis regularisation steering
  - Version 17.1, add scan type RhoSquare, small bug fixes with useAxisBinning
  - Version 17.0, support for density regularisation, complex binning schemes, tau scan
*/

#include "TUnfoldDensity.h"
#include <TMath.h>
#include <TVectorD.h>
#include "TGraph.h"
#include <iostream>
#include <map>

//#define DEBUG

#ifdef DEBUG
using namespace std;
#endif

ClassImp(TUnfoldDensity);

TUnfoldDensity::~TUnfoldDensity(void)
{
   // clean up
   if(fOwnedOutputBins) delete fOwnedOutputBins;
   if(fOwnedInputBins) delete fOwnedInputBins;
   if(fRegularisationConditions) delete fRegularisationConditions;
}

////////////////////////////////////////////////////////////////////////////////
/// Only for use by root streamer or derived classes.

TUnfoldDensity::TUnfoldDensity(void)
{
   fConstOutputBins=0;
   fConstInputBins=0;
   fOwnedOutputBins=0;
   fOwnedInputBins=0;
   fRegularisationConditions=0;
}

////////////////////////////////////////////////////////////////////////////////
/// Eet up response matrix A, uncorrelated uncertainties of A,
/// regularisation scheme and binning schemes.
///
/// \param[in] hist_A matrix that describes the migrations
/// \param[in] histmap mapping of the histogram axes to the unfolding output
/// \param[in] regmode (default=kRegModeSize) global regularisation mode
/// \param[in] constraint (default=kEConstraintArea) type of constraint
/// \param[in] densityMode (default=kDensityModeBinWidthAndUser)
/// regularisation scale factors to construct the matrix L
/// \param[in] outputBins (default=0) binning scheme for truth (unfolding output)
/// \param[in] inputBins (default=0) binning scheme for measurement (unfolding
/// input)
/// \param[in] regularisationDistribution (default=0) selection of
/// regularized distribution
/// \param[in] regularisationAxisSteering (default=0) detailed
/// regularisation steering for selected distribution
///
/// The parameters <b>hist_A, histmap, constraint</b> are
/// explained with the TUnfoldSys constructor.
///
/// The parameters <b>outputBins,inputBins</b> set the binning
/// schemes. If these arguments are zero, simple binning schemes are
/// constructed which correspond to the axes of the histogram
/// <b>hist_A</b>.
///
/// The parameters
/// <b>regmode, densityMode, regularisationDistribution, regularisationAxisSteering</b>
/// together control how the initial matrix L of regularisation conditions
/// is constructed. as explained in RegularizeDistribution().

TUnfoldDensity::TUnfoldDensity
(const TH2 *hist_A, EHistMap histmap,ERegMode regmode,EConstraint constraint,
 EDensityMode densityMode,const TUnfoldBinning *outputBins,
 const TUnfoldBinning *inputBins,const char *regularisationDistribution,
 const char *regularisationAxisSteering) :
   TUnfoldSys(hist_A,histmap,kRegModeNone,constraint)
{
   fRegularisationConditions=0;
   // set up binning schemes
   fConstOutputBins = outputBins;
   fOwnedOutputBins = 0;
   TAxis const *genAxis,*detAxis;
   if(histmap==kHistMapOutputHoriz) {
      genAxis=hist_A->GetXaxis();
      detAxis=hist_A->GetYaxis();
   } else {
      genAxis=hist_A->GetYaxis();
      detAxis=hist_A->GetXaxis();
   }
   if(!fConstOutputBins) {
      // underflow and overflow are included in the binning scheme
      // They may be used on generator level
      fOwnedOutputBins =
         new TUnfoldBinning(*genAxis,1,1);
      fConstOutputBins = fOwnedOutputBins;
   }
   // check whether binning scheme is valid
   if(fConstOutputBins->GetParentNode()) {
      Error("TUnfoldDensity",
            "Invalid output binning scheme (node is not the root node)");
   }
   fConstInputBins = inputBins;
   fOwnedInputBins = 0;
   if(!fConstInputBins) {
      // underflow and overflow are not included in the binning scheme
      // They are still used to count events which have not been reconstructed
      fOwnedInputBins =
         new TUnfoldBinning(*detAxis,0,0);
      fConstInputBins = fOwnedInputBins;
   }
   if(fConstInputBins->GetParentNode()) {
      Error("TUnfoldDensity",
            "Invalid input binning scheme (node is not the root node)");
   }
   // check whether binning scheme matches with the histogram
   // in terms of total number of bins
   Int_t nOut=genAxis->GetNbins();
   Int_t nOutMappedT=TMath::Abs(fConstOutputBins->GetTH1xNumberOfBins(kTRUE));
   Int_t nOutMappedF=TMath::Abs(fConstOutputBins->GetTH1xNumberOfBins
                                (fOwnedOutputBins));
   if((nOutMappedT!= nOut)&&(nOutMappedF!=nOut)) {
      Error("TUnfoldDensity",
            "Output binning incompatible number of bins: axis %d binning scheme %d (%d)",
            nOut,nOutMappedT,nOutMappedF);
   }
   // check whether binning scheme matches with the histogram
   Int_t nInput=detAxis->GetNbins();
   Int_t nInputMappedT=TMath::Abs(fConstInputBins->GetTH1xNumberOfBins(kTRUE));
   Int_t nInputMappedF=TMath::Abs(fConstInputBins->GetTH1xNumberOfBins
                                  (fOwnedInputBins));
   if((nInputMappedT!= nInput)&&(nInputMappedF!= nInput)) {
      Error("TUnfoldDensity",
            "Input binning incompatible number of bins:axis %d binning scheme %d (%d) ",
            nInput,nInputMappedT,nInputMappedF);
   }

   // report detailed list of excluded bins
   for (Int_t ix = 0; ix <= nOut+1; ix++) {
      if(fHistToX[ix]<0) {
   Info("TUnfold","*NOT* unfolding bin %s",(char const *)GetOutputBinName(ix));
      }
   }

   // set up the regularisation here
   if(regmode !=kRegModeNone) {
      RegularizeDistribution
   (regmode,densityMode,regularisationDistribution,
    regularisationAxisSteering);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get bin name of an output bin.
///
/// \param[in] iBinX bin number
///
/// Return value: name of the bin. The name is constructed from the
/// entries in the binning scheme and includes information about the
/// bin borders etc.

TString TUnfoldDensity::GetOutputBinName(Int_t iBinX) const {
   if(!fConstOutputBins) return TUnfold::GetOutputBinName(iBinX);
   else return fConstOutputBins->GetBinName(iBinX);
}

////////////////////////////////////////////////////////////////////////////////
/// Density correction factor for a given bin.
///
/// \param[in]  densityMode type of factor to calculate
/// \param[in]  iBin  bin number
///
/// return a multiplicative factor, for scaling the regularisation
/// conditions from this bin.
///
/// For densityMode=kDensityModeNone the factor is set to unity.
/// For densityMode=kDensityModeBinWidth
/// the factor is set to 1/binArea
/// where the binArea is the product of the bin widths in all
/// dimensions. If the width of a bin is zero or can not be
/// determined, the factor is set to zero.
/// For densityMode=kDensityModeUser the factor is determined from the
///  method TUnfoldBinning::GetBinFactor().
/// For densityMode=kDensityModeBinWidthAndUser, the results of
/// kDensityModeBinWidth and kDensityModeUser are multiplied.

Double_t TUnfoldDensity::GetDensityFactor
(EDensityMode densityMode,Int_t iBin) const
{
   Double_t factor=1.0;
   if((densityMode == kDensityModeBinWidth)||
      (densityMode == kDensityModeBinWidthAndUser)) {
      Double_t binSize=fConstOutputBins->GetBinSize(iBin);
      if(binSize>0.0) factor /= binSize;
      else factor=0.0;
   }
   if((densityMode == kDensityModeUser)||
      (densityMode == kDensityModeBinWidthAndUser)) {
      factor *= fConstOutputBins->GetBinFactor(iBin);
   }
   return factor;
}

////////////////////////////////////////////////////////////////////////////////
/// Set up regularisation conditions.
///
/// \param[in] regmode basic regularisation mode (see class TUnfold)
/// \param[in] densityMode how to apply bin-wise factors
/// \param[in] distribution  name of the TUnfoldBinning node for which
/// the regularisation conditions shall be set (zero matches all nodes)
/// \param[in] axisSteering  regularisation fine-tuning
///
/// <b>axisSteering</b> is a string with several tokens, separated by
/// a semicolon: `"axisName[options];axisName[options];..."`.
///
///  - <b>axisName</b>:
///    the name of an axis. The special name * matches all.
///    So the argument <b>distribution</b> selects one (or all)
///    distributions. Within the selected distribution(s), steering options may be
///    specified for each axis (or for all axes) to define the
///    regularisation conditions.
///  - <b>options</b>
///    one or several character as follows:
///    - u : exclude underflow bin from derivatives along this axis
///    - o : exclude overflow bin from derivatives along this axis
///    - U : exclude underflow bin
///    - O : exclude overflow bin
///    - b : use bin width for derivative calculation
///    - B : same as 'b', in addition normalize to average bin width
///    - N : completely exclude derivatives along this axis
///    - p : axis is periodic (e.g. azimuthal angle), so
///          include derivatives built from combinations involving bins at
///          both ends of the axis "wrap around"
///
/// example:  <b>axisSteering</b>=`"*[UOB]"` uses bin widths to calculate
/// derivatives but underflow/overflow bins are not regularized

void TUnfoldDensity::RegularizeDistribution
(ERegMode regmode,EDensityMode densityMode,const char *distribution,
 const char *axisSteering)
{

   RegularizeDistributionRecursive(GetOutputBinning(),regmode,densityMode,
                                   distribution,axisSteering);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively add regularisation conditions for this node and its children.
///
/// \param[in] binning current node
/// \param[in] regmode regularisation mode
/// \param[in] densityMode type of regularisation scaling
/// \param[in] distribution target distribution(s) name
/// \param[in] axisSteering steering within the target distribution(s)

void TUnfoldDensity::RegularizeDistributionRecursive
(const TUnfoldBinning *binning,ERegMode regmode,
 EDensityMode densityMode,const char *distribution,const char *axisSteering) {
   if((!distribution)|| !TString(distribution).CompareTo(binning->GetName())) {
      RegularizeOneDistribution(binning,regmode,densityMode,axisSteering);
   }
   for(const TUnfoldBinning *child=binning->GetChildNode();child;
       child=child->GetNextNode()) {
      RegularizeDistributionRecursive(child,regmode,densityMode,distribution,
                                      axisSteering);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Regularize the distribution of the given node.
///
/// \param[in] binning current node
/// \param[in] regmode regularisation mode
/// \param[in] densityMode type of regularisation scaling
/// \param[in] axisSteering detailed steering for the axes of the distribution

void TUnfoldDensity::RegularizeOneDistribution
(const TUnfoldBinning *binning,ERegMode regmode,
 EDensityMode densityMode,const char *axisSteering)
{
#ifdef DEBUG
   cout<<"TUnfoldDensity::RegularizeOneDistribution node="
       <<binning->GetName()<<" "<<regmode<<" "<<densityMode
       <<" "<<(axisSteering ? axisSteering : "")<<"\n";
#endif
   if(!fRegularisationConditions)
      fRegularisationConditions=new TUnfoldBinning("regularisation");

   TUnfoldBinning *thisRegularisationBinning=
      fRegularisationConditions->AddBinning(binning->GetName());

   // decode steering
   Int_t isOptionGiven[8];
   binning->DecodeAxisSteering(axisSteering,"uUoObBpN",isOptionGiven);
   // U implies u
   isOptionGiven[0] |= isOptionGiven[1];
   // O implies o
   isOptionGiven[2] |= isOptionGiven[3];
   // B implies b
   isOptionGiven[4] |= isOptionGiven[5];
   // option N is removed if any other option is on
   for(Int_t i=0;i<7;i++) {
      isOptionGiven[7] &= ~isOptionGiven[i];
   }
   // option "c" does not work with options UuOo
   if(isOptionGiven[6] & (isOptionGiven[0]|isOptionGiven[2]) ) {
      Error("RegularizeOneDistribution",
            "axis steering %s is not valid",axisSteering);
   }
#ifdef DEBUG
   cout<<" "<<isOptionGiven[0]
       <<" "<<isOptionGiven[1]
       <<" "<<isOptionGiven[2]
       <<" "<<isOptionGiven[3]
       <<" "<<isOptionGiven[4]
       <<" "<<isOptionGiven[5]
       <<" "<<isOptionGiven[6]
       <<" "<<isOptionGiven[7]
       <<"\n";
#endif
   Info("RegularizeOneDistribution","regularizing %s regMode=%d"
   " densityMode=%d axisSteering=%s",
   binning->GetName(),(Int_t) regmode,(Int_t)densityMode,
   axisSteering ? axisSteering : "");
   Int_t startBin=binning->GetStartBin();
   Int_t endBin=startBin+ binning->GetDistributionNumberOfBins();
   std::vector<Double_t> factor(endBin-startBin);
   Int_t nbin=0;
   for(Int_t bin=startBin;bin<endBin;bin++) {
      factor[bin-startBin]=GetDensityFactor(densityMode,bin);
      if(factor[bin-startBin] !=0.0) nbin++;
   }
#ifdef DEBUG
   cout<<"initial number of bins "<<nbin<<"\n";
#endif
   Int_t dimension=binning->GetDistributionDimension();

   // decide whether to skip underflow/overflow bins
   nbin=0;
   for(Int_t bin=startBin;bin<endBin;bin++) {
      Int_t uStatus,oStatus;
      binning->GetBinUnderflowOverflowStatus(bin,&uStatus,&oStatus);
      if(uStatus & isOptionGiven[1]) factor[bin-startBin]=0.;
      if(oStatus & isOptionGiven[3]) factor[bin-startBin]=0.;
      if(factor[bin-startBin] !=0.0) nbin++;
   }
#ifdef DEBUG
   cout<<"after underflow/overflow bin removal "<<nbin<<"\n";
#endif
   if(regmode==kRegModeSize) {
      Int_t nRegBins=0;
      // regularize all bins of the distribution, possibly excluding
      // underflow/overflow bins
      for(Int_t bin=startBin;bin<endBin;bin++) {
         if(factor[bin-startBin]==0.0) continue;
         if(AddRegularisationCondition(bin,factor[bin-startBin])) {
            nRegBins++;
         }
      }
      if(nRegBins) {
         thisRegularisationBinning->AddBinning("size",nRegBins);
      }
   } else if((regmode==kRegModeDerivative)||(regmode==kRegModeCurvature)) {
      for(Int_t direction=0;direction<dimension;direction++) {
         // for each direction
         Int_t nRegBins=0;
         Int_t directionMask=(1<<direction);
         if(isOptionGiven[7] & directionMask) {
#ifdef DEBUG
            cout<<"skip direction "<<direction<<"\n";
#endif
            continue;
         }
         Double_t binDistanceNormalisation=
            (isOptionGiven[5] & directionMask)  ?
            binning->GetDistributionAverageBinSize
            (direction,isOptionGiven[0] & directionMask,
             isOptionGiven[2] & directionMask) : 1.0;
         for(Int_t bin=startBin;bin<endBin;bin++) {
            // check whether bin is excluded
            if(factor[bin-startBin]==0.0) continue;
            // for each bin, find the neighbour bins
            Int_t iPrev,iNext;
            Double_t distPrev,distNext;
            Int_t error=binning->GetBinNeighbours
               (bin,direction,&iPrev,&distPrev,&iNext,&distNext,
                isOptionGiven[6] & directionMask);
            if(error) {
               Error("RegularizeOneDistribution",
                     "invalid option %s (isPeriodic) for axis %s"
                     " (has underflow or overflow)",axisSteering,
                     binning->GetDistributionAxisLabel(direction).Data());
            }
            if((regmode==kRegModeDerivative)&&(iNext>=0)) {
               Double_t f0 = -factor[bin-startBin];
               Double_t f1 = factor[iNext-startBin];
               if(isOptionGiven[4] & directionMask) {
                  if(distNext>0.0) {
                     f0 *= binDistanceNormalisation/distNext;
                     f1 *= binDistanceNormalisation/distNext;
                  } else {
                     f0=0.;
                     f1=0.;
                  }
               }
               if((f0==0.0)||(f1==0.0)) continue;
               if(AddRegularisationCondition(bin,f0,iNext,f1)) {
                  nRegBins++;
#ifdef DEBUG
                  std::cout<<"Added Reg: bin "<<bin<<" "<<f0
                           <<" next: "<<iNext<<" "<<f1<<"\n";
#endif
               }
            } else if((regmode==kRegModeCurvature)&&(iPrev>=0)&&(iNext>=0)) {
               Double_t f0 = factor[iPrev-startBin];
               Double_t f1 = -factor[bin-startBin];
               Double_t f2 = factor[iNext-startBin];
               if(isOptionGiven[4] & directionMask) {
                  if((distPrev<0.)&&(distNext>0.)) {
                     distPrev= -distPrev;
                     Double_t f=TMath::Power(binDistanceNormalisation,2.)/
                        (distPrev+distNext);
                     f0 *= f/distPrev;
                     f1 *= f*(1./distPrev+1./distNext);
                     f2 *= f/distNext;
                  } else {
                     f0=0.;
                     f1=0.;
                     f2=0.;
                  }
               }
               if((f0==0.0)||(f1==0.0)||(f2==0.0)) continue;
               if(AddRegularisationCondition(iPrev,f0,bin,f1,iNext,f2)) {
                  nRegBins++;
#ifdef DEBUG
                  std::cout<<"Added Reg: prev "<<iPrev<<" "<<f0
                           <<" bin: "<<bin<<" "<<f1
                           <<" next: "<<iNext<<" "<<f2<<"\n";
#endif
               }
            }
         }
         if(nRegBins) {
            TString name;
            if(regmode==kRegModeDerivative) {
               name="derivative_";
            } else if(regmode==kRegModeCurvature) {
               name="curvature_";
            }
            name +=  binning->GetDistributionAxisLabel(direction);
            thisRegularisationBinning->AddBinning(name,nRegBins);
         }
      }
   }
#ifdef DEBUG
   //fLsquared->Print();
#endif
}

///////////////////////////////////////////////////////////////////////
/// retrieve unfolding result as a new histogram
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// return value: pointer to a new histogram.  If
/// <b>useAxisBinning</b> is set and if the selected distribution fits
/// into a root histogram (1,2,3 dimensions) then return a histogram
/// with the proper binning on each axis. Otherwise, return a 1D
/// histogram with equidistant binning. If the histogram title is
/// zero, a title is assigned automatically.
///
/// The <b>axisSteering</b> is defines as follows: "axis[mode];axis[mode];..."
/// where:
///
///   - axis = name of an axis or *
///   - mode = any combination of the letters CUO0123456789
///
///   - C collapse axis into one bin (add up results). If
///     any of the numbers 0-9 are given in addition, only these bins are added up.
///     Here bins are counted from zero, whereas in root, bins are counted
///     from 1. Obviously, this only works for up to 10 bins.
///   - U discard underflow bin
///   - O discard overflow bin
///
/// examples: imagine the binning has two axis, named x and y.
///
///   - "*[UO]" exclude underflow and overflow bins for all axis.
///     So here a TH2 is returned but all undeflow and overflow bins are empty
///   - "x[UOC123]" integrate over the variable x but only using the
///     bins 1,2,3 and not the underflow and overflow in x.
///     Here a TH1 is returned, the axis is labelled "y" and
///     the underflow and overflow (in y) are filled. However only the x-bins
///     1,2,3 are used to determine the content.
///   - "x[C];y[UO]" integrate over the variable x, including
///     underflow and overflow but exclude underflow and overflow in y.
///     Here a TH1 is returned, the axis is labelled "y". The underflow
///     and overflow in y are empty.

TH1 *TUnfoldDensity::GetOutput
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning) const
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
     TUnfoldSys::GetOutput(r,binMap);
   }
   if(binMap) {
     delete [] binMap;
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve bias vector as a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetBias
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning) const
{
  TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetBias(r,binMap);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve unfolding result folded back as a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
/// \param[in] addBgr (default=false) if true, include the background
/// contribution (useful for direct comparison to data)
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetFoldedOutput
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning,Bool_t addBgr) const
{
   TUnfoldBinning const *binning=fConstInputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetFoldedOutput(r,binMap);
      if(addBgr) {
         TUnfoldSys::GetBackground(r,0,binMap,0,kFALSE);
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve a background source in a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] bgrSource the background source to retrieve
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
/// \param[in] includeError (default=3) type of background errors to
/// be included (+1 uncorrelated bgr errors, +2 correlated bgr errors)
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetBackground
(const char *histogramName,const char *bgrSource,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,Bool_t useAxisBinning,
 Int_t includeError) const
{
   TUnfoldBinning const *binning=fConstInputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetBackground(r,bgrSource,binMap,includeError,kTRUE);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve input distribution in a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetInput
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning) const
{
   TUnfoldBinning const *binning=fConstInputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetInput(r,binMap);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve global correlation coefficients including all uncertainty sources.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
/// \param[out] ematInv (default=0) to return the inverse covariance matrix
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments. The inverse of the covariance matrix
/// is stored in a new histogram returned by <b>ematInv</b> if that
/// pointer is non-zero.

TH1 *TUnfoldDensity::GetRhoItotal
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning,TH2 **ematInv) {
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TH2 *invEmat=0;
      if(ematInv) {
         if(r->GetDimension()==1) {
            TString ematName(histogramName);
            ematName += "_inverseEMAT";
            Int_t *binMap2D=0;
            invEmat=binning->CreateErrorMatrixHistogram
               (ematName,useAxisBinning,&binMap2D,histogramTitle,
                axisSteering);
            if(binMap2D) delete [] binMap2D;
         } else {
            Error("GetRhoItotal",
                  "can not return inverse of error matrix for this binning");
         }
      }
      TUnfoldSys::GetRhoItotal(r,binMap,invEmat);
      if(invEmat) {
         *ematInv=invEmat;
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve global correlation coefficients including input
/// (statistical) and background uncertainties.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
/// \param[out] ematInv (default=0) to return the inverse covariance matrix
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments. The inverse of the covariance matrix
/// is stored in a new histogram returned by <b>ematInv</b> if that
/// pointer is non-zero.

TH1 *TUnfoldDensity::GetRhoIstatbgr
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning,TH2 **ematInv) {
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TH2 *invEmat=0;
      if(ematInv) {
         if(r->GetDimension()==1) {
            TString ematName(histogramName);
            ematName += "_inverseEMAT";
            Int_t *binMap2D=0;
            invEmat=binning->CreateErrorMatrixHistogram
               (ematName,useAxisBinning,&binMap2D,histogramTitle,
                axisSteering);
            if(binMap2D) delete [] binMap2D;
         } else {
            Error("GetRhoItotal",
                  "can not return inverse of error matrix for this binning");
         }
      }
      TUnfoldSys::GetRhoI(r,binMap,invEmat);
      if(invEmat) {
         *ematInv=invEmat;
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve a correlated systematic 1-sigma shift.
///
/// \param[in] source identifier of the systematic uncertainty source
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetDeltaSysSource
(const char *source,const char *histogramName,
 const char *histogramTitle,const char *distributionName,
 const char *axisSteering,Bool_t useAxisBinning) {
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      if(!TUnfoldSys::GetDeltaSysSource(r,source,binMap)) {
         delete r;
         r=0;
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve systematic 1-sigma shift corresponding to a background
/// scale uncertainty.
///
/// \param[in] bgrSource identifier of the background
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetDeltaSysBackgroundScale
(const char *bgrSource,const char *histogramName,
 const char *histogramTitle,const char *distributionName,
 const char *axisSteering,Bool_t useAxisBinning) {
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      if(!TUnfoldSys::GetDeltaSysBackgroundScale(r,bgrSource,binMap)) {
         delete r;
         r=0;
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve 1-sigma shift corresponding to the previously specified uncertainty
/// on tau.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH1 *TUnfoldDensity::GetDeltaSysTau
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,Bool_t useAxisBinning)
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH1 *r=binning->CreateHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      if(!TUnfoldSys::GetDeltaSysTau(r,binMap)) {
         delete r;
         r=0;
      }
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve correlation coefficients, including all uncertainties.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH2 *TUnfoldDensity::GetRhoIJtotal
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning)
{
   TH2 *r=GetEmatrixTotal
      (histogramName,histogramTitle,distributionName,
       axisSteering,useAxisBinning);
   if(r) {
      for(Int_t i=0;i<=r->GetNbinsX()+1;i++) {
         Double_t e_i=r->GetBinContent(i,i);
         if(e_i>0.0) e_i=TMath::Sqrt(e_i);
         else e_i=0.0;
         for(Int_t j=0;j<=r->GetNbinsY()+1;j++) {
      if(i==j) continue;
            Double_t e_j=r->GetBinContent(j,j);
            if(e_j>0.0) e_j=TMath::Sqrt(e_j);
            else e_j=0.0;
            Double_t e_ij=r->GetBinContent(i,j);
            if((e_i>0.0)&&(e_j>0.0)) {
               r->SetBinContent(i,j,e_ij/e_i/e_j);
            } else {
               r->SetBinContent(i,j,0.0);
            }
         }
      }
      for(Int_t i=0;i<=r->GetNbinsX()+1;i++) {
   if(r->GetBinContent(i,i)>0.0) {
     r->SetBinContent(i,i,1.0);
   } else {
     r->SetBinContent(i,i,0.0);
   }
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve covariance contribution from uncorrelated (statistical)
/// uncertainties of the response matrix.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH2 *TUnfoldDensity::GetEmatrixSysUncorr
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning)
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH2 *r=binning->CreateErrorMatrixHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetEmatrixSysUncorr(r,binMap);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve covariance contribution from uncorrelated background uncertainties.
///
/// \param[in] bgrSource identifier of the background
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH2 *TUnfoldDensity::GetEmatrixSysBackgroundUncorr
(const char *bgrSource,const char *histogramName,
 const char *histogramTitle,const char *distributionName,
 const char *axisSteering,Bool_t useAxisBinning)
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH2 *r=binning->CreateErrorMatrixHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetEmatrixSysBackgroundUncorr(r,bgrSource,binMap,kFALSE);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get covariance contribution from the input uncertainties (data
/// statistical uncertainties).
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH2 *TUnfoldDensity::GetEmatrixInput
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning)
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH2 *r=binning->CreateErrorMatrixHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetEmatrixInput(r,binMap);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get matrix of probabilities in a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. if histogramTitle is null, choose a title
/// automatically.

TH2 *TUnfoldDensity::GetProbabilityMatrix
(const char *histogramName,const char *histogramTitle,
 Bool_t useAxisBinning) const
{
   TH2 *r=TUnfoldBinning::CreateHistogramOfMigrations
      (fConstOutputBins,fConstInputBins,histogramName,
       useAxisBinning,useAxisBinning,histogramTitle);
   TUnfold::GetProbabilityMatrix(r,kHistMapOutputHoriz);
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get covariance matrix including all contributions.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] distributionName (default=0) identifier of the distribution to be extracted
/// \param[in] axisSteering (default=0) detailed steering within the extracted
/// distribution
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. See method GetOutput() for a detailed
/// description of the arguments

TH2 *TUnfoldDensity::GetEmatrixTotal
(const char *histogramName,const char *histogramTitle,
 const char *distributionName,const char *axisSteering,
 Bool_t useAxisBinning)
{
   TUnfoldBinning const *binning=fConstOutputBins->FindNode(distributionName);
   Int_t *binMap=0;
   TH2 *r=binning->CreateErrorMatrixHistogram
      (histogramName,useAxisBinning,&binMap,histogramTitle,axisSteering);
   if(r) {
      TUnfoldSys::GetEmatrixTotal(r,binMap);
   }
   if(binMap) delete [] binMap;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Access matrix of regularisation conditions in a new histogram.
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
/// \param[in] useAxisBinning (default=true) if set to true, try to extract a histogram with
/// proper binning and axis labels
///
/// returns a new histogram. if histogramTitle is null, choose a title
/// automatically.

TH2 *TUnfoldDensity::GetL
(const char *histogramName,const char *histogramTitle,Bool_t useAxisBinning)
{
   if(fRegularisationConditions &&
      (fRegularisationConditions->GetEndBin()-
       fRegularisationConditions->GetStartBin()!= fL->GetNrows())) {
      Warning("GetL",
              "remove invalid scheme of regularisation conditions %d %d",
              fRegularisationConditions->GetEndBin(),fL->GetNrows());
      delete fRegularisationConditions;
      fRegularisationConditions=0;
   }
   if(!fRegularisationConditions) {
      fRegularisationConditions=new TUnfoldBinning("regularisation",fL->GetNrows());
      Warning("GetL","create flat regularisation conditions scheme");
   }
   TH2 *r=TUnfoldBinning::CreateHistogramOfMigrations
      (fConstOutputBins,fRegularisationConditions,histogramName,
       useAxisBinning,useAxisBinning,histogramTitle);
   TUnfold::GetL(r);
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get regularisation conditions multiplied by result vector minus bias
///   L(x-biasScale*biasVector).
///
/// \param[in] histogramName name of the histogram
/// \param[in] histogramTitle (default=0) title of the histogram
///
/// returns a new histogram.
/// This is a measure of the level of regularization required per
/// regularisation condition.
/// If there are (negative or positive) spikes,
/// these regularisation conditions dominate
/// over the other regularisation conditions and may introduce
/// the largest biases.

TH1 *TUnfoldDensity::GetLxMinusBias
(const char *histogramName,const char *histogramTitle)
{
   TMatrixD dx(*GetX(), TMatrixD::kMinus, fBiasScale * (*fX0));
   TMatrixDSparse *Ldx=MultiplyMSparseM(fL,&dx);
   if(fRegularisationConditions &&
      (fRegularisationConditions->GetEndBin()-
       fRegularisationConditions->GetStartBin()!= fL->GetNrows())) {
      Warning("GetLxMinusBias",
              "remove invalid scheme of regularisation conditions %d %d",
              fRegularisationConditions->GetEndBin(),fL->GetNrows());
      delete fRegularisationConditions;
      fRegularisationConditions=0;
   }
   if(!fRegularisationConditions) {
      fRegularisationConditions=new TUnfoldBinning("regularisation",fL->GetNrows());
      Warning("GetLxMinusBias","create flat regularisation conditions scheme");
   }
   TH1 *r=fRegularisationConditions->CreateHistogram
      (histogramName,kFALSE,0,histogramTitle);
   const Int_t *Ldx_rows=Ldx->GetRowIndexArray();
   const Double_t *Ldx_data=Ldx->GetMatrixArray();
   for(Int_t row=0;row<Ldx->GetNrows();row++) {
      if(Ldx_rows[row]<Ldx_rows[row+1]) {
         r->SetBinContent(row+1,Ldx_data[Ldx_rows[row]]);
      }
   }
   delete Ldx;
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a binning node for the input (measured) quantities.
///
/// \param[in] distributionName (default=0) distribution to look
/// for. if zero, return the root node
///
/// returns: pointer to a TUnfoldBinning object or zero if not found

const TUnfoldBinning *TUnfoldDensity::GetInputBinning
(const char *distributionName) const
{
   // find binning scheme, input bins
   //   distributionName : the distribution to locate
   return fConstInputBins->FindNode(distributionName);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a binning node for the unfolded (truth level) quantities.
///
/// \param[in] distributionName (default=0) distribution to look
/// for. if zero, return the root node
///
/// returns: pointer to a TUnfoldBinning object or zero if not found

const TUnfoldBinning *TUnfoldDensity::GetOutputBinning
(const char *distributionName) const
{
   // find binning scheme, output bins
   //   distributionName : the distribution to locate
   return fConstOutputBins->FindNode(distributionName);
}

////////////////////////////////////////////////////////////////////////////////
/// Scan a function wrt tau and determine the minimum.
///
/// \param[in] nPoint number of points to be scanned
/// \param[in] tauMin smallest tau value to study
/// \param[in] tauMax largest tau value to study
/// \param[out] scanResult the scanned function wrt log(tau)
/// \param[in] mode 1st parameter for the scan function
/// \param[in] distribution 2nd parameter for the scan function
/// \param[in] axisSteering 3rd parameter for the scan function
/// \param[out] lCurvePlot for monitoring, shows the L-curve
/// \param[out] logTauXPlot for monitoring, L-curve(X) as a function of log(tau)
/// \param[out] logTauYPlot for monitoring, L-curve(Y) as a function of log(tau)
///
/// Return value: the coordinate number on the curve <b>scanResult</b>
/// which corresponds to the minimum
///
/// The function is scanned by repeating the following steps <b>nPoint</b>
/// times
///
///   1. Choose a value of tau
///   2. Perform the unfolding for this choice of tau, DoUnfold(tau)
///   3. Determine the scan variable GetScanVariable()
///
/// The method  GetScanVariable() defines scans of correlation
/// coefficients, where <b>mode</b> is chosen from the enum
/// EScanTauMode. In addition one may set <b>distribution</b>
/// and/or <b>projectionMode</b> to refine the calculation of
/// correlations (e.g. restrict the calculation to the signal
/// distribution and/or exclude underflow and overflow bins).
/// See the documentation of GetScanVariable() for details.
/// Alternative scan variables may be defined by overriding the
/// GetScanVariable() method.
///
/// Automatic choice of scan range: if (tauMin,tauMax) do not
/// correspond to a valid tau range (e.g. tauMin=tauMax=0.0) then
/// the tau range is determined automatically. Use with care!

Int_t TUnfoldDensity::ScanTau
(Int_t nPoint,Double_t tauMin,Double_t tauMax,TSpline **scanResult,
 Int_t mode,const char *distribution,const char *axisSteering,
 TGraph **lCurvePlot,TSpline **logTauXPlot,TSpline **logTauYPlot)
{
   typedef std::map<Double_t,Double_t> TauScan_t;
   typedef std::map<Double_t,std::pair<Double_t,Double_t> > LCurve_t;
   TauScan_t curve;
   LCurve_t lcurve;

   //==========================================================
   // algorithm:
   //  (1) do the unfolding for nPoint-1 points
   //      and store the results in the map
   //        curve
   //    (1a) store minimum and maximum tau to curve
   //    (1b) insert additional points, until nPoint-1 values
   //          have been calculated
   //
   //  (2) determine the best choice of tau
   //      do the unfolding for this point
   //      and store the result in
   //        curve
   //  (3) return the result in scanResult

   //==========================================================
   //  (1) do the unfolding for nPoint-1 points
   //      and store the results in
   //        curve
   //    (1a) store minimum and maximum tau to curve

   if((tauMin<=0)||(tauMax<=0.0)||(tauMin>=tauMax)) {
      // here no range is given, has to be determined automatically
      // the maximum tau is determined from the chi**2 values
      // observed from unfolding without regularization

      // first unfolding, without regularisation
      DoUnfold(0.0);

      // if the number of degrees of freedom is too small, create an error
      if(GetNdf()<=0) {
         Error("ScanTau","too few input bins, NDF<=0 %d",GetNdf());
      }
      Double_t X0=GetLcurveX();
      Double_t Y0=GetLcurveY();
      Double_t y0=GetScanVariable(mode,distribution,axisSteering);
      Info("ScanTau","logtau=-Infinity y=%lf X=%lf Y=%lf",y0,X0,Y0);
      {
         // unfolding guess maximum tau and store it
         Double_t logTau=
            0.5*(TMath::Log10(GetChi2A()+3.*TMath::Sqrt(GetNdf()+1.0))
                 -GetLcurveY());
         DoUnfold(TMath::Power(10.,logTau));
         if((!TMath::Finite(GetLcurveX())) ||(!TMath::Finite(GetLcurveY()))) {
            Fatal("ScanTau","problem (missing regularisation?) X=%f Y=%f",
                  GetLcurveX(),GetLcurveY());
         }
         Double_t y=GetScanVariable(mode,distribution,axisSteering);
         curve[logTau]=y;
         lcurve[logTau]=std::make_pair(GetLcurveX(),GetLcurveY());
         Info("ScanTau","logtau=%lf y=%lf X=%lf Y=%lf",logTau,y,
              GetLcurveX(),GetLcurveY());
      }
      // minimum tau is chosen such that it is less than
      // 1% different from the case of no regularisation
      // here, several points are inserted as needed
      while(((int)curve.size()<nPoint-1)&&
            ((TMath::Abs(GetLcurveX()-X0)>0.00432)||
             (TMath::Abs(GetLcurveY()-Y0)>0.00432))) {
         Double_t logTau=(*curve.begin()).first-0.5;
         DoUnfold(TMath::Power(10.,logTau));
         Double_t y=GetScanVariable(mode,distribution,axisSteering);
         curve[logTau]=y;
         lcurve[logTau]=std::make_pair(GetLcurveX(),GetLcurveY());
         Info("ScanTay","logtau=%lf y=%lf X=%lf Y=%lf",logTau,y,
              GetLcurveX(),GetLcurveY());
      }
   } else {
      Double_t logTauMin=TMath::Log10(tauMin);
      Double_t logTauMax=TMath::Log10(tauMax);
      if(nPoint>1) {
         // insert maximum tau
         DoUnfold(TMath::Power(10.,logTauMax));
         Double_t y=GetScanVariable(mode,distribution,axisSteering);
         curve[logTauMax]=y;
         lcurve[logTauMax]=std::make_pair(GetLcurveX(),GetLcurveY());
         Info("ScanTau","logtau=%lf y=%lf X=%lf Y=%lf",logTauMax,y,
              GetLcurveX(),GetLcurveY());
      }
      // insert minimum tau
      DoUnfold(TMath::Power(10.,logTauMin));
      Double_t y=GetScanVariable(mode,distribution,axisSteering);
      curve[logTauMin]=y;
      lcurve[logTauMin]=std::make_pair(GetLcurveX(),GetLcurveY());
      Info("ScanTau","logtau=%lf y=%lf X=%lf Y=%lf",logTauMin,y,
           GetLcurveX(),GetLcurveY());
   }

   //==========================================================
   //    (1b) insert additional points, until nPoint-1 values
   //          have been calculated
   while((int)curve.size()<nPoint-1) {
      // insert additional points
      // points are inserted such that the largest interval in log(tau)
      // is divided into two smaller intervals
      // however, there is a penalty term for subdividing intervals
      // which are far away from the minimum
      TauScan_t::const_iterator i0,i1;
      i0=curve.begin();
      // locate minimum
      Double_t logTauYMin=(*i0).first;
      Double_t yMin=(*i0).second;
      for (; i0 != curve.end(); ++i0) {
         if((*i0).second<yMin) {
            yMin=(*i0).second;
            logTauYMin=(*i0).first;
         }
      }
      // insert additional points such that large log(tau) intervals
      // near the minimum rho are divided into two
      i0=curve.begin();
      i1=i0;
      Double_t distMax=0.0;
      Double_t logTau=0.0;
      for (++i1; i1 != curve.end(); ++i1) {
         Double_t dist;
         // check size of rho interval
         dist=TMath::Abs((*i0).first-(*i1).first)
            // penalty term if distance from rhoMax is large
            +0.25*TMath::Power(0.5*((*i0).first+(*i1).first)-logTauYMin,2.)/
            ((*curve.rbegin()).first-(*curve.begin()).first)/nPoint;
         if((dist<=0.0)||(dist>distMax)) {
            distMax=dist;
            logTau=0.5*((*i0).first+(*i1).first);
         }
         i0=i1;
      }
      DoUnfold(TMath::Power(10.,logTau));
      Double_t y=GetScanVariable(mode,distribution,axisSteering);
      curve[logTau]=y;
      lcurve[logTau]=std::make_pair(GetLcurveX(),GetLcurveY());
      Info("ScanTau","logtau=%lf y=%lf X=%lf Y=%lf",logTau,y,
           GetLcurveX(),GetLcurveY());
   }

   //==========================================================
   //  (2) determine the best choice of tau
   //      do the unfolding for this point
   //      and store the result in
   //        curve

   Double_t cTmin=0.0;
   {
   Double_t *cTi=new Double_t[curve.size()];
   Double_t *cCi=new Double_t[curve.size()];
   Int_t n=0;
   for (TauScan_t::const_iterator i = curve.begin(); i != curve.end(); ++i) {
      cTi[n]=(*i).first;
      cCi[n]=(*i).second;
      n++;
   }
   // create rho Spline
   TSpline3 *splineC=new TSpline3("L curve curvature",cTi,cCi,n);
   // find the maximum of the curvature
   // if the parameter iskip is non-zero, then iskip points are
   // ignored when looking for the largest curvature
   // (there are problems with the curvature determined from the first
   //  few points of splineX,splineY in the algorithm above)
   Int_t iskip=0;
   if(n>3) iskip=1;
   if(n>6) iskip=2;
   Double_t cCmin=cCi[iskip];
   cTmin=cTi[iskip];
   for(Int_t i=iskip;i<n-1-iskip;i++) {
      // find minimum on this spline section
      // check boundary conditions for x[i+1]
      Double_t xMin=cTi[i+1];
      Double_t yMin=cCi[i+1];
      if(cCi[i]<yMin) {
         yMin=cCi[i];
         xMin=cTi[i];
      }
      // find minimum for x[i]<x<x[i+1]
      // get spline coefficients and solve equation
      //   derivative(x)==0
      Double_t x,y,b,c,d;
      splineC->GetCoeff(i,x,y,b,c,d);
      // coefficients of quadratic equation
      Double_t m_p_half=-c/(3.*d);
      Double_t q=b/(3.*d);
      Double_t discr=m_p_half*m_p_half-q;
      if(discr>=0.0) {
         // solution found
         discr=TMath::Sqrt(discr);
         Double_t xx;
         if(m_p_half>0.0) {
            xx = m_p_half + discr;
         } else {
            xx = m_p_half - discr;
         }
         Double_t dx=cTi[i+1]-x;
         // check first solution
         if((xx>0.0)&&(xx<dx)) {
            y=splineC->Eval(x+xx);
            if(y<yMin) {
               yMin=y;
               xMin=x+xx;
            }
         }
         // second solution
         if(xx !=0.0) {
            xx= q/xx;
         } else {
            xx=0.0;
         }
         // check second solution
         if((xx>0.0)&&(xx<dx)) {
            y=splineC->Eval(x+xx);
            if(y<yMin) {
               yMin=y;
               xMin=x+xx;
            }
         }
      }
      // check whether this local minimum is a global minimum
      if(yMin<cCmin) {
         cCmin=yMin;
         cTmin=xMin;
      }
   }
   delete splineC;
   delete[] cTi;
   delete[] cCi;
   }
   Double_t logTauFin=cTmin;
   DoUnfold(TMath::Power(10.,logTauFin));
   {
      Double_t y=GetScanVariable(mode,distribution,axisSteering);
      curve[logTauFin]=y;
      lcurve[logTauFin]=std::make_pair(GetLcurveX(),GetLcurveY());
      Info("ScanTau","Result logtau=%lf y=%lf X=%lf Y=%lf",logTauFin,y,
           GetLcurveX(),GetLcurveY());
   }
   //==========================================================
   //  (3) return the result in
   //       scanResult lCurve logTauX logTauY

   Int_t bestChoice=-1;
   if(curve.size()>0) {
      Double_t *y=new Double_t[curve.size()];
      Double_t *logT=new Double_t[curve.size()];
      int n=0;
      for (TauScan_t::const_iterator i = curve.begin(); i != curve.end(); ++i) {
         if(logTauFin==(*i).first) {
            bestChoice=n;
         }
         y[n]=(*i).second;
         logT[n]=(*i).first;
         n++;
      }
      if(scanResult) {
         TString name;
         name = TString::Format("scan(%d,",mode);
         if(distribution) name+= distribution;
         name += ",";
         if(axisSteering) name += axisSteering;
         name +=")";
         (*scanResult)=new TSpline3(name+"%log(tau)",logT,y,n);
      }
      delete[] y;
      delete[] logT;
   }
   if(lcurve.size()) {
      Double_t *logT=new Double_t[lcurve.size()];
      Double_t *x=new Double_t[lcurve.size()];
      Double_t *y=new Double_t[lcurve.size()];
      Int_t n=0;
      for (LCurve_t::const_iterator i = lcurve.begin(); i != lcurve.end(); ++i) {
         logT[n]=(*i).first;
         x[n]=(*i).second.first;
         y[n]=(*i).second.second;
         //cout<<logT[n]<<" "<< x[n]<<" "<<y[n]<<"\n";
         n++;
      }
      if(lCurvePlot) {
         *lCurvePlot=new TGraph(n,x,y);
         (*lCurvePlot)->SetTitle("L curve");
      }
      if(logTauXPlot)
         *logTauXPlot=new TSpline3("log(chi**2)%log(tau)",logT,x,n);
      if(logTauYPlot)
         *logTauYPlot=new TSpline3("log(reg.cond)%log(tau)",logT,y,n);
      delete [] y;
      delete [] x;
      delete [] logT;
   }
   return bestChoice;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the function for ScanTau().
///
/// \param[in] mode the variable to be calculated
/// \param[in] distribution distribution for which the variable
///              is to be calculated
/// \param[in] axisSteering detailed steering for selecting bins on
/// the axes of the distribution (see method GetRhoItotal())
///
/// return value: the scan result for the given choice of tau (for
/// which the unfolding was performed prior to calling this method)
///
/// In ScanTau() the unfolding is repeated for various choices of tau.
/// For each tau, after unfolding, GetScanVariable() is called to
/// determine the scan result for this choice of tau.
///
/// the following modes are implemented
///
///   - kEScanTauRhoAvg : average (stat+bgr) global correlation
///   - kEScanTauRhoSquaredAvg : average (stat+bgr) global correlation squared
///   - kEScanTauRhoMax : maximum (stat+bgr) global correlation
///   - kEScanTauRhoAvgSys : average (stat+bgr+sys) global correlation
///   - kEScanTauRhoAvgSquaredSys : average (stat+bgr+sys) global correlation squared
///   - kEScanTauRhoMaxSys : maximum (stat+bgr+sys) global correlation

Double_t TUnfoldDensity::GetScanVariable
(Int_t mode,const char *distribution,const char *axisSteering)
{
   Double_t r=0.0;
   TString name("GetScanVariable(");
   name += TString::Format("%d,",mode);
   if(distribution) name += distribution;
   name += ",";
   if(axisSteering) name += axisSteering;
   name += ")";
   TH1 *rhoi=0;
   if((mode==kEScanTauRhoAvg)||(mode==kEScanTauRhoMax)||
      (mode==kEScanTauRhoSquareAvg)) {
      rhoi=GetRhoIstatbgr(name,0,distribution,axisSteering,kFALSE);
   } else if((mode==kEScanTauRhoAvgSys)||(mode==kEScanTauRhoMaxSys)||
             (mode==kEScanTauRhoSquareAvgSys)) {
      rhoi=GetRhoItotal(name,0,distribution,axisSteering,kFALSE);
   }
   if(rhoi) {
      Double_t sum=0.0;
      Double_t sumSquare=0.0;
      Double_t rhoMax=0.0;
      Int_t n=0;
      for(Int_t i=0;i<=rhoi->GetNbinsX()+1;i++) {
         Double_t c=rhoi->GetBinContent(i);
         if(c>=0.) {
            if(c>rhoMax) rhoMax=c;
            sum += c;
            sumSquare += c*c;
            n ++;
         }
      }
      if((mode==kEScanTauRhoAvg)||(mode==kEScanTauRhoAvgSys)) {
         r=sum/n;
      } else if((mode==kEScanTauRhoSquareAvg)||
                (mode==kEScanTauRhoSquareAvgSys)) {
         r=sum/n;
      } else {
         r=rhoMax;
      }
      // cout<<r<<" "<<GetRhoAvg()<<" "<<GetRhoMax()<<"\n";
      delete rhoi;
   } else {
      Fatal("GetScanVariable","mode %d not implemented",mode);
   }
   return r;
}
