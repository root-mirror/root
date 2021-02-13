// @(#)root/unfold:$Id$
// Author: Stefan Schmitt DESY, 10/08/11

/** \class TUnfoldBinning
\ingroup Unfold
Binning schemes for use with the unfolding algorithm TUnfoldDensity.

Binning schemes are used to map analysis bins on a single histogram
axis and back. The analysis bins may include unconnected bins (e.g
nuisances for background normalisation) or various multidimensional
histograms (signal bins, differential background normalisation bins, etc).

If you use this software, please consider the following citation

<b>S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]</b>

Detailed documentation and updates are available on
http://www.desy.de/~sschmitt

### Functionality

The TUnfoldBinning objects are connected by a tree-like structure.
The structure does not hold any data, but is only responsible for
arranging the analysis bins in the proper order.
Each node of the tree is responsible for a group of bins. That group
may consist of

  -  several unconnected bins, each with a dedicated name.
  -  bins organized in a multidimensional distribution, defined by a
set of axes. The axes are defined by a number of bins N and by (N+1)
bin borders. In addition to the N bins inside there may be an underflow and an
overflow bin

Each bin has a "global" bin number, which can be found using the
GetGlobalBinNumber() methods. The global bin number 0 is reserved and
corresponds to the case where no bin is found in the
TUnfoldBinning tree.

### Use in the analysis
Booking histograms:

  - Define binning schemes on detector level and on truth level. This
can be done using the XML language, use the class TUnfoldBinningXML to
read the binning scheme. The TUnfoldBinning objects can be written to
a root file, preferentially together with the corresponding histograms.
  - For Monte Carlo, book histograms for the response matrix (detector
vs truth level) using the
method CreateHistogramOfMigrations()
  - For data and background, book histograms using the
"detector level" binning scheme and the method CreateHistogram()
  - (if required) for the data covariance matrix, book a histogram using the
"detector level" binning scheme and the method CreateErrorMatrixHistogram()
  - For truth histograms, book histograms using the
"truth level" binning scheme and the method CreateHistogram()

The histograms which are booked have all analysis bins arranged on one
axis (global bin number). TUnfoldBinning provides methods to locate
the global bin number:

  - Use the method FindNode() to locate a group of bins (e.g. signal,
control distribution, etc) by their name, then:
  - Use the method GetGlobalBinNumber() to locate a bin in a
distribution, then:
  - Use the TH1::Fill() method and the bin number to fill the
appropriate bin in one of the histograms booked above.

Unfolding: Specify the response matrix and the binning schemes when
constructing a TUnfoldDensity object. Tell TUnfoldDensity about the
data, background, systematic error histograms using the corresponding
methods of class TUnfoldDensity. Then run the unfolding. Use the
GetXXX() methods to retrieve the unfolding results into properly
binned multidimensional histograms.

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

<b>Version 17.6, bug fix to avoid possible crash in method
CreateHistogramOfMigrations(). Possible bug fix with NaN in GetGlobalBinNUmber() </b>

#### History:
  - Version 17.5, in parallel to changes in TUnfold
  - Version 17.4, bug fix with error handling
  - Version 17.3, bug fix with underflow/overflow bins
  - Version 17.2, with XML support, bug fix with bin map creation, isPeriodic option for neighbour bins
  - Version 17.1, in parallel to changes in TUnfold
  - Version 17.0, initial version, numbered in parallel to TUnfold

*/

#include "TUnfoldBinningXML.h"
#include <TVectorD.h>
#include <TAxis.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TIterator.h>
#include <iomanip>

// #define DEBUG

using namespace std;

ClassImp(TUnfoldBinning);

////////////////////////////////////////////////////////////////////////////////
// Destructor.

TUnfoldBinning::~TUnfoldBinning(void)
{
   // delete all children
   while(childNode) delete childNode;
   // remove this node from the tree
   if(GetParentNode() && (GetParentNode()->GetChildNode()==this)) {
      parentNode->childNode=nextNode;
   }
   if(GetPrevNode()) prevNode->nextNode=nextNode;
   if(GetNextNode()) nextNode->prevNode=prevNode;
   delete fAxisList;
   delete fAxisLabelList;
   if(fBinFactorFunction) {
      if(!dynamic_cast<TF1 *>(fBinFactorFunction))
         delete fBinFactorFunction;
   }
}

/********************* setup **************************/

////////////////////////////////////////////////////////////////////////////////
/// Initialize variables for a given number of bins.

void TUnfoldBinning::Initialize(Int_t nBins)
{
   parentNode=0;
   childNode=0;
   nextNode=0;
   prevNode=0;
   fAxisList=new TObjArray();
   fAxisLabelList=new TObjArray();
   fAxisList->SetOwner();
   fAxisLabelList->SetOwner();
   fHasUnderflow=0;
   fHasOverflow=0;
   fDistributionSize=nBins;
   fBinFactorFunction=0;
   fBinFactorConstant=1.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Update fFirstBin and fLastBin members of this node and its children.
///
/// \param[in] startWithRootNode if true, start the update with the root node

Int_t TUnfoldBinning::UpdateFirstLastBin(Bool_t startWithRootNode)
{
   if(startWithRootNode) {
      return GetRootNode()->UpdateFirstLastBin(kFALSE);
   }
   if(GetPrevNode()) {
      // if this is not the first node in a sequence,
      // start with the end bin of the previous node
      fFirstBin=GetPrevNode()->GetEndBin();
   } else if(GetParentNode()) {
      // if this is the first node in a sequence but has a parent,
      // start with the end bin of the parent's distribution
      fFirstBin=GetParentNode()->GetStartBin()+
         GetParentNode()->GetDistributionNumberOfBins();
   } else {
      // if this is the top level node, the first bin number is 1
     fFirstBin=1;
     //  ... unless the top level node is the only node
     //  ... with dimension=1
     //  ... and there are no child nodes
     //  ... and there is an underflow bin
     if((!GetChildNode())&&(GetDistributionDimension()==1)&&
        (fHasUnderflow==1)) {
        fFirstBin=0;
     }
   }
   fLastBin=fFirstBin+fDistributionSize;
   // now update count for all children
   for(TUnfoldBinning *node=childNode;node;node=node->nextNode) {
      fLastBin=node->UpdateFirstLastBin(kFALSE);
   }
   return fLastBin;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new node without axis.
///
/// \param[in] name identifier of the node
/// \param[in] nBins number of unconnected bins (could be zero)
/// \param[in] binNames (optional) names of the bins separated by ';'

TUnfoldBinning::TUnfoldBinning
(const char *name,Int_t nBins,const char *binNames)
   : TNamed(name ? name : "",name ? name : "")
{
   Initialize(nBins);
   if(binNames) {
      TString nameString(binNames);
      delete fAxisLabelList;
      fAxisLabelList=nameString.Tokenize(";");
   }
   UpdateFirstLastBin();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new node containing a distribution with one axis.
///
/// \param[in] axis the axis to represent
/// \param[in] includeUnderflow  true if underflow bin should be included
/// \param[in] includeOverflow true if overflow bin should be included

TUnfoldBinning::TUnfoldBinning
(const TAxis &axis,Int_t includeUnderflow,Int_t includeOverflow)
   : TNamed(axis.GetName(),axis.GetTitle())
{
   Initialize(0);
   AddAxis(axis,includeUnderflow,includeOverflow);
   UpdateFirstLastBin();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new binning node as last last child of this node.
///
/// \param[in] name name of the node
/// \param[in] nBins number of extra bins
/// \param[in] binNames (optional) names of the bins separated by ';'
///
/// this is a shortcut for AddBinning(new TUnfoldBinning(name,nBins,binNames))

TUnfoldBinning *TUnfoldBinning::AddBinning
(const char *name,Int_t nBins,const char *binNames)
{
  return AddBinning(new TUnfoldBinning(name,nBins,binNames));
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TUnfoldBinning as the last child of this node.
///
/// \param[in] binning the new binning to be added
///
/// return value: if succeeded, return "binning"
///               otherwise return 0

TUnfoldBinning *TUnfoldBinning::AddBinning(TUnfoldBinning *binning)
{
  TUnfoldBinning *r=0;
  if(binning->GetParentNode()) {
     Error("AddBinning",
           "binning \"%s\" already has parent \"%s\", can not be added to %s",
     (char *)binning->GetName(),
     (char *)binning->GetParentNode()->GetName(),
     (char *)GetName());
  } else if(binning->GetPrevNode()) {
    Error("AddBinning",
          "binning \"%s\" has previous node \"%s\", can not be added to %s",
     (char *)binning->GetName(),
     (char *)binning->GetPrevNode()->GetName(),
     (char *)GetName());
  } else if(binning->GetNextNode()) {
    Error("AddBinning",
          "binning \"%s\" has next node \"%s\", can not be added to %s",
     (char *)binning->GetName(),
     (char *)binning->GetNextNode()->GetName(),
     (char *)GetName());
  } else {
    r=binning;
    binning->parentNode=this;
    if(childNode) {
      TUnfoldBinning *child=childNode;
      // find last child
      while(child->nextNode) {
   child=child->nextNode;
      }
      // add as last child
      child->nextNode=r;
      r->prevNode=child;
    } else {
      childNode=r;
    }
    UpdateFirstLastBin();
    r=binning;
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an axis with equidistant bins.
///
/// \param[in] name name of the axis
/// \param[in] nBin number of bins
/// \param[in] xMin lower edge of the first bin
/// \param[in] xMax upper edge of the last bin
/// \param[in] hasUnderflow decide whether the axis has an underflow bin
/// \param[in] hasOverflow decide whether the axis has an overflow bin
///
/// returns true if the axis has been added

Bool_t TUnfoldBinning::AddAxis
(const char *name,Int_t nBin,Double_t xMin,Double_t xMax,
 Bool_t hasUnderflow,Bool_t hasOverflow)
{
  Bool_t r=kFALSE;
   if(nBin<=0) {
      Fatal("AddAxis","number of bins %d is not positive",
            nBin);
   } else if((!TMath::Finite(xMin))||(!TMath::Finite(xMax))||
      (xMin>=xMax)) {
      Fatal("AddAxis","xmin=%f required to be smaller than xmax=%f",
            xMin,xMax);
   } else {
     Double_t *binBorders=new Double_t[nBin+1];
     Double_t x=xMin;
     Double_t dx=(xMax-xMin)/nBin;
     for(Int_t i=0;i<=nBin;i++) {
       binBorders[i]=x+i*dx;
     }
     r=AddAxis(name,nBin,binBorders,hasUnderflow,hasOverflow);
     delete [] binBorders;
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an axis to the distribution, using the TAxis as blueprint.
///
/// \param[in] axis blueprint of the axis
/// \param[in] hasUnderflow decide whether the underflow bin should be included
/// \param[in] hasOverflow decide whether the overflow bin should be included
///
/// returns true if the axis has been added
///
/// Note: axis labels are not imported

Bool_t TUnfoldBinning::AddAxis
(const TAxis &axis,Bool_t hasUnderflow,Bool_t hasOverflow)
{
  Int_t nBin=axis.GetNbins();
  Double_t *binBorders=new Double_t[nBin+1];
  for(Int_t i=0;i<nBin;i++) {
    binBorders[i]=axis.GetBinLowEdge(i+1);
  }
  binBorders[nBin]=axis.GetBinUpEdge(nBin);
  Bool_t r=AddAxis(axis.GetTitle(),nBin,binBorders,hasUnderflow,hasOverflow);
  delete [] binBorders;
  return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an axis with the specified bin borders.
///
/// \param[in] name name of the axis
/// \param[in] nBin number of bins
/// \param[in] binBorders array of bin borders, with nBin+1 elements
/// \param[in] hasUnderflow decide whether the axis has an underflow bin
/// \param[in] hasOverflow decide whether the axis has an overflow bin
///
/// returns true if the axis has been added

Bool_t TUnfoldBinning::AddAxis
(const char *name,Int_t nBin,const Double_t *binBorders,
 Bool_t hasUnderflow,Bool_t hasOverflow)
{
  Bool_t r=kFALSE;
  if(HasUnconnectedBins()) {
    Fatal("AddAxis","node already has %d bins without axis",
     GetDistributionNumberOfBins());
  } else if(nBin<=0) {
    Fatal("AddAxis","number of bins %d is not positive",
     nBin);
  } else {
    TVectorD *bins=new TVectorD(nBin+1);
    r=kTRUE;
    for(Int_t i=0;i<=nBin;i++) {
      (*bins)(i)=binBorders[i];
      if(!TMath::Finite((*bins)(i))) {
   Fatal("AddAxis","bin border %d is not finite",i);
   r=kFALSE;
      } else if((i>0)&&((*bins)(i)<=(*bins)(i-1))) {
   Fatal("AddAxis","bins not in order x[%d]=%f <= %f=x[%d]",
         i,(*bins)(i),(*bins)(i-1),i-1);
   r=kFALSE;
      }
    }
    if(r) {
      Int_t axis=fAxisList->GetEntriesFast();
      Int_t bitMask=1<<axis;
      Int_t nBinUO=nBin;
      if(hasUnderflow) {
   fHasUnderflow |= bitMask;
   nBinUO++;
      } else {
   fHasUnderflow &= ~bitMask;
      }
      if(hasOverflow) {
   fHasOverflow |= bitMask;
   nBinUO++;
      } else {
   fHasOverflow &= ~bitMask;
      }
      fAxisList->AddLast(bins);
      fAxisLabelList->AddLast(new TObjString(name));
      if(!fDistributionSize) fDistributionSize=1;
      fDistributionSize *= nBinUO;
      UpdateFirstLastBin();
    }
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Print some information about this binning tree.
///
/// \param[out] out stream to write to
/// \param[in] indent initial indentation (sub-trees have indent+1)
/// \param[in] debug if debug&gt;0 print more information

void TUnfoldBinning::PrintStream(ostream &out,Int_t indent,int debug)
   const {
   for(Int_t i=0;i<indent;i++) out<<"  ";
   out<<"TUnfoldBinning \""<<GetName()<<"\" has ";
   Int_t nBin=GetEndBin()-GetStartBin();
   if(nBin==1) {
      out<<"1 bin";
   } else {
      out<<nBin<<" bins";
   }
   out<<" ["
      <<GetStartBin()<<","<<GetEndBin()<<"] nTH1x="
      <<GetTH1xNumberOfBins()
      <<"\n";
   if(GetDistributionNumberOfBins()) {
      for(Int_t i=0;i<indent;i++) out<<"  ";
      out<<" distribution: "<<GetDistributionNumberOfBins()<<" bins\n";
      if(fAxisList->GetEntriesFast()) {
         /* for(Int_t i=0;i<indent;i++) out<<"  ";
            out<<" axes:\n"; */
          for(Int_t axis=0;axis<GetDistributionDimension();axis++) {
             for(Int_t i=0;i<indent;i++) out<<"  ";
             out<<"  \""
                 <<GetDistributionAxisLabel(axis)
                 <<"\" nbin="<<GetDistributionBinning(axis)->GetNrows()-1;
             if(HasUnderflow(axis)) out<<" plus underflow";
             if(HasOverflow(axis)) out<<" plus overflow";
             out<<"\n";
          }
      } else {
         for(Int_t i=0;i<indent;i++) out<<"  ";
         out<<" no axis\n";
         for(Int_t i=0;i<indent;i++) out<<"  ";
         out<<" names: ";
         for(Int_t ibin=0;(ibin<GetDistributionNumberOfBins())&&
                (ibin<fAxisLabelList->GetEntriesFast());ibin++) {
            if(ibin) out<<";";
            if(GetDistributionAxisLabel(ibin)) {
               out<<GetDistributionAxisLabel(ibin);
            }
         }
         out<<"\n";
      }
      if(debug>0) {
         // print all bins with full name, size, status, user factor
         for(int iBin=GetStartBin();iBin<GetEndBin();iBin++) {
            for(Int_t i=0;i<indent;i++) out<<"  ";
            out<<GetBinName(iBin)
               <<" size="<<GetBinSize(iBin)
               <<" factor="<<GetBinFactor(iBin);
            out<<"\n";
         }
      }
   }
   TUnfoldBinning const *child=GetChildNode();
   if(child) {
      while(child) {
         child->PrintStream(out,indent+1,debug);
         child=child->GetNextNode();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set normalisation factors which are used in calls to GetBinFactor().
///
/// \param[in] normalisation normalisation factor
/// \param[in] binfactor object which defines a factor for each bin
///
/// In the present implementation, <b>binfactor</b> can be a TF1 or a
/// TVectorD. The TF1 is evaluated a the bin centers of the
/// relevant axes. The TVectorD is indexed by the global bin number
/// minus the start bin number of this node.

void TUnfoldBinning::SetBinFactor
(Double_t normalisation,TObject *binfactor) {
   fBinFactorConstant=normalisation;
   if(fBinFactorFunction) {
      if(!dynamic_cast<TF1 *>(fBinFactorFunction))
         delete fBinFactorFunction;
   }
   fBinFactorFunction=binfactor;
}

////////////////////////////////////////////////////////////////////////////////
/// Set normalisation factor and function which are used in calls to GetBinFactor().
///
/// \param[in] normalisation normalisation factor
/// \param[in] userFunc function evaluated at the (multi-dimensional)
/// bin centres

void TUnfoldBinning::SetBinFactorFunction
(Double_t normalisation,TF1 *userFunc) {
   SetBinFactor(normalisation,userFunc);
}

/********************* Navigation **********************/

////////////////////////////////////////////////////////////////////////////////
/// Traverse the tree and return the first node which matches the given name.
///
/// \param[in] name the identifier of the node to find (zero matches
/// the first node)
///
/// returns the node found or zero

TUnfoldBinning const *TUnfoldBinning::FindNode(char const *name) const
{
   TUnfoldBinning const *r=0;
   if((!name)||(!TString(GetName()).CompareTo(name))) {
      r=this;
   }
   for(TUnfoldBinning const *child=GetChildNode();
       (!r) && child;child=child->GetNextNode()) {
      r=child->FindNode(name);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return root node.

TUnfoldBinning *TUnfoldBinning::GetRootNode(void)
{
   TUnfoldBinning *node=this;
   while(node->GetParentNode()) node=node->parentNode;
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Return root node.

TUnfoldBinning const *TUnfoldBinning::GetRootNode(void) const
{
   TUnfoldBinning const *node=this;
   while(node->GetParentNode()) node=node->GetParentNode();
   return node;
}

/********************* Create THxx histograms **********/

////////////////////////////////////////////////////////////////////////////////
/// Construct a title.
///
/// \param[in] histogramName distribution name
/// \param[in] histogramTitle default title
/// \param[in] axisList array indicating which axis of this node is
/// mapped to which histogram axis
///
/// if histogramTitle!=0 this title is used. Otherwise, the title is
/// composed as:
/// histogramName;axisname[axisList[0]];axisname[axisList[1]];...

TString TUnfoldBinning::BuildHistogramTitle
(const char *histogramName,const char *histogramTitle,Int_t const *axisList)
   const
{
   TString r;
   if(histogramTitle) {
      r=histogramTitle;
   } else {
      r=histogramName;
      Int_t iEnd;
      for(iEnd=2;iEnd>0;iEnd--) {
         if(axisList[iEnd]>=0) break;
      }
      for(Int_t i=0;i<=iEnd;i++) {
         r += ";";
         if(axisList[i]<0) {
            r += GetName();
         } else {
            r += GetNonemptyNode()->GetDistributionAxisLabel(axisList[i]);
         }
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a histogram title for a 2D histogram with different
/// binning schemes on x and y axis.
///
/// \param[in] histogramName distribution name
/// \param[in] histogramTitle default title
/// \param[in] xAxis indicates which x-axis name to use
/// \param[in] yAxisBinning binning scheme for y-axis
/// \param[in] yAxis indicates which y-axis name to use
///
/// build a title
///  - input:
///    histogramTitle : if this is non-zero, use that title
///   - otherwise:
///     title=histogramName;x;y
///  - xAxis :
///      - -1 no title for this axis
///      - >=0 use name of the corresponding axis

TString TUnfoldBinning::BuildHistogramTitle2D
(const char *histogramName,const char *histogramTitle,
 Int_t xAxis,const TUnfoldBinning *yAxisBinning,Int_t yAxis) const
{
   TString r;
   if(histogramTitle) {
      r=histogramTitle;
   } else {
      r=histogramName;
      r += ";";
      if(xAxis==-1) {
         r += GetName();
      } else if(xAxis>=0) {
         r += GetNonemptyNode()->GetDistributionAxisLabel(xAxis);
      }
      r+= ";";
      if(yAxis==-1) {
         r += yAxisBinning->GetName();
      } else if(yAxis>=0) {
         r += yAxisBinning->GetNonemptyNode()->GetDistributionAxisLabel(yAxis);
      }

   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of histogram bins required when storing
/// this binning in a one-dimensional histogram.
///
/// \param[in] originalAxisBinning  if true, try to have the histogram
/// axis reflect precisely the relevant axis of the binning scheme
/// \param[in] axisSteering steering to integrate over axis and/or
/// skip underflow and overflow bins
///
/// returns the number of bins of the TH1, where the underflow/overflow
/// are not used, unless the distribution has only one axis and
/// originalAxisBinning=true)
///
/// axisSteering is a string as follows:
/// "axis[options];axis[options];..."
/// where: axis = name or * is an identifier of an axis (* matches all)
/// and: options is any combination of the letters C,U,O  (other
/// letters are ignored).
///
/// The letter C means that the corresponding axis is collapsed into
/// one bin, i.e. one dimension is removed from the counting.
/// The letters U,O remove for the matching axis the underflow.overflow
/// bins from the counting

Int_t TUnfoldBinning::GetTH1xNumberOfBins
(Bool_t originalAxisBinning,const char *axisSteering) const
{
   Int_t axisBins[3],axisList[3];
   GetTHxxBinning(originalAxisBinning ? 1 : 0,axisBins,axisList,
                  axisSteering);
   return axisBins[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Create a THxx histogram capable to hold the bins of this binning
/// node and its children.
///
/// \param[in] histogramName name of the histogram which is created
/// \param[in] originalAxisBinning if true, try to preserve the axis binning
/// \param[out] binMap (default=0) mapping of global bins to histogram bins.
/// if(binMap==0), no binMap is created
/// \param[in] histogramTitle (default=0) title of the histogram. If zero, a title
/// is selected automatically
/// \param[in] axisSteering (default=0) steer the handling of underflow/overflow
/// and projections
///
/// returns a new histogram (TH1D, TH2D or TH3D)
///
/// if the parameter <b>originalAxisBinning</b> parameter is true, the
/// resulting histogram has bin widths and histogram dimension (TH1D,
/// TH2D, TH3D) in parallel to this binning node, if possible.
///
/// The <b>binMap</b> is an array which translates global bin numbers to bin
/// numbers in the histogram returned by this method. The global bin
/// numbers correspond to the bin numbers in a histogram created by
/// calling GetRootNode()->CreateHistogram(name,false,0,0,0)
///
/// The <b>axisSteering</b> is a string to steer whether underflow and
/// overflow bins are included in the bin map. Furthermore, it is
/// possible to "collapse" axes, such that their content is summed and
/// the axis does not show up in the created histogram.
///
/// The string looks like this: "axis[options];axis[options];..." where
///
///   - axis is the name of an axis or equal to *, the latter matches
/// all axes
///   - options is a combination of characters chosen from
/// OUC0123456789
///
///   - if O is included, the overflow bin of that axis is discarded
///   - if U is included, the underflow bin of that axis is discarded
///   - if C is included, the bins on that axes are collapsed,
/// i.e. the corresponding histogram axis is not present in the output.
/// The corresponding bin contents are added up
/// (projected onto the remaining axes). Using the characters O and U
/// one can decide to exclude underflow or overflow from the
/// projection. Using a selection of the characters 0123456789 one can
/// restrict the sum further to only include the corresponding
/// bins. In this counting, the first non-underflow bin corresponds to
/// the character 0. This obviously only works for up to ten
/// bins.

TH1 *TUnfoldBinning::CreateHistogram
(const char *histogramName,Bool_t originalAxisBinning,Int_t **binMap,
 const char *histogramTitle,const char *axisSteering) const
{
   Int_t nBin[3],axisList[3];
   Int_t nDim=GetTHxxBinning(originalAxisBinning ? 3 : 0,nBin,axisList,
                             axisSteering);
   const TUnfoldBinning *neNode=GetNonemptyNode();
   TString title=BuildHistogramTitle(histogramName,histogramTitle,axisList);
   TH1 *r=0;
   if(nDim>0) {
      const TVectorD *axisBinsX=
         neNode->GetDistributionBinning(axisList[0]);
      if(nDim>1) {
         const TVectorD *axisBinsY=
            neNode->GetDistributionBinning(axisList[1]);
         if(nDim>2) {
            const TVectorD *axisBinsZ=
               neNode->GetDistributionBinning(axisList[2]);
            r=new TH3D(histogramName,title,
                       nBin[0],axisBinsX->GetMatrixArray(),
                       nBin[1],axisBinsY->GetMatrixArray(),
                       nBin[2],axisBinsZ->GetMatrixArray());
         } else {
            r=new TH2D(histogramName,title,
                       nBin[0],axisBinsX->GetMatrixArray(),
                       nBin[1],axisBinsY->GetMatrixArray());
         }
      } else {
         r=new TH1D(histogramName,title,nBin[0],axisBinsX->GetMatrixArray());
      }
   } else {
      if(originalAxisBinning) {
         Warning("CreateHistogram",
       "Original binning can not be represented as THxx");
      }
      r=new TH1D(histogramName,title,nBin[0],0.5,nBin[0]+0.5);
      nDim=0;
   }
   if(binMap) {
      *binMap=CreateBinMap(r,nDim,axisList,axisSteering);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TH2D histogram capable to hold a covariance matrix.
///
///
/// \param[in] histogramName name of the histogram which is created
/// \param[in] originalAxisBinning if true, try to preserve the axis binning
/// \param[out] binMap (default=0) mapping of global bins to histogram bins.
/// if(binMap==0), no binMap is created
/// \param[in] histogramTitle (default=0) title of the histogram. If zero, a title
/// is selected automatically
/// \param[in] axisSteering (default=0) steer the handling of underflow/overflow
/// and projections
///
/// returns a new TH2D. The options are described in greater detail
/// with the CreateHistogram() method.

TH2D *TUnfoldBinning::CreateErrorMatrixHistogram
(const char *histogramName,Bool_t originalAxisBinning,Int_t **binMap,
 const char *histogramTitle,const char *axisSteering) const
{
   Int_t nBin[3],axisList[3];
   Int_t nDim=GetTHxxBinning(originalAxisBinning ? 1 : 0,nBin,axisList,
                             axisSteering);
   TString title=BuildHistogramTitle(histogramName,histogramTitle,axisList);
   TH2D *r=0;
   if(nDim==1) {
      const TVectorD *axisBinsX=(TVectorD const *)
         GetNonemptyNode()->fAxisList->At(axisList[0]);
      r=new TH2D(histogramName,title,nBin[0],axisBinsX->GetMatrixArray(),
                 nBin[0],axisBinsX->GetMatrixArray());
   } else {
      if(originalAxisBinning) {
         Info("CreateErrorMatrixHistogram",
              "Original binning can not be represented on one axis");
      }
      r=new TH2D(histogramName,title,nBin[0],0.5,nBin[0]+0.5,
                 nBin[0],0.5,nBin[0]+0.5);
      nDim=0;
   }
   if(binMap) {
      *binMap=CreateBinMap(r,nDim,axisList,axisSteering);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TH2D histogram capable to hold the bins of the two
/// input binning schemes on the x and y axes, respectively.
///
/// \param[in] xAxis binning scheme for the x axis
/// \param[in] yAxis binning scheme for the y axis
/// \param[in] histogramName name of the histogram which is created
/// \param[in] originalXAxisBinning preserve x-axis bin widths if possible
/// \param[in] originalXAxisBinning preserve y-axis bin widths if possible
/// \param[in] histogramTitle if is non-zero, it is taken as histogram title
///                     otherwise, the title is created automatically
///
/// returns a new TH2D.

TH2D *TUnfoldBinning::CreateHistogramOfMigrations
(TUnfoldBinning const *xAxis,TUnfoldBinning const *yAxis,
 char const *histogramName,Bool_t originalXAxisBinning,
 Bool_t originalYAxisBinning,char const *histogramTitle)
{
   Int_t nBinX[3],axisListX[3];
   Int_t nDimX=
      xAxis->GetTHxxBinning(originalXAxisBinning ? 1 : 0,nBinX,axisListX,0);
   const TUnfoldBinning *neNodeX=xAxis->GetNonemptyNode();
   Int_t nBinY[3],axisListY[3];
   Int_t nDimY=
      yAxis->GetTHxxBinning(originalYAxisBinning ? 1 : 0,nBinY,axisListY,0);
   const TUnfoldBinning *neNodeY=yAxis->GetNonemptyNode();
   TString title=xAxis->BuildHistogramTitle2D
      (histogramName,histogramTitle,axisListX[0],yAxis,axisListY[0]);
   if(nDimX==1) {
      const TVectorD *axisBinsX=(TVectorD const *)
            neNodeX->fAxisList->At(axisListX[0]);
      if(nDimY==1) {
         const TVectorD *axisBinsY=(TVectorD const *)
            neNodeY->fAxisList->At(axisListY[0]);
         return new TH2D(histogramName,title,
                         nBinX[0],axisBinsX->GetMatrixArray(),
                         nBinY[0],axisBinsY->GetMatrixArray());
      } else {
         return new TH2D(histogramName,title,
                         nBinX[0],axisBinsX->GetMatrixArray(),
                         nBinY[0],0.5,0.5+nBinY[0]);
      }
   } else {
      if(nDimY==1) {
         const TVectorD *axisBinsY=(TVectorD const *)
             neNodeY->fAxisList->At(axisListY[0]);
         return new TH2D(histogramName,title,
                         nBinX[0],0.5,0.5+nBinX[0],
                         nBinY[0],axisBinsY->GetMatrixArray());
      } else {
         return new TH2D(histogramName,title,
                         nBinX[0],0.5,0.5+nBinX[0],
                         nBinY[0],0.5,0.5+nBinY[0]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate properties of a THxx histogram to store this binning.
///
/// \param[in] maxDim  maximum dimension of the THxx (0 or 1..3)
/// maxDim==0 is used to indicate that the histogram should be
///    dimensional with all bins mapped on one axis,
///    bin centers equal to bin numbers
/// \param[in] axisSteering see method CreateHistogram()
/// \param[out] axisBins[3] number of bins on the THxx axes
/// \param[out] axisList[3] TUnfoldBinning axis number corresponding
/// to the THxx axis
///
/// returns 1-3 dimension of THxx or 0 for 1-dim THxx with equidistant bins

Int_t TUnfoldBinning::GetTHxxBinning
(Int_t maxDim,Int_t *axisBins,Int_t *axisList,
 const char *axisSteering) const
{
   for(Int_t i=0;i<3;i++) {
      axisBins[i]=0;
      axisList[i]=-1;
   }
   const TUnfoldBinning *theNode=GetNonemptyNode();
   if(theNode) {
      Int_t r=theNode->GetTHxxBinningSingleNode
         (maxDim,axisBins,axisList,axisSteering);
      return r;
   } else {
      axisBins[0]=GetTHxxBinsRecursive(axisSteering);
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find a node which has non-empty distributions
/// if there is none or if there are many, return zero.

const TUnfoldBinning *TUnfoldBinning::GetNonemptyNode(void) const
{
   const TUnfoldBinning *r=GetDistributionNumberOfBins()>0 ? this : 0;
   for(TUnfoldBinning const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      const TUnfoldBinning *c=child->GetNonemptyNode();
      if(!r) {
         // new candidate found
         r=c;
      } else {
         if(c) {
            // multiple nodes found
            r=0;
            break;
         }
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the properties of a histogram capable to hold the distribution
/// attached to this node.
///
/// \param[in] maxDim maximum dimension of the THxx (0 or 1..3)
///              maxDim==0 is used to indicate that the histogram should
///              1-dimensional with all bins mapped on one axis
/// \param[out] axisBins[3] number of bins on the THxx axes
/// \param[out] axisList[3] TUnfoldBinning axis numbers
///              corresponding to the THxx axis
/// \param[in] axisSteering  see method CreateHistogram()
/// and projection
///
/// returns 1-3 dimension of THxx or use 1-dim THxx, binning structure
/// is not preserved

Int_t TUnfoldBinning::GetTHxxBinningSingleNode
(Int_t maxDim,Int_t *axisBins,Int_t *axisList,const char *axisSteering) const
{
   // decode axisSteering
   //   isOptionGiven[0] ('C'): bit vector which axes to collapse
   //   isOptionGiven[1] ('U'): bit vector to discard underflow bins
   //   isOptionGiven[2] ('O'): bit vector to discard overflow bins
   Int_t isOptionGiven[3];
   DecodeAxisSteering(axisSteering,"CUO",isOptionGiven);
   // count number of axes after projecting
   Int_t numDimension=GetDistributionDimension();
   Int_t r=0;
   for(Int_t i=0;i<numDimension;i++) {
      if(isOptionGiven[0] & (1<<i)) continue;
      r++;
   }
   if((r>0)&&(r<=maxDim)) {
      // 0<r<=maxDim
      //
      // -> preserve the original binning
      //    axisList[] and axisBins[] are overwritten
      r=0;
      for(Int_t i=0;i<numDimension;i++) {
         if(isOptionGiven[0] & (1<<i)) continue;
         axisList[r]=i;
         axisBins[r]=GetDistributionBinning(i)->GetNrows()-1;
         r++;
      }
   } else {
      // map everything on one axis
      //  axisBins[0] is the number of bins
      if(HasUnconnectedBins() || (GetDistributionNumberOfBins()<=0)) {
         axisBins[0] = GetDistributionNumberOfBins();
      } else {
         Int_t nBin=1;
         for(Int_t i=0;i<numDimension;i++) {
            Int_t mask=(1<<i);
            if(isOptionGiven[0] & mask) continue;
            Int_t nBinI=GetDistributionBinning(i)->GetNrows()-1;
            if((fHasUnderflow & mask)&& !(isOptionGiven[1] & mask)) nBinI++;
            if((fHasOverflow & mask)&& !(isOptionGiven[2] & mask)) nBinI++;
            nBin *= nBinI;
         }
         axisBins[0] = nBin;
      }
      r=0;
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate number of bins required to store this binning with the
/// given axisSteering.
///
/// \param[in] axisSteering see method CreateHistogram()
///
/// returns the number of bins

Int_t TUnfoldBinning::GetTHxxBinsRecursive(const char *axisSteering) const
{

   Int_t r=0;
   for(TUnfoldBinning const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      r +=child->GetTHxxBinsRecursive(axisSteering);
   }
   // here: process distribution of this node
   Int_t axisBins[3],axisList[3];
   GetTHxxBinningSingleNode(0,axisBins,axisList,axisSteering);
   r += axisBins[0];
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an empty bin map, useful together with the getter methods of
/// class TUnfold and TUnfoldSys.
///
/// returns: a new Int array of the proper size, all elements set to -1

Int_t *TUnfoldBinning::CreateEmptyBinMap(void) const {
   // create empty bin map which can be manipulated by
   //  MapGlobalBin()
   Int_t nMax=GetRootNode()->GetEndBin()+1;
   Int_t *r=new Int_t[nMax];
   for(Int_t i=0;i<nMax;i++) {
         r[i]=-1;
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Set one entry in a bin map.
///
/// \param[out] binMap to be used with TUnfoldSys::GetOutput() etc
/// \param[in] globalBin source bin, global bin number in this binning scheme
/// \param[in] destBin destination bin in the output histogram

void TUnfoldBinning::SetBinMapEntry
(Int_t *binMap,Int_t globalBin,Int_t destBin) const {
   Int_t nMax=GetRootNode()->GetEndBin()+1;
   if((globalBin<0)||(globalBin>=nMax)) {
      Error("SetBinMapEntry","global bin number %d outside range (max=%d)",
            globalBin,nMax);
   } else {
      binMap[globalBin]=destBin;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map all global bins referenced by this node to the one-dimensional
/// histogram destHist, starting with bin firstBinX
///
/// \param[out] binMap to be used with TUnfoldSys::GetOutput() etc
/// \param[in] axisSteering steering for underflow/overflow/projections
/// \param[in] firstBinX first bin of destination histogram to be filled
///
/// returns: highest bin number in destination histogram plus 1
/// The parameter <b>axisSteering</b> is explained with the
/// method CreateHistogram()

Int_t TUnfoldBinning::FillBinMap1D
(Int_t *binMap,const char *axisSteering,Int_t firstBinX) const {
   Int_t r=firstBinX;
   Int_t axisBins[3],axisList[3];
   Int_t nDim=GetTHxxBinningSingleNode(3,axisBins,axisList,axisSteering);
   if((nDim==1)|| !GetDistributionDimension()) {
      r+=FillBinMapSingleNode(0,r,0,0,axisSteering,binMap);
   } else {
      Error("FillBinMap1D","distribution %s with steering=%s is not 1D",
            (char *)GetName(),axisSteering);
   }
   for(TUnfoldBinning const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      r =child->FillBinMap1D(binMap,axisSteering,r);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Create mapping from global bin number to a histogram for this node.
///
/// \param[in] hist destination histogram
/// \param[in] nDim target dimension
/// \param[in] axisList map axes in the binning scheme to histogram axes
/// \param[in] axisSteering steering for underflow/overflow/projections
///
/// The <b>axisSteering</b> is explained with the method CreateHistogram()
/// create mapping from global bin number to a histogram for this node
/// global bins are the bins of the root node binning scheme
/// when projecting them on a TH1 histogram "hRootNode" without special
/// axis steering and without attempting to preserve the axis binning
///
/// The bin map is an array of size hRootNode->GetNbinsX()+2
/// For each bin of the "hRootNode" histogram it holds the target bin in
/// "hist" or the number -1 if the corresponding "hRootNode" bin is not
/// represented in "hist"
///
/// input
///  - hist : the histogram (to calculate root bin numbers)
///  - nDim : target dimension of the TUnfoldBinning
///          if(nDim==0) all bins are mapped linearly
///  - axisSteering:
///       "pattern1;pattern2;...;patternN"
///       patternI = axis[mode]
///       axis = name or *
///       mode = C|U|O
///       - C: collapse axis into one bin
///       - U: discard underflow bin
///       - O: discard overflow bin
///
///  - input used only if nDim>0:
///    axisList : for each THxx axis give the TUnfoldBinning axis number
///
///  - return value:
///    an new array which holds the bin mapping
///      - r[0] : to which THxx bin to map global bin number 0
///      - r[1] : to which THxx bin to map global bin number 1
///      ...
///      - r[nmax]
///        where nmax=GetRootNode()->GetEndBin()+1

Int_t *TUnfoldBinning::CreateBinMap
(const TH1 *hist,Int_t nDim,const Int_t *axisList,const char *axisSteering)
   const
{
   Int_t *r=CreateEmptyBinMap();
   Int_t startBin=GetRootNode()->GetStartBin();
   if(nDim>0) {
     const TUnfoldBinning *nonemptyNode=GetNonemptyNode();
     if(nonemptyNode) {
       nonemptyNode->
          FillBinMapSingleNode(hist,startBin,nDim,axisList,axisSteering,r);
     } else {
       Fatal("CreateBinMap","called with nDim=%d but GetNonemptyNode()=0",
        nDim);
     }
   } else {
     FillBinMapRecursive(startBin,axisSteering,r);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively fill bin map.
///
/// \param[in] startBin first histogram bin
/// \param[in] axisSteering see CreateHistogram() method
/// \param[out]  binMap the bin mapping which is to be filled
///
/// the positions
///      - binMap[GetStartBin()]...binMap[GetEndBin()-1]
/// are filled

Int_t TUnfoldBinning::FillBinMapRecursive
(Int_t startBin,const char *axisSteering,Int_t *binMap) const
{
   Int_t nbin=0;
   nbin = FillBinMapSingleNode(0,startBin,0,0,axisSteering,binMap);
   for(TUnfoldBinning const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      nbin += child->FillBinMapRecursive(startBin+nbin,axisSteering,binMap);
   }
   return nbin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill bin map for a single node.
///
/// \param[in] hist the histogram representing this node (used if nDim>0)
/// \param[in] startBin start bin in the bin map
/// \param[in] nDim number of dimensions to resolve
/// \param[in] axisList[3] TUnfoldBinning axis numbers corresponding
/// to the axes of <b>hist</b>
/// \param[in] axisSteering see documentation of CreateHistogram()
/// \param[out] binMap the bin map to fill
///
/// returns the number of bins mapped.
///
/// The result depends on the parameter <b>nDim</b> as follows
///
///   -  nDim==0: bins are mapped in linear order, ignore hist and
/// axisList
///   -  nDim==hist->GetDimension():
/// bins are mapped to "hist" bin numbers
/// the corresponding TUnfoldBinning axes are taken from
/// axisList[]
///   - nDim=1 and hist->GetDimension()>1:
/// bins are mapped to the x-axis of "hist"
/// the corresponding TUnfoldBinning axis is taken from
/// axisList[0]

Int_t TUnfoldBinning::FillBinMapSingleNode
(const TH1 *hist,Int_t startBin,Int_t nDim,const Int_t *axisList,
 const char *axisSteering,Int_t *binMap) const
{
   // first, decode axisSteering
   //   isOptionGiven[0] ('C'): bit vector which axes to collapse
   //   isOptionGiven[1] ('U'): bit vector to discard underflow bins
   //   isOptionGiven[2] ('O'): bit vector to discard overflow bins
   Int_t isOptionGiven[3+10];
   DecodeAxisSteering(axisSteering,"CUO0123456789",isOptionGiven);
   Int_t haveSelectedBin=0;
   for(Int_t i=3;i<3+10;i++) {
      haveSelectedBin |= isOptionGiven[i];
   }

   Int_t axisBins[MAXDIM];
   Int_t dimension=GetDistributionDimension();
   Int_t axisNbin[MAXDIM];
   for(Int_t i=0;i<dimension;i++) {
      const TVectorD *binning=GetDistributionBinning(i);
      axisNbin[i]=binning->GetNrows()-1;
   };
   for(Int_t i=0;i<GetDistributionNumberOfBins();i++) {
      Int_t globalBin=GetStartBin()+i;
      const TUnfoldBinning *dest=ToAxisBins(globalBin,axisBins);
      if(dest!=this) {
         if(!dest) {
            Fatal("FillBinMapSingleNode",
                  "bin %d outside binning scheme",
                  globalBin);
         } else {
            Fatal("FillBinMapSingleNode",
                  "bin %d located in %s %d-%d rather than %s %d=%d",
                  i,(const char *)dest->GetName(),
                  dest->GetStartBin(),dest->GetEndBin(),
                  (const char *)GetName(),GetStartBin(),GetEndBin());
         }
      }
      // check whether this bin has to be skipped
      Bool_t skip=kFALSE;
      for(Int_t axis=0;axis<dimension;axis++) {
         Int_t mask=(1<<axis);
         // underflow/overflow excluded by steering
         if(((axisBins[axis]<0)&&(isOptionGiven[1] & mask))||
            ((axisBins[axis]>=axisNbin[axis])&&(isOptionGiven[2] & mask)))
            skip=kTRUE;
         // only certain bins selected by steering
         if((axisBins[axis]>=0)&&(axisBins[axis]<axisNbin[axis])&&
            (haveSelectedBin & mask)) {
            if(!(isOptionGiven[3+axisBins[axis]] & mask)) skip=kTRUE;
         }
      }
      if(skip) {
         continue;
      }

      if(nDim>0) {
         // get bin number from THxx function(s)
         if(nDim==hist->GetDimension()) {
            Int_t ibin[3];
            ibin[0]=ibin[1]=ibin[2]=0;
            for(Int_t hdim=0;hdim<nDim;hdim++) {
               Int_t axis=axisList[hdim];
               ibin[hdim]=axisBins[axis]+1;
            }
            binMap[globalBin]=hist->GetBin(ibin[0],ibin[1],ibin[2]);
         } else if(nDim==1) {
            // histogram has more dimensions than the binning scheme
            // and the binning scheme has one axis only
            //
            // special case: error histogram is 2-d
            //   create nor error if ndim==1 && hist->GetDimension()==2
            if((nDim!=1)||( hist->GetDimension()!=2)) {
               // -> use the first valid axis only
               Error("FillBinMapSingleNode","inconsistent dimensions %d %d",nDim,
                     hist->GetDimension());
            }
            for(Int_t ii=0;ii<hist->GetDimension();ii++) {
               if(axisList[ii]>=0) {
                  binMap[globalBin]=axisBins[axisList[ii]]+1;
                  break;
               }
            }
         } else {
            Fatal("FillBinMapSingleNode","inconsistent dimensions %d %d",nDim,
                  hist->GetDimension());
         }
      } else {
         // order all bins in sequence
         // calculation in parallel to ToGlobalBin()
         // but take care of
         //   startBin,collapseAxis,discardeUnderflow,discardeOverflow
         if(dimension>0) {
            Int_t r=0;
            for(Int_t axis=dimension-1;axis>=0;axis--) {
               Int_t mask=(1<<axis);
               if(isOptionGiven[0] & mask) {
                  // bins on this axis are integrated over
                  continue;
               }
               Int_t iBin=axisBins[axis];
               Int_t nMax=axisNbin[axis];
               if((fHasUnderflow & ~isOptionGiven[1]) & mask) {
                  nMax +=1;
                  iBin +=1;
               }
               if((fHasOverflow & ~isOptionGiven[2]) & mask) {
                  nMax += 1;
               }
               r = r*nMax +iBin;
            }
            binMap[globalBin] = startBin + r;
         } else {
      binMap[globalBin] = startBin + axisBins[0];
         }
      }
   }
   Int_t nbin;
   if(dimension>0) {
     nbin=1;
     for(Int_t axis=dimension-1;axis>=0;axis--) {
       Int_t mask=(1<<axis);
       if(isOptionGiven[0] & mask) {
    // bins on this axis are integrated over
    continue;
       }
       Int_t nMax=axisNbin[axis];
       if((fHasUnderflow & ~isOptionGiven[1]) & mask) {
    nMax +=1;
       }
       if((fHasOverflow & ~isOptionGiven[2]) & mask) {
    nMax += 1;
       }
       nbin = nbin*nMax;
     }
   } else {
     nbin=GetDistributionNumberOfBins();
   }
   return nbin;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract a distribution from the given set of global bins.
///
/// input:
///   - histogramName : name of the histogram which ic created
///   - globalBins : histogram with all bins
///   - globalBinsEmatrix : corresponding error matrix
///                 if this pointer is zero, only diagonal errors
///                 are considered
///   - originalAxisBinning :  extract  histogram with proper binning
///                          (if possible)
///   - axisSteering
///      - "pattern1;pattern2;...;patternN"
///      - patternI = axis[mode]
///      - axis = name or *
///      - mode = C|U|O
///         - C: collapse axis into one bin
///         - U: discard underflow bin
///         - O: discard overflow bin

TH1 *TUnfoldBinning::ExtractHistogram
(const char *histogramName,const TH1 *globalBins,
 const TH2 *globalBinsEmatrix,Bool_t originalAxisBinning,
 const char *axisSteering) const
{
   Int_t *binMap=0;
   TH1 *r=CreateHistogram(histogramName,originalAxisBinning,&binMap,0,
                          axisSteering);
   if(!r) return 0;
   TUnfoldBinning const *root=GetRootNode();
   Int_t nMax=-1;
   for(Int_t iSrc=root->GetStartBin();iSrc<root->GetEndBin();iSrc++) {
      if(binMap[iSrc]>nMax) nMax=binMap[iSrc];
   }
   if(nMax<0) {
      delete r;
      r=0;
   } else {
      TVectorD eSquared(nMax+1);
      for(Int_t iSrc=root->GetStartBin();iSrc<root->GetEndBin();iSrc++) {
         Int_t iDest=binMap[iSrc];
         if(iDest>=0) {
            Double_t c=r->GetBinContent(iDest);
            r->SetBinContent(iDest,c+globalBins->GetBinContent(iSrc));
            if(!globalBinsEmatrix) {
               eSquared(iDest)+=TMath::Power(globalBins->GetBinError(iSrc),2.);
            } else {
               for(Int_t jSrc=root->GetStartBin();jSrc<root->GetEndBin();
                   jSrc++) {
                  if(binMap[jSrc]==iDest) {
                     eSquared(iDest) +=
                        TMath::Power(globalBins->GetBinError(jSrc),2.);
                  }
               }
            }
         }
      }
      for(Int_t i=0;i<nMax;i++) {
         Double_t e2=eSquared(i);
         if(e2>0.0) {
            r->SetBinError(i,TMath::Sqrt(e2));
         }
      }
   }
   delete[] binMap;
   return r;
}

/********************* Calculate global bin number ******/

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a one-dimensional distribution.
///
/// \param[in] x coordinate
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme

Int_t TUnfoldBinning::GetGlobalBinNumber(Double_t x) const
{
   if(GetDistributionDimension()!=1) {
      Fatal("GetBinNumber",
            "called with 1 argument for %d dimensional distribution",
            GetDistributionDimension());
   }
   return GetGlobalBinNumber(&x);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a two-dimensional distribution.
///
/// \param[in] x coordinate on first axis
/// \param[in] y coordinate on second axis
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme

Int_t TUnfoldBinning::GetGlobalBinNumber(Double_t x,Double_t y) const
{
   if(GetDistributionDimension()!=2) {
      Fatal("GetBinNumber",
            "called with 2 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[2];
   xx[0]=x;
   xx[1]=y;
   return GetGlobalBinNumber(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a three-dimensional distribution.
///
/// \param[in] x coordinate on first axis
/// \param[in] y coordinate on second axis
/// \param[in] z coordinate on third axis
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme
///
/// locate bin on a three-dimensional distribution
///  - input:
///    x,y,z: coordinates to locate

Int_t TUnfoldBinning::GetGlobalBinNumber
(Double_t x,Double_t y,Double_t z) const
{
   if(GetDistributionDimension()!=3) {
      Fatal("GetBinNumber",
            "called with 3 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   return GetGlobalBinNumber(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a four-dimensional distribution.
///
/// \param[in] x0 coordinate on first axis
/// \param[in] x1 coordinate on second axis
/// \param[in] x2 coordinate on third axis
/// \param[in] x3 coordinate on fourth axis
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme
///
/// locate bin on a four-dimensional distribution
///   - input:
///     x0,x1,x2,x3: coordinates to locate

Int_t TUnfoldBinning::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3) const
{
   if(GetDistributionDimension()!=4) {
      Fatal("GetBinNumber",
            "called with 4 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[4];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   return GetGlobalBinNumber(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a five-dimensional distribution.
///
/// \param[in] x0 coordinate on first axis
/// \param[in] x1 coordinate on second axis
/// \param[in] x2 coordinate on third axis
/// \param[in] x3 coordinate on fourth axis
/// \param[in] x4 coordinate on fifth axis
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme
///
/// locate bin on a five-dimensional distribution
///  - input:
///    x0,x1,x2,x3,x4: coordinates to locate

Int_t TUnfoldBinning::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3,Double_t x4) const
{
   if(GetDistributionDimension()!=5) {
      Fatal("GetBinNumber",
            "called with 5 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[5];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   xx[4]=x4;
   return GetGlobalBinNumber(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// Locate a bin in a six-dimensional distribution.
///
/// \param[in] x0 coordinate on first axis
/// \param[in] x1 coordinate on second axis
/// \param[in] x2 coordinate on third axis
/// \param[in] x3 coordinate on fourth axis
/// \param[in] x4 coordinate on fifth axis
/// \param[in] x5 coordinate on sixth axis
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme
///
/// locate bin on a five-dimensional distribution
///  - input:
///    x0,x1,x2,x3,x4,x5: coordinates to locate

Int_t TUnfoldBinning::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3,Double_t x4,Double_t x5) const
{
   if(GetDistributionDimension()!=6) {
      Fatal("GetBinNumber",
            "called with 6 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[6];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   xx[4]=x4;
   xx[5]=x5;
   return GetGlobalBinNumber(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// locate a bin in an N-dimensional distribution
///
/// \param[in] x array of coordinates
/// \param[out] isBelow pointer to an integer (bit vector) to indicate
/// coordinates which do not fit in the binning scheme
/// \param[out] isAbove pointer to an integer (bit vector) to indicate
/// coordinates which do not fit in the binning scheme
///
/// returns the global bin number within the distribution attached to
/// this node. The global bin number is valid for the root node of the
/// binning scheme. If some coordinates do not fit, zero is returned.
/// The integers pointed to by isBelow and isAbove are set to zero.
/// However, if coordinate i is below
/// the lowest bin border and there is no underflow bin, the bin i is
/// set in (*isBelow). Overflows are handled in a similar manner with
/// (*isAbove).
///
/// If a coordinate is NaN, the result is undefined for TUnfold
/// Version<17.6. As of version 17.6, NaN is expected to end up in the
/// underflow or by setting the corresponding bit in (*isBelow).
///
/// locate bin on a n-dimensional distribution
///  - input
///    x[]: coordinates to locate
///  - output:
///    isBelow,isAbove: bit vectors,
///       indicating which cut on which axis failed

Int_t TUnfoldBinning::GetGlobalBinNumber
(const Double_t *x,Int_t *isBelow,Int_t *isAbove) const
{
   if(!GetDistributionDimension()) {
      Fatal("GetBinNumber",
            "no axes are defined for node %s",
            (char const *)GetName());
   }
   Int_t iAxisBins[MAXDIM] = {0};
   for(Int_t dim=0;dim<GetDistributionDimension();dim++) {
      TVectorD const *bins=(TVectorD const *) fAxisList->At(dim);
      Int_t i0=0;
      Int_t i1=bins->GetNrows()-1;
      Int_t iBin= 0;
      if(!(x[dim]>=(*bins)[i0])) {
         // underflow or NaN
         iBin += i0-1;
      } else if(!(x[dim]<(*bins)[i1])) {
         // overflow
         iBin += i1;
      } else {
         while(i1-i0>1) {
            Int_t i2=(i0+i1)/2;
            if(x[dim]<(*bins)[i2]) {
               i1=i2;
            } else {
               i0=i2;
            }
         }
         iBin += i0;
      }
      iAxisBins[dim]=iBin;
   }
   Int_t r=ToGlobalBin(iAxisBins,isBelow,isAbove);
   if(r<0) r=0;
   return r;
}

/********************* access by global bin number ******/

////////////////////////////////////////////////////////////////////////////////
/// Get the name of a bin.
///
/// \param[in] iBin global bin number
///
/// returns a string describing the bin
///
/// Get the name of a bin in the given tree
///    iBin: bin number

TString TUnfoldBinning::GetBinName(Int_t iBin) const
{
   Int_t axisBins[MAXDIM];
   TString r=TString::Format("#%d",iBin);
   TUnfoldBinning const *distribution=ToAxisBins(iBin,axisBins);
   if(distribution) {
      r +=" (";
      r += distribution->GetName();
      Int_t dimension=distribution->GetDistributionDimension();
      if(dimension>0) {
         TString axisString;
         for(Int_t axis=0;axis<dimension;axis++) {
            TString thisAxisString=
               distribution->GetDistributionAxisLabel(axis);
            TVectorD const *bins=distribution->GetDistributionBinning(axis);
            Int_t i=axisBins[axis];
            if(i<0) thisAxisString += "[ufl]";
            else if(i>=bins->GetNrows()-1) thisAxisString += "[ofl]";
            else {
               thisAxisString +=
                  TString::Format("[%.3g,%.3g]",(*bins)[i],(*bins)[i+1]);
            }
            axisString = ":"+thisAxisString+axisString;
         }
         r += axisString;
      } else {
         // extra bins
         Int_t i=axisBins[0];
         if((i>=0)&&(i<distribution->fAxisLabelList->GetEntriesFast())) {
            r += distribution->GetDistributionAxisLabel(i);
         } else {
            r += TString::Format(" %d",i);
         }
      }
      r +=")";
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get N-dimensional bin size.
///
/// \param[in] iBin global bin number
///
///    includeUO : include underflow/overflow bins or not

Double_t TUnfoldBinning::GetBinSize(Int_t iBin) const
{
   Int_t axisBins[MAXDIM];
   TUnfoldBinning const *distribution=ToAxisBins(iBin,axisBins);
   Double_t r=0.0;
   if(distribution) {
      if(distribution->GetDistributionDimension()>0) r=1.0;
      for(Int_t axis=0;axis<distribution->GetDistributionDimension();axis++) {
         TVectorD const *bins=distribution->GetDistributionBinning(axis);
         Int_t pos=axisBins[axis];
         if(pos<0) {
       r *= distribution->GetDistributionUnderflowBinWidth(axis);
         } else if(pos>=bins->GetNrows()-1) {
       r *= distribution->GetDistributionOverflowBinWidth(axis);
         } else {
            r *= (*bins)(pos+1)-(*bins)(pos);
         }
         if(r<=0.) break;
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether there is only a global scaling factor for  this node.

Bool_t TUnfoldBinning::IsBinFactorGlobal(void) const {
   return fBinFactorFunction ? kFALSE : kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return global scaling factor for  this node.

Double_t TUnfoldBinning::GetGlobalFactor(void) const {
   return fBinFactorConstant;
}

////////////////////////////////////////////////////////////////////////////////
/// Return scaling factor for the given global bin number.
///
/// \param[in] iBin global bin number
///
/// returns the scaling factor for this bin.
/// The scaling factors can be set using the method SetBinFactorFunction()
///
/// return user factor for a bin
///    iBin : global bin number

Double_t TUnfoldBinning::GetBinFactor(Int_t iBin) const
{
   Int_t axisBins[MAXDIM];
   TUnfoldBinning const *distribution=ToAxisBins(iBin,axisBins);
   Double_t r=distribution->fBinFactorConstant;
   if((r!=0.0) && distribution->fBinFactorFunction) {
      TF1 *function=dynamic_cast<TF1 *>(distribution->fBinFactorFunction);
      if(function) {
         Double_t x[MAXDIM];
         Int_t dimension=distribution->GetDistributionDimension();
         if(dimension>0) {
            for(Int_t  axis=0;axis<dimension;axis++) {
               x[axis]=distribution->GetDistributionBinCenter
                  (axis,axisBins[axis]);
            }
            r *= function->EvalPar(x,function->GetParameters());
         } else {
            x[0]=axisBins[0];
            r *= function->Eval(x[0]);
         }
      } else {
         TVectorD *vect=dynamic_cast<TVectorD *>
            (distribution->fBinFactorFunction);
         if(vect) {
            r=(*vect)[iBin-GetStartBin()];
         } else {
            Error("GetBinFactor",
                  "internal error: user function is neither TF1 or TVectorD");
         }
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get neighbour bins along the specified axis.
///
/// \param[in] bin global bin number
/// \param[in] axis axis number of interest
/// \param[out] prev bin number of previous bin or -1 if not existing
/// \param[out] distPrev distance between bin centres
/// \param[out] next bin number of next bin or -1 if not existing
/// \param[out] distNext distance between bin centres
/// \param[in] isPeriodic (default=false) if true, the first bin is counted as neighbour of the last bin
///
/// return code
///
///   - 0 everything is fine
///   - 1,2,3 isPeriodic option was reset to false, because underflow/overflow
/// bins are present
///
///   - +1 invalid isPeriodic option was specified with underflow bin
///   - +2 invalid isPeriodic option was specified with overflow bin

Int_t TUnfoldBinning::GetBinNeighbours
(Int_t bin,Int_t axis,Int_t *prev,Double_t *distPrev,
 Int_t *next,Double_t *distNext,Bool_t isPeriodic) const
{
   Int_t axisBins[MAXDIM];
   TUnfoldBinning const *distribution=ToAxisBins(bin,axisBins);
   Int_t dimension=distribution->GetDistributionDimension();
   *prev=-1;
   *next=-1;
   *distPrev=0.;
   *distNext=0.;
   Int_t r=0;
   if((axis>=0)&&(axis<dimension)) {
     //TVectorD const *bins=distribution->GetDistributionBinning(axis);
      //Int_t nBin=bins->GetNrows()-1;
      Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
      Int_t centerBin= axisBins[axis];
      axisBins[axis] =centerBin-1;
      if(isPeriodic) {
         if(HasUnderflow(axis)) {
            r +=1;
         } else if((axisBins[axis]<0)&&(nMax>=3)) {
            axisBins[axis]=nMax-1;
         }
      }
      *prev=ToGlobalBin(axisBins);
      if(*prev>=0) {
   *distPrev=distribution->GetDistributionBinCenter(axis,axisBins[axis])-
     distribution->GetDistributionBinCenter(axis,centerBin);
      }
      axisBins[axis] =centerBin+1;
      if(isPeriodic) {
         if(HasOverflow(axis)) {
            r +=2;
         } else if((axisBins[axis]==nMax)&&(nMax>=3)) {
            axisBins[axis]=0;
         }
      }
      *next=ToGlobalBin(axisBins);
      if(*next>=0) {
   *distNext=distribution->GetDistributionBinCenter(axis,axisBins[axis])-
     distribution->GetDistributionBinCenter(axis,centerBin);
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return bit maps indicating underflow and overflow status.
///
/// \param[in] iBin global bin number
/// \param[out] uStatus bit map indicating whether the bin is underflow
/// \param[out] oStatus bit map indicating whether the bin is overflow

void TUnfoldBinning::GetBinUnderflowOverflowStatus
(Int_t iBin,Int_t *uStatus,Int_t *oStatus) const
{
  Int_t axisBins[MAXDIM];
  TUnfoldBinning const *distribution=ToAxisBins(iBin,axisBins);
  Int_t dimension=distribution->GetDistributionDimension();
  *uStatus=0;
  *oStatus=0;
  for(Int_t axis=0;axis<dimension;axis++) {
    TVectorD const *bins=distribution->GetDistributionBinning(axis);
    Int_t nBin=bins->GetNrows()-1;
    if(axisBins[axis]<0) *uStatus |= (1<<axis);
    if(axisBins[axis]>=nBin) *oStatus |= (1<<axis);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether there are bins but no axis.

Bool_t TUnfoldBinning::HasUnconnectedBins(void) const
{
   return (!GetDistributionDimension())&&(GetDistributionNumberOfBins()>0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the bin names of unconnected bins.
///
/// \param[in] bin local bin number

const TObjString *TUnfoldBinning::GetUnconnectedBinName(Int_t bin) const {
   TObjString const *r = nullptr;
   if(HasUnconnectedBins()) {
      if(bin<fAxisLabelList->GetEntriesFast()) {
         r = ((TObjString const *)fAxisLabelList->At(bin));
      }
   }
   return r;
}


////////////////////////////////////////////////////////////////////////////////
/// Get average bin size on the specified axis.
///
/// \param[in] axis axis number
/// \param[in] includeUnderflow whether to include the underflow bin
/// \param[in] includeOverflow whether to include the overflow bin

Double_t TUnfoldBinning::GetDistributionAverageBinSize
(Int_t axis,Bool_t includeUnderflow,Bool_t includeOverflow) const
{
   Double_t r=0.0;
   if((axis>=0)&&(axis<GetDistributionDimension())) {
      TVectorD const *bins=GetDistributionBinning(axis);
      Double_t d=(*bins)[bins->GetNrows()-1]-(*bins)[0];
      Double_t nBins=bins->GetNrows()-1;
      if(includeUnderflow && HasUnderflow(axis)) {
         Double_t w=GetDistributionUnderflowBinWidth(axis);
         if(w>0) {
            nBins++;
            d += w;
         }
      }
      if(includeOverflow && HasOverflow(axis)) {
         Double_t w=GetDistributionOverflowBinWidth(axis);
         if(w>0.0) {
            nBins++;
            d += w;
         }
      }
      if(nBins>0) {
         r=d/nBins;
      }
   } else {
      Error("GetDistributionAverageBinSize","axis %d does not exist",axis);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin width assigned to the underflow bin.
///
/// \param[in] axis axis number
///
/// the bin width of the first bin is returned.
/// The method is virtual, so this behaviour can be adjusted.
///
/// return width of the underflow bin
///   axis: axis number

Double_t TUnfoldBinning::GetDistributionUnderflowBinWidth(Int_t axis) const
{
   TVectorD const *bins=GetDistributionBinning(axis);
   return (*bins)[1]-(*bins)[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin width assigned to the overflow bin.
///
/// \param[in] axis axis number
///
/// the bin width of the last bin is returned.
/// The method is virtual, so this behaviour can be adjusted.
///
/// return width of the underflow bin
///   axis: axis number

Double_t TUnfoldBinning::GetDistributionOverflowBinWidth(Int_t axis) const
{
   TVectorD const *bins=GetDistributionBinning(axis);
   return (*bins)[bins->GetNrows()-1]-(*bins)[bins->GetNrows()-2];
}

////////////////////////////////////////////////////////////////////////////////
/// return bin center for a given axis and bin number
///
/// \param[in] axis axis number
/// \param[in] bin local bin number on the specified axis
///
/// returns the geometrical bin center.
/// for underflow and overflow, the calculation is using the
/// GetDistributionUnderflowBinWidth() and
/// GetDistributionOverflowBinWidth() methods.
///
/// position of the bin center
///  - input:
///     - axis : axis number
///     - bin : bin number on the axis

Double_t TUnfoldBinning::GetDistributionBinCenter
(Int_t axis,Int_t bin) const
{
   TVectorD const *bins=GetDistributionBinning(axis);
   Double_t r=0.0;
   if(bin<0) {
      // underflow bin
      r=(*bins)[0]-0.5*GetDistributionUnderflowBinWidth(axis);
   } else if(bin>=bins->GetNrows()-1) {
      // overflow bin
      r=(*bins)[bins->GetNrows()-1]+0.5*GetDistributionOverflowBinWidth(axis);
   } else {
      r=0.5*((*bins)[bin+1]+(*bins)[bin]);
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Get global bin number, given axis bin numbers.
///
/// \param[in] axisBins[] bin numbers on each axis
/// \param[out] isBelow indicates bins are in underflow but there is
/// no undeflow bin
/// \param[out] isAbove indicates bins are in overflow but there is
/// no overflow bin
///
/// return: global bin number or -1 if not matched.

Int_t TUnfoldBinning::ToGlobalBin
(Int_t const *axisBins,Int_t *isBelow,Int_t *isAbove) const
{
  Int_t dimension=GetDistributionDimension();
  Int_t r=0;
  if(isBelow) *isBelow=0;
  if(isAbove) *isAbove=0;
  if(dimension>0) {
    for(Int_t axis=dimension-1;axis>=0;axis--) {
      Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
      Int_t i=axisBins[axis];
      if(HasUnderflow(axis)) {
   nMax +=1;
   i +=1;
      }
      if(HasOverflow(axis)) nMax +=1;
      if((i>=0)&&(i<nMax)) {
         if(r>=0) r = r*nMax +i;
      } else {
         r=-1;
         if((i<0)&&(isBelow)) *isBelow |= 1<<axis;
         if((i>=nMax)&&(isAbove)) *isAbove |= 1<<axis;
      }
    }
    if(r>=0) {
      r += GetStartBin();
    }
  } else {
    if((axisBins[0]>=0)&&(axisBins[0]<GetDistributionNumberOfBins()))
      r=GetStartBin()+axisBins[0];
    else
       Fatal("ToGlobalBin","bad input %d for dimensionless binning %s %d",
             axisBins[0],(const char *)GetName(),
             GetDistributionNumberOfBins());
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return distribution in which the bin is located
/// and bin numbers on the corresponding axes.
///
/// \param[in] globalBin global bin number
/// \param[out] axisBins local bin numbers of the distribution's axes
///
/// returns the distribution in which the globalBin is located
///  or 0 if the globalBin is outside this node and its children

TUnfoldBinning const *TUnfoldBinning::ToAxisBins
(Int_t globalBin,Int_t *axisBins) const
{
   TUnfoldBinning const *r=0;
   if((globalBin>=GetStartBin())&&(globalBin<GetEndBin())) {
      TUnfoldBinning const *node;
      for(node=GetChildNode();node && !r; node=node->GetNextNode()) {
         r=node->ToAxisBins(globalBin,axisBins);
      }
      if(!r) {
         r=this;
         Int_t i=globalBin-GetStartBin();
         Int_t dimension=GetDistributionDimension();
         if(dimension>0) {
            for(int axis=0;axis<dimension;axis++) {
               Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
               axisBins[axis]=0;
               if(HasUnderflow(axis)) {
                  axisBins[axis] =-1;
                  nMax += 1;
               }
               if(HasOverflow(axis)) nMax +=1;
               axisBins[axis] += i % nMax;
               i /= nMax;
            }
         } else {
            axisBins[0]=i;
         }
      }
   }
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode axis steering.
///
/// \param[in] axisSteering the steering to decode
/// \param[in] options the allowed options to extract
/// \param[out] isOptionGiven array of decoded steering options,
///            the dimension equal to the number of
///            characters in <b>options</b>
///
/// the axis steering is given in the form
///  "axis[option];axis[option];..."
///  axis : the name of the axis for which the optionlist is relevant
/// the character * matches all axes
/// option : a list of characters taken from <b>options</b>
/// for each match the corresponding bit number corresponding to
///  the axis number is set in
///  <b>isOptionGiven</b>[i], where i is the position of the matching option
///  character in <b>options</b>

void TUnfoldBinning::DecodeAxisSteering
(const char *axisSteering,const char *options,Int_t *isOptionGiven) const
{
  Int_t nOpt=TString(options).Length();
  for(Int_t i=0;i<nOpt;i++) isOptionGiven[i]=0;
  if(axisSteering) {
     TObjArray *patterns=TString(axisSteering).Tokenize(";");
     Int_t nPattern=patterns->GetEntries();
     Int_t nAxis=fAxisLabelList->GetEntries();
     for(Int_t i=0;i<nPattern;i++) {
        TString const &pattern=((TObjString const *)patterns->At(i))
           ->GetString();
        Int_t bracketBegin=pattern.Last('[');
        Int_t len=pattern.Length();
        if((bracketBegin>0)&&(pattern[len-1]==']')) {
     TString axisId=pattern(0,bracketBegin);
     Int_t mask=0;
     if((axisId[0]=='*')&&(axisId.Length()==1)) {
       // turn all bins on
       mask=(1<<nAxis)-1;
     } else {
       // if axis is there, turn its bit on
       for(Int_t j=0;j<nAxis;j++) {
         if(!axisId.CompareTo(GetDistributionAxisLabel(j))) {
      mask|= (1<<j);
         }
       }
     }
     for(Int_t o=0;o<nOpt;o++) {
       if(pattern.Last(options[o])>bracketBegin) {
         isOptionGiven[o] |= mask;
       }
     }
        } else {
           Error("DecodeAxisSteering",
                 "steering \"%s\" does not end with [options]",
       (const char *)pattern);
        }
     }
  }
}

