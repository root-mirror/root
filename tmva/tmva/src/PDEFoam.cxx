// @(#)root/tmva $Id$
// Author: S.Jadach, Tancredi Carli, Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoam                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementations                                                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Tancredi Carli   - CERN, Switzerland                                      *
 *      Dominik Dannheim - CERN, Switzerland                                      *
 *      S. Jadach        - Institute of Nuclear Physics, Cracow, Poland           *
 *      Alexander Voigt  - TU Dresden, Germany                                    *
 *      Peter Speckmayer - CERN, Switzerland                                      *
 *                                                                                *
 * Copyright (c) 2008, 2010:                                                      *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::PDEFoam
\ingroup TMVA

Implementation of PDEFoam

The PDEFoam method is an extension of the PDERS method, which uses
self-adapting binning to divide the multi-dimensional phase space
in a finite number of hyper-rectangles (boxes).

For a given number of boxes, the binning algorithm adjusts the size
and position of the boxes inside the multidimensional phase space,
minimizing the variance of the signal and background densities inside
the boxes. The binned density information is stored in binary trees,
allowing for a very fast and memory-efficient classification of
events.

The implementation of the PDEFoam is based on the monte-carlo
integration package TFoam included in the analysis package ROOT.

The class TMVA::PDEFoam defines the default interface for the
PDEFoam variants:

  - PDEFoamEvent
  - PDEFoamDiscriminant
  - PDEFoamTarget
  - PDEFoamMultiTarget
  - PDEFoamDecisionTree

Per default PDEFoam stores in the cells the number of events (event
weights) and therefore acts as an event density estimator.
However, the above listed derived classes override this behaviour
to implement certain PDEFoam variations.

In order to use PDEFoam the user has to set the density estimator
of the type TMVA::PDEFoamDensityBase, which is used to during the foam
build-up.  The default PDEFoam should be used with
PDEFoamEventDensity.
*/

#include "TMVA/PDEFoam.h"

#include "TMVA/Event.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDEFoamKernelBase.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TStyle.h"
#include "TObject.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVectorT.h"
#include "TRandom3.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TObjArray.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

ClassImp(TMVA::PDEFoam);

static const Float_t kHigh= FLT_MAX;
static const Float_t kVlow=-FLT_MAX;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for streamer, user should not use it.

TMVA::PDEFoam::PDEFoam() :
   fName("PDEFoam"),
   fDim(0),
   fNCells(0),
   fNBin(5),
   fNSampl(2000),
   fEvPerBin(0),
   fMaskDiv(0),
   fInhiDiv(0),
   fNoAct(1),
   fLastCe(-1),
   fCells(0),
   fHistEdg(0),
   fRvec(0),
   fPseRan(new TRandom3(4356)),
   fAlpha(0),
   fFoamType(kSeparate),
   fXmin(0),
   fXmax(0),
   fNElements(0),
   fNmin(100),
   fMaxDepth(0),
   fVolFrac(1.0/15.0),
   fFillFoamWithOrigWeights(kFALSE),
   fDTSeparation(kFoam),
   fPeekMax(kTRUE),
   fDistr(NULL),
   fTimer(new Timer(0, "PDEFoam", kTRUE)),
   fVariableNames(new TObjArray()),
   fLogger(new MsgLogger("PDEFoam"))
{
   // fVariableNames may delete it's heap-based content
   if (fVariableNames)
      fVariableNames->SetOwner(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// User constructor, to be employed by the user

TMVA::PDEFoam::PDEFoam(const TString& name) :
   fName(name),
   fDim(0),
   fNCells(1000),
   fNBin(5),
   fNSampl(2000),
   fEvPerBin(0),
   fMaskDiv(0),
   fInhiDiv(0),
   fNoAct(1),
   fLastCe(-1),
   fCells(0),
   fHistEdg(0),
   fRvec(0),
   fPseRan(new TRandom3(4356)),
   fAlpha(0),
   fFoamType(kSeparate),
   fXmin(0),
   fXmax(0),
   fNElements(0),
   fNmin(100),
   fMaxDepth(0),
   fVolFrac(1.0/15.0),
   fFillFoamWithOrigWeights(kFALSE),
   fDTSeparation(kFoam),
   fPeekMax(kTRUE),
   fDistr(NULL),
   fTimer(new Timer(1, "PDEFoam", kTRUE)),
   fVariableNames(new TObjArray()),
   fLogger(new MsgLogger("PDEFoam"))
{
   if(strlen(name) > 128)
      Log() << kFATAL << "Name too long " << name.Data() << Endl;

   // fVariableNames may delete it's heap-based content
   if (fVariableNames)
      fVariableNames->SetOwner(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor

TMVA::PDEFoam::~PDEFoam()
{
   delete fVariableNames;
   delete fTimer;
   if (fDistr)  delete fDistr;
   if (fPseRan) delete fPseRan;
   if (fXmin) { delete [] fXmin;  fXmin=0; }
   if (fXmax) { delete [] fXmax;  fXmax=0; }

   ResetCellElements();
   if(fCells!= 0) {
      for(Int_t i=0; i<fNCells; i++) delete fCells[i]; // PDEFoamCell*[]
      delete [] fCells;
   }
   delete [] fRvec;    //double[]
   delete [] fAlpha;   //double[]
   delete [] fMaskDiv; //int[]
   delete [] fInhiDiv; //int[]

   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor  NOT IMPLEMENTED (NEVER USED)

TMVA::PDEFoam::PDEFoam(const PDEFoam &from) :
   TObject(from)
   , fDim(0)
   , fNCells(0)
   , fNBin(0)
   , fNSampl(0)
   , fEvPerBin(0)
   , fMaskDiv(0)
   , fInhiDiv(0)
   , fNoAct(0)
   , fLastCe(0)
   , fCells(0)
   , fHistEdg(0)
   , fRvec(0)
   , fPseRan(0)
   , fAlpha(0)
   , fFoamType(kSeparate)
   , fXmin(0)
   , fXmax(0)
   , fNElements(0)
   , fNmin(0)
   , fMaxDepth(0)
   , fVolFrac(1.0/15.0)
   , fFillFoamWithOrigWeights(kFALSE)
   , fDTSeparation(kFoam)
   , fPeekMax(kTRUE)
   , fDistr(0)
   , fTimer(0)
   , fVariableNames(0)
   , fLogger(new MsgLogger(*from.fLogger))
{
   Log() << kFATAL << "COPY CONSTRUCTOR NOT IMPLEMENTED" << Endl;

   // fVariableNames may delete it's heap-based content
   if (fVariableNames)
      fVariableNames->SetOwner(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets dimension of cubical space

void TMVA::PDEFoam::SetDim(Int_t kDim)
{
   if (kDim < 1)
      Log() << kFATAL << "<SetDim>: Dimension is zero or negative!" << Endl;

   fDim = kDim;
   if (fXmin) delete [] fXmin;
   if (fXmax) delete [] fXmax;
   fXmin = new Double_t[GetTotDim()];
   fXmax = new Double_t[GetTotDim()];
}

////////////////////////////////////////////////////////////////////////////////
/// set lower foam bound in dimension idim

void TMVA::PDEFoam::SetXmin(Int_t idim, Double_t wmin)
{
   if (idim<0 || idim>=GetTotDim())
      Log() << kFATAL << "<SetXmin>: Dimension out of bounds!" << Endl;

   fXmin[idim]=wmin;
}

////////////////////////////////////////////////////////////////////////////////
/// set upper foam bound in dimension idim

void TMVA::PDEFoam::SetXmax(Int_t idim, Double_t wmax)
{
   if (idim<0 || idim>=GetTotDim())
      Log() << kFATAL << "<SetXmax>: Dimension out of bounds!" << Endl;

   fXmax[idim]=wmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Basic initialization of FOAM invoked by the user.
/// IMPORTANT: Random number generator and the distribution object has to be
/// provided using SetPseRan and SetRho prior to invoking this initializator!
///
/// After the foam is grown, space for 2 variables is reserved in
/// every cell.  They are used for filling the foam cells.

void TMVA::PDEFoam::Create()
{
   Bool_t addStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   if(fPseRan==0) Log() << kFATAL << "Random number generator not set" << Endl;
   if(fDistr==0)  Log() << kFATAL << "Distribution function not set" << Endl;
   if(fDim==0)    Log() << kFATAL << "Zero dimension not allowed" << Endl;

   /////////////////////////////////////////////////////////////////////////
   //                   ALLOCATE SMALL LISTS                              //
   //  it is done globally, not for each cell, to save on allocation time //
   /////////////////////////////////////////////////////////////////////////
   fRvec = new Double_t[fDim];   // Vector of random numbers
   if(fRvec==0)  Log() << kFATAL << "Cannot initialize buffer fRvec" << Endl;

   if(fDim>0){
      fAlpha = new Double_t[fDim];    // sum<1 for internal parametrization of the simplex
      if(fAlpha==0)  Log() << kFATAL << "Cannot initialize buffer fAlpha" << Endl;
   }

   //====== List of directions inhibited for division
   if(fInhiDiv == 0){
      fInhiDiv = new Int_t[fDim];
      for(Int_t i=0; i<fDim; i++) fInhiDiv[i]=0;
   }
   //====== Dynamic mask used in Explore for edge determination
   if(fMaskDiv == 0){
      fMaskDiv = new Int_t[fDim];
      for(Int_t i=0; i<fDim; i++) fMaskDiv[i]=1;
   }
   //====== Initialize list of histograms
   fHistEdg = new TObjArray(fDim);           // Initialize list of histograms
   for(Int_t i=0; i<fDim; i++){
      TString hname, htitle;
      hname   = fName+TString("_HistEdge_");
      hname  += i;
      htitle  = TString("Edge Histogram No. ");
      htitle += i;
      (*fHistEdg)[i] = new TH1D(hname.Data(),htitle.Data(),fNBin,0.0, 1.0); // Initialize histogram for each edge
      ((TH1D*)(*fHistEdg)[i])->Sumw2();
   }

   // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
   //                     BUILD-UP of the FOAM                            //
   // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //

   // prepare PDEFoam for growing
   ResetCellElements(); // reset all cell elements

   // Define and explore root cell(s)
   InitCells();
   Grow();

   TH1::AddDirectory(addStatus);

   // prepare PDEFoam for the filling with events
   ResetCellElements(); // reset all cell elements
} // Create

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// It initializes "root part" of the FOAM of the tree of cells.

void TMVA::PDEFoam::InitCells()
{
   fLastCe =-1;                             // Index of the last cell
   if(fCells!= 0) {
      for(Int_t i=0; i<fNCells; i++) delete fCells[i];
      delete [] fCells;
   }

   fCells = new(nothrow) PDEFoamCell*[fNCells];
   if (!fCells) {
      Log() << kFATAL << "not enough memory to create " << fNCells
            << " cells" << Endl;
   }
   for(Int_t i=0; i<fNCells; i++){
      fCells[i]= new PDEFoamCell(fDim); // Allocate BIG list of cells
      fCells[i]->SetSerial(i);
   }

   /////////////////////////////////////////////////////////////////////////////
   //              Single Root Hypercube                                      //
   /////////////////////////////////////////////////////////////////////////////
   CellFill(1,   0);  //  0-th cell ACTIVE

   // Exploration of the root cell(s)
   for(Long_t iCell=0; iCell<=fLastCe; iCell++){
      Explore( fCells[iCell] );    // Exploration of root cell(s)
   }
}//InitCells

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// It initializes content of the newly allocated active cell.

Int_t TMVA::PDEFoam::CellFill(Int_t status, PDEFoamCell *parent)
{
   PDEFoamCell *cell;
   if (fLastCe==fNCells){
      Log() << kFATAL << "Too many cells" << Endl;
   }
   fLastCe++;   // 0-th cell is the first

   cell = fCells[fLastCe];

   cell->Fill(status, parent, 0, 0);

   cell->SetBest( -1);         // pointer for planning division of the cell
   cell->SetXdiv(0.5);         // factor for division
   Double_t xInt2,xDri2;
   if(parent!=0){
      xInt2  = 0.5*parent->GetIntg();
      xDri2  = 0.5*parent->GetDriv();
      cell->SetIntg(xInt2);
      cell->SetDriv(xDri2);
   }else{
      cell->SetIntg(0.0);
      cell->SetDriv(0.0);
   }
   return fLastCe;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// It explores newly defined cell with help of special short MC sampling.
/// As a result, estimates of kTRUE and drive volume is defined/determined
/// Average and dispersion of the weight distribution will is found along
/// each edge and the best edge (minimum dispersion, best maximum weight)
/// is memorized for future use.
/// The optimal division point for eventual future cell division is
/// determined/recorded. Recorded are also minimum and maximum weight etc.
/// The volume estimate in all (inactive) parent cells is updated.
/// Note that links to parents and initial volume = 1/2 parent has to be
/// already defined prior to calling this routine.
///
/// If fNmin > 0 then the total number of (training) events found in
/// the cell during the exploration is stored in the cell.  This
/// information is used within PeekMax() to avoid splitting cells
/// which contain less than fNmin events.

void TMVA::PDEFoam::Explore(PDEFoamCell *cell)
{
   Double_t wt, dx, xBest=0, yBest;
   Double_t intOld, driOld;

   Long_t iev;
   Double_t nevMC;
   Int_t i, j, k;
   Int_t nProj, kBest;
   Double_t ceSum[5], xproj;

   Double_t event_density = 0;
   Double_t totevents     = 0;
   Double_t toteventsOld  = 0;

   PDEFoamVect  cellSize(fDim);
   PDEFoamVect  cellPosi(fDim);

   cell->GetHcub(cellPosi,cellSize);

   PDEFoamCell  *parent;

   Double_t *xRand = new Double_t[fDim];

   Double_t *volPart=0;

   // calculate volume scale
   Double_t vol_scale = 1.0;
   for (Int_t idim = 0; idim < fDim; ++idim)
      vol_scale *= fXmax[idim] - fXmin[idim];

   cell->CalcVolume();
   dx = cell->GetVolume() * vol_scale;
   intOld = cell->GetIntg(); //memorize old values,
   driOld = cell->GetDriv(); //will be needed for correcting parent cells
   toteventsOld = GetCellElement(cell, 0);

   /////////////////////////////////////////////////////
   //    Special Short MC sampling to probe cell      //
   /////////////////////////////////////////////////////
   ceSum[0]=0;
   ceSum[1]=0;
   ceSum[2]=0;
   ceSum[3]=kHigh;  //wtmin
   ceSum[4]=kVlow;  //wtmax

   for (i=0;i<fDim;i++) ((TH1D *)(*fHistEdg)[i])->Reset(); // Reset histograms

   Double_t nevEff=0.;
   // ||||||||||||||||||||||||||BEGIN MC LOOP|||||||||||||||||||||||||||||
   for (iev=0;iev<fNSampl;iev++){
      MakeAlpha();               // generate uniformly vector inside hypercube

      if (fDim>0) for (j=0; j<fDim; j++) xRand[j]= cellPosi[j] +fAlpha[j]*(cellSize[j]);

      wt         = dx*Eval(xRand, event_density);
      totevents += event_density;

      nProj = 0;
      if (fDim>0) {
         for (k=0; k<fDim; k++) {
            xproj =fAlpha[k];
            ((TH1D *)(*fHistEdg)[nProj])->Fill(xproj,wt);
            nProj++;
         }
      }

      ceSum[0] += wt;    // sum of weights
      ceSum[1] += wt*wt; // sum of weights squared
      ceSum[2]++;        // sum of 1
      if (ceSum[3]>wt) ceSum[3]=wt;  // minimum weight;
      if (ceSum[4]<wt) ceSum[4]=wt;  // maximum weight
      // test MC loop exit condition
      if (ceSum[1]>0) nevEff = ceSum[0]*ceSum[0]/ceSum[1];
      else            nevEff = 0;
      if ( nevEff >= fNBin*fEvPerBin) break;
   }   // ||||||||||||||||||||||||||END MC LOOP|||||||||||||||||||||||||||||
   totevents *= dx;

   if (fNSampl>0) totevents /= fNSampl;

   // make sure that, if root cell is explored, more than zero
   // events were found.
   if (cell==fCells[0] && ceSum[0]<=0.0){
      if (ceSum[0]==0.0)
         Log() << kFATAL << "No events were found during exploration of "
               << "root cell.  Please check PDEFoam parameters nSampl "
               << "and VolFrac." << Endl;
      else
         Log() << kWARNING << "Negative number of events found during "
               << "exploration of root cell" << Endl;
   }

   //------------------------------------------------------------------
   //---  predefine logics of searching for the best division edge ---
   for (k=0; k<fDim;k++){
      fMaskDiv[k] =1;                       // default is all
      if ( fInhiDiv[k]==1) fMaskDiv[k] =0; // inhibit some...
   }
   kBest=-1;

   /////////////////////////////////////////////////////////////////////////////

   nevMC            = ceSum[2];
   Double_t intTrue = ceSum[0]/(nevMC+0.000001);
   Double_t intDriv=0.;

   if (kBest == -1) Varedu(ceSum,kBest,xBest,yBest); // determine the best edge,
   intDriv =sqrt(ceSum[1]/nevMC) -intTrue; // Foam build-up, sqrt(<w**2>) -<w>

   //=================================================================================
   cell->SetBest(kBest);
   cell->SetXdiv(xBest);
   cell->SetIntg(intTrue);
   cell->SetDriv(intDriv);
   SetCellElement(cell, 0, totevents);

   // correct/update integrals in all parent cells to the top of the tree
   Double_t  parIntg, parDriv;
   for (parent = cell->GetPare(); parent!=0; parent = parent->GetPare()){
      parIntg = parent->GetIntg();
      parDriv = parent->GetDriv();
      parent->SetIntg( parIntg   +intTrue -intOld );
      parent->SetDriv( parDriv   +intDriv -driOld );
      SetCellElement( parent, 0, GetCellElement(parent, 0) + totevents - toteventsOld);
   }
   delete [] volPart;
   delete [] xRand;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// In determines the best edge candidate and the position of the cell division plane
/// in case of the variance reduction for future cell division,
/// using results of the MC exploration run stored in fHistEdg

void TMVA::PDEFoam::Varedu(Double_t ceSum[5], Int_t &kBest, Double_t &xBest, Double_t &yBest)
{
   Double_t nent   = ceSum[2];
   // Double_t swAll  = ceSum[0];
   Double_t sswAll = ceSum[1];
   Double_t ssw    = sqrt(sswAll)/sqrt(nent);
   //
   Double_t sswIn,sswOut,xLo,xUp;
   kBest =-1;
   xBest =0.5;
   yBest =1.0;
   Double_t maxGain=0.0;
   // Now go over all projections kProj
   for(Int_t kProj=0; kProj<fDim; kProj++) {
      if( fMaskDiv[kProj]) {
         // initialize search over bins
         // Double_t sigmIn =0.0; Double_t sigmOut =0.0;
         Double_t sswtBest = kHigh;
         Double_t gain =0.0;
         Double_t xMin=0.0; Double_t xMax=0.0;
         // Double loop over all pairs jLo<jUp
         for(Int_t jLo=1; jLo<=fNBin; jLo++) {
            Double_t aswIn=0;  Double_t asswIn=0;
            for(Int_t jUp=jLo; jUp<=fNBin;jUp++) {
               aswIn  +=     ((TH1D *)(*fHistEdg)[kProj])->GetBinContent(jUp);
               asswIn += Sqr(((TH1D *)(*fHistEdg)[kProj])->GetBinError(  jUp));
               xLo=(jLo-1.0)/fNBin;
               xUp=(jUp*1.0)/fNBin;
               // swIn  =        aswIn/nent;
               // swOut = (swAll-aswIn)/nent;
               if ( (xUp-xLo) < std::numeric_limits<double>::epsilon()) sswIn=0.;
               else sswIn = sqrt(asswIn)       /sqrt(nent*(xUp-xLo))     *(xUp-xLo);
               if ( (1.0-xUp+xLo) < std::numeric_limits<double>::epsilon()) sswOut=0.;
               else if ( sswAll-asswIn < std::numeric_limits<double>::epsilon()) sswOut=0.;
               else sswOut= sqrt(sswAll-asswIn)/sqrt(nent*(1.0-xUp+xLo)) *(1.0-xUp+xLo);
               if( (sswIn+sswOut) < sswtBest) {
                  sswtBest = sswIn+sswOut;
                  gain     = ssw-sswtBest;
                  // sigmIn   = sswIn -swIn;  // Debug
                  // sigmOut  = sswOut-swOut; // Debug
                  xMin    = xLo;
                  xMax    = xUp;
               }
            }//jUp
         }//jLo
         Int_t iLo = (Int_t) (fNBin*xMin);
         Int_t iUp = (Int_t) (fNBin*xMax);

         if(gain>=maxGain) {
            maxGain=gain;
            kBest=kProj; // <--- !!!!! The best edge
            xBest=xMin;
            yBest=xMax;
            if(iLo == 0)     xBest=yBest; // The best division point
            if(iUp == fNBin) yBest=xBest; // this is not really used
         }
      }
   } //kProj

   if( (kBest >= fDim) || (kBest<0) )
      Log() << kFATAL << "Something wrong with kBest" << Endl;
}          //PDEFoam::Varedu

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// Provides random vector Alpha  0< Alpha(i) < 1

void TMVA::PDEFoam::MakeAlpha()
{
   // simply generate and load kDim uniform random numbers
   fPseRan->RndmArray(fDim,fRvec);   // kDim random numbers needed
   for(Int_t k=0; k<fDim; k++) fAlpha[k] = fRvec[k];
} //MakeAlpha

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.  It finds cell with maximal
/// driver integral for the purpose of the division.  This function
/// is overridden by the PDEFoam Class to apply cuts on the number
/// of events in the cell (fNmin) and the cell tree depth
/// (GetMaxDepth() > 0) during cell buildup.

Long_t TMVA::PDEFoam::PeekMax()
{
   Long_t iCell = -1;

   Long_t  i;
   Double_t  drivMax, driv, xDiv;
   Bool_t bCutNmin = kTRUE;
   Bool_t bCutMaxDepth = kTRUE;
   //   drivMax = kVlow;
   drivMax = 0;  // only split cells if gain>0 (this also avoids splitting at cell boundary)
   for(i=0; i<=fLastCe; i++) {//without root
      if( fCells[i]->GetStat() == 1 ) {
         // if driver integral < numeric limit, skip cell
         driv = fCells[i]->GetDriv();
         if (driv < std::numeric_limits<float>::epsilon())
            continue;

         // do not split cell at the edges
         xDiv = TMath::Abs(fCells[i]->GetXdiv());
         if (xDiv <= std::numeric_limits<Double_t>::epsilon() ||
             xDiv >= 1.0 - std::numeric_limits<Double_t>::epsilon())
            continue;

         // apply cut on depth
         if (GetMaxDepth() > 0)
            bCutMaxDepth = fCells[i]->GetDepth() < GetMaxDepth();

         // apply Nmin-cut
         if (GetNmin() > 0)
            bCutNmin = GetCellElement(fCells[i], 0) > GetNmin();

         // choose cell
         if(driv > drivMax && bCutNmin && bCutMaxDepth) {
            drivMax = driv;
            iCell = i;
         }
      }
   }

   if (iCell == -1){
      if (!bCutNmin)
         Log() << kVERBOSE << "Warning: No cell with more than "
               << GetNmin() << " events found!" << Endl;
      else if (!bCutMaxDepth)
         Log() << kVERBOSE << "Warning: Maximum depth reached: "
               << GetMaxDepth() << Endl;
      else
         Log() << kWARNING << "<PDEFoam::PeekMax>: no more candidate cells (drivMax>0) found for further splitting." << Endl;
   }

   return(iCell);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// It divides cell iCell into two daughter cells.
/// The iCell is retained and tagged as inactive, daughter cells are appended
/// at the end of the buffer.
/// New vertex is added to list of vertices.
/// List of active cells is updated, iCell removed, two daughters added
/// and their properties set with help of MC sampling (PDEFoam_Explore)
/// Returns Code RC=-1 of buffer limit is reached,  fLastCe=fnBuf.

Int_t TMVA::PDEFoam::Divide(PDEFoamCell *cell)
{
   // Double_t xdiv;
   Int_t   kBest;

   if(fLastCe+1 >= fNCells) Log() << kFATAL << "Buffer limit is reached, fLastCe=fnBuf" << Endl;

   cell->SetStat(0); // reset cell as inactive
   fNoAct++;

   // xdiv  = cell->GetXdiv();
   kBest = cell->GetBest();
   if( kBest<0 || kBest>=fDim ) Log() << kFATAL << "Wrong kBest" << Endl;

   //////////////////////////////////////////////////////////////////
   //           define two daughter cells (active)                 //
   //////////////////////////////////////////////////////////////////

   Int_t d1 = CellFill(1,   cell);
   Int_t d2 = CellFill(1,   cell);
   cell->SetDau0((fCells[d1]));
   cell->SetDau1((fCells[d2]));

   Explore( (fCells[d1]) );
   Explore( (fCells[d2]) );

   return 1;
} // PDEFoam_Divide

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram.
/// Evaluates (training) distribution.

Double_t TMVA::PDEFoam::Eval(Double_t *xRand, Double_t &event_density)
{
   // Transform variable xRand, since Foam boundaries are [0,1] and
   // fDistr is filled with events which range in [fXmin,fXmax]
   //
   // Transformation:  [0, 1] --> [xmin, xmax]
   std::vector<Double_t> xvec;
   xvec.reserve(GetTotDim());
   for (Int_t idim = 0; idim < GetTotDim(); ++idim)
      xvec.push_back( VarTransformInvers(idim, xRand[idim]) );

   return GetDistr()->Density(xvec, event_density);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram used by Create.
/// It grow new cells by the binary division process.
/// This function is overridden by the PDEFoam class to stop the foam buildup process
/// if one of the cut conditions stop the cell split.

void TMVA::PDEFoam::Grow()
{
   fTimer->Init(fNCells);

   Long_t iCell;
   PDEFoamCell* newCell;

   while ( (fLastCe+2) < fNCells ) {  // this condition also checked inside Divide
      iCell = PeekMax(); // peek up cell with maximum driver integral

      if ( (iCell<0) || (iCell>fLastCe) ) {
         Log() << kVERBOSE << "Break: "<< fLastCe+1 << " cells created" << Endl;
         // remove remaining empty cells
         for (Long_t jCell=fLastCe+1; jCell<fNCells; jCell++)
            delete fCells[jCell];
         fNCells = fLastCe+1;
         break;
      }
      newCell = fCells[iCell];

      OutputGrow();

      if ( Divide( newCell )==0) break;  // and divide it into two
   }
   OutputGrow( kTRUE );
   CheckAll(1);   // set arg=1 for more info

   Log() << kVERBOSE << GetNActiveCells() << " active cells created" << Endl;
}// Grow

////////////////////////////////////////////////////////////////////////////////
/// This can be called before Create, after setting kDim
/// It defines which variables are excluded in the process of the cell division.
/// For example 'FoamX->SetInhiDiv(1, 1);' inhibits division of y-variable.

void  TMVA::PDEFoam::SetInhiDiv(Int_t iDim, Int_t inhiDiv)
{
   if(fDim==0) Log() << kFATAL << "SetInhiDiv: fDim=0" << Endl;
   if(fInhiDiv == 0) {
      fInhiDiv = new Int_t[ fDim ];
      for(Int_t i=0; i<fDim; i++) fInhiDiv[i]=0;
   }
   //
   if( ( 0<=iDim) && (iDim<fDim)) {
      fInhiDiv[iDim] = inhiDiv;
   } else
      Log() << kFATAL << "Wrong iDim" << Endl;
}//SetInhiDiv

////////////////////////////////////////////////////////////////////////////////
///  User utility, miscellaneous and debug.
///  Checks all pointers in the tree of cells. This is useful auto-diagnostic.
///  level=0, no printout, failures causes STOP
///  level=1, printout, failures lead to WARNINGS only

void TMVA::PDEFoam::CheckAll(Int_t level)
{
   Int_t errors, warnings;
   PDEFoamCell *cell;
   Long_t iCell;

   errors = 0; warnings = 0;
   if (level==1) Log() << kVERBOSE <<  "Performing consistency checks for created foam" << Endl;
   for(iCell=1; iCell<=fLastCe; iCell++) {
      cell = fCells[iCell];
      //  checking general rules
      if( ((cell->GetDau0()==0) && (cell->GetDau1()!=0) ) ||
          ((cell->GetDau1()==0) && (cell->GetDau0()!=0) ) ) {
         errors++;
         if (level==1) Log() << kFATAL << "ERROR: Cell's no %d has only one daughter " << iCell << Endl;
      }
      if( (cell->GetDau0()==0) && (cell->GetDau1()==0) && (cell->GetStat()==0) ) {
         errors++;
         if (level==1) Log() << kFATAL << "ERROR: Cell's no %d  has no daughter and is inactive " << iCell << Endl;
      }
      if( (cell->GetDau0()!=0) && (cell->GetDau1()!=0) && (cell->GetStat()==1) ) {
         errors++;
         if (level==1) Log() << kFATAL << "ERROR: Cell's no %d has two daughters and is active " << iCell << Endl;
      }

      // checking parents
      if( (cell->GetPare())!=fCells[0] ) { // not child of the root
         if ( (cell != cell->GetPare()->GetDau0()) && (cell != cell->GetPare()->GetDau1()) ) {
            errors++;
            if (level==1) Log() << kFATAL << "ERROR: Cell's no %d parent not pointing to this cell " << iCell << Endl;
         }
      }

      // checking daughters
      if(cell->GetDau0()!=0) {
         if(cell != (cell->GetDau0())->GetPare()) {
            errors++;
            if (level==1)  Log() << kFATAL << "ERROR: Cell's no %d daughter 0 not pointing to this cell " << iCell << Endl;
         }
      }
      if(cell->GetDau1()!=0) {
         if(cell != (cell->GetDau1())->GetPare()) {
            errors++;
            if (level==1) Log() << kFATAL << "ERROR: Cell's no %d daughter 1 not pointing to this cell " << iCell << Endl;
         }
      }
      if(cell->GetVolume()<1E-50) {
         errors++;
         if(level==1) Log() << kFATAL << "ERROR: Cell no. " << iCell << " has Volume of <1E-50" << Endl;
      }
   }// loop after cells;

   // Check for cells with Volume=0
   for(iCell=0; iCell<=fLastCe; iCell++) {
      cell = fCells[iCell];
      if( (cell->GetStat()==1) && (cell->GetVolume()<1E-11) ) {
         errors++;
         if(level==1) Log() << kFATAL << "ERROR: Cell no. " << iCell << " is active but Volume is 0 " <<  Endl;
      }
   }
   // summary
   if(level==1){
      Log() << kVERBOSE << "Check has found " << errors << " errors and " << warnings << " warnings." << Endl;
   }
   if(errors>0){
      Info("CheckAll","Check - found total %d  errors \n",errors);
   }
} // Check

////////////////////////////////////////////////////////////////////////////////
/// Prints geometry of and elements of 'iCell', as well as relations
/// to parent and daughter cells.

void TMVA::PDEFoam::PrintCell(Long_t iCell)
{
   if (iCell < 0 || iCell > fLastCe) {
      Log() << kWARNING << "<PrintCell(iCell=" << iCell
            << ")>: cell number " << iCell << " out of bounds!"
            << Endl;
      return;
   }

   PDEFoamVect cellPosi(fDim), cellSize(fDim);
   fCells[iCell]->GetHcub(cellPosi,cellSize);
   Int_t    kBest = fCells[iCell]->GetBest();
   Double_t xBest = fCells[iCell]->GetXdiv();

   Log() << "Cell[" << iCell << "]={ ";
   Log() << "  " << fCells[iCell] << "  " << Endl;  // extra DEBUG
   Log() << " Xdiv[abs. coord.]="
         << VarTransformInvers(kBest,cellPosi[kBest] + xBest*cellSize[kBest])
         << Endl;
   Log() << " Abs. coord. = (";
   for (Int_t idim=0; idim<fDim; idim++) {
      Log() << "dim[" << idim << "]={"
            << VarTransformInvers(idim,cellPosi[idim]) << ","
            << VarTransformInvers(idim,cellPosi[idim] + cellSize[idim])
            << "}";
      if (idim < fDim-1)
         Log() << ", ";
   }
   Log() << ")" << Endl;
   fCells[iCell]->Print("1");
   // print the cell elements
   Log() << "Elements: [";
   TVectorD *vec = (TVectorD*)fCells[iCell]->GetElement();
   if (vec != NULL){
      for (Int_t i=0; i<vec->GetNrows(); i++){
         if (i>0) Log() << ", ";
         Log() << GetCellElement(fCells[iCell], i);
      }
   } else
      Log() << "not set";
   Log() << "]" << Endl;
   Log()<<"}"<<Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints geometry of ALL cells of the FOAM

void TMVA::PDEFoam::PrintCells(void)
{
   for(Long_t iCell=0; iCell<=fLastCe; iCell++)
      PrintCell(iCell);
}

////////////////////////////////////////////////////////////////////////////////
/// This function fills a weight 'wt' into the PDEFoam cell, which
/// corresponds to the given event 'ev'.  Per default cell element 0
/// is filled with the weight 'wt', and cell element 1 is filled
/// with the squared weight.  This function can be overridden by a
/// subclass in order to change the values stored in the foam cells.

void TMVA::PDEFoam::FillFoamCells(const Event* ev, Float_t wt)
{
   // find corresponding foam cell
   std::vector<Float_t> values  = ev->GetValues();
   std::vector<Float_t> tvalues = VarTransform(values);
   PDEFoamCell *cell = FindCell(tvalues);

   // 0. Element: Sum of weights 'wt'
   // 1. Element: Sum of weights 'wt' squared
   SetCellElement(cell, 0, GetCellElement(cell, 0) + wt);
   SetCellElement(cell, 1, GetCellElement(cell, 1) + wt*wt);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the cell elements from all cells.

void TMVA::PDEFoam::ResetCellElements()
{
   if (!fCells) return;

   Log() << kVERBOSE << "Delete cell elements" << Endl;
   for (Long_t iCell = 0; iCell < fNCells; ++iCell) {
      TObject* elements = fCells[iCell]->GetElement();
      if (elements) {
         delete elements;
         fCells[iCell]->SetElement(NULL);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true, if the value of the given cell is undefined.
/// Default value: kFALSE.  This function can be overridden by
/// sub-classes.

Bool_t TMVA::PDEFoam::CellValueIsUndefined( PDEFoamCell* /* cell */ )
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// This function finds the cell, which corresponds to the given
/// untransformed event vector 'xvec' and return its value, which is
/// given by the parameter 'cv'.  If kernel != NULL, then
/// PDEFoamKernelBase::Estimate() is called on the transformed event
/// variables.
///
/// Parameters:
///
/// - xvec - event vector (untransformed, [fXmin,fXmax])
///
/// - cv - the cell value to return
///
/// - kernel - PDEFoam kernel estimator.  If NULL is given, than the
///   pure cell value is returned
///
/// Return:
///
/// The cell value, corresponding to 'xvec', estimated by the given
/// kernel.

Float_t TMVA::PDEFoam::GetCellValue(const std::vector<Float_t> &xvec, ECellValue cv, PDEFoamKernelBase *kernel)
{
   std::vector<Float_t> txvec(VarTransform(xvec));
   if (kernel == NULL)
      return GetCellValue(FindCell(txvec), cv);
   else
      return kernel->Estimate(this, txvec, cv);
}

////////////////////////////////////////////////////////////////////////////////
/// This function finds all cells, which corresponds to the given
/// (incomplete) untransformed event vector 'xvec' and returns the
/// cell values, according to the parameter 'cv'.
///
/// Parameters:
///
/// - xvec - map for the untransformed vector.  The key (Int_t) is
///   the dimension, and the value (Float_t) is the event
///   coordinate.  Note that not all coordinates have to be
///   specified.
///
/// - cv - cell values to return
///
/// Return:
///
/// cell values from all cells that were found

std::vector<Float_t> TMVA::PDEFoam::GetCellValue( const std::map<Int_t,Float_t>& xvec, ECellValue cv )
{
   // transformed event
   std::map<Int_t,Float_t> txvec;
   for (std::map<Int_t,Float_t>::const_iterator it=xvec.begin(); it!=xvec.end(); ++it)
      txvec.insert(std::pair<Int_t, Float_t>(it->first, VarTransform(it->first, it->second)));

   // find all cells, which correspond to the transformed event
   std::vector<PDEFoamCell*> cells = FindCells(txvec);

   // get the cell values
   std::vector<Float_t> cell_values;
   cell_values.reserve(cells.size());
   for (std::vector<PDEFoamCell*>::const_iterator cell_it=cells.begin();
        cell_it != cells.end(); ++cell_it)
      cell_values.push_back(GetCellValue(*cell_it, cv));

   return cell_values;
}

////////////////////////////////////////////////////////////////////////////////
/// Find cell that contains 'xvec' (in foam coordinates [0,1]).
///
/// Loop to find cell that contains 'xvec' starting at root cell,
/// and traversing binary tree to find the cell quickly.  Note, that
/// if 'xvec' lies outside the foam, the cell which is nearest to
/// 'xvec' is returned.  (The returned pointer should never be
/// NULL.)
///
/// Parameters:
///
/// - xvec - event vector (in foam coordinates [0,1])
///
/// Return:
///
/// PDEFoam cell corresponding to 'xvec'

TMVA::PDEFoamCell* TMVA::PDEFoam::FindCell( const std::vector<Float_t> &xvec ) const
{
   PDEFoamVect  cellPosi0(GetTotDim()), cellSize0(GetTotDim());
   PDEFoamCell *cell, *cell0;

   cell=fCells[0]; // start with root cell
   Int_t idim=0;
   while (cell->GetStat()!=1) { //go down binary tree until cell is found
      idim=cell->GetBest();  // dimension that changed
      cell0=cell->GetDau0();
      cell0->GetHcub(cellPosi0,cellSize0);

      if (xvec.at(idim)<=cellPosi0[idim]+cellSize0[idim])
         cell=cell0;
      else
         cell=(cell->GetDau1());
   }
   return cell;
}

////////////////////////////////////////////////////////////////////////////////
/// This is a helper function for std::vector<PDEFoamCell*>
/// FindCells(...) and a generalisation of PDEFoamCell* FindCell().
/// It saves in 'cells' all cells, which contain the coordinates
/// specifies in 'txvec'.  Note, that not all coordinates have to be
/// specified in 'txvec'.
///
/// Parameters:
///
/// - txvec - event vector in foam coordinates [0,1].  The key is
///   the dimension and the value is the event coordinate.  Note,
///   that not all coordinates have to be specified.
///
/// - cell - cell to start searching with (usually root cell
///   fCells[0])
///
/// - cells - list of cells that were found

void TMVA::PDEFoam::FindCells(const std::map<Int_t, Float_t> &txvec, PDEFoamCell* cell, std::vector<PDEFoamCell*> &cells) const
{
   PDEFoamVect  cellPosi0(GetTotDim()), cellSize0(GetTotDim());
   PDEFoamCell *cell0;
   Int_t idim=0;

   while (cell->GetStat()!=1) { //go down binary tree until cell is found
      idim=cell->GetBest();  // dimension that changed

      // check if dimension 'idim' is specified in 'txvec'
      map<Int_t, Float_t>::const_iterator it = txvec.find(idim);

      if (it != txvec.end()){
         // case 1: cell is splitten in a dimension which is specified
         // in txvec
         cell0=cell->GetDau0();
         cell0->GetHcub(cellPosi0,cellSize0);
         // check, whether left daughter cell contains txvec
         if (it->second <= cellPosi0[idim] + cellSize0[idim])
            cell=cell0;
         else
            cell=cell->GetDau1();
      } else {
         // case 2: cell is splitten in target dimension
         FindCells(txvec, cell->GetDau0(), cells);
         FindCells(txvec, cell->GetDau1(), cells);
         return;
      }
   }
   cells.push_back(cell);
}

////////////////////////////////////////////////////////////////////////////////
/// Find all cells, that contain txvec.  This function can be used,
/// when the dimension of the foam is greater than the dimension of
/// txvec.  E.g. this is the case for multi-target regression.
///
/// Parameters:
///
/// - txvec - event vector of variables, transformed into foam
///   coordinates [0,1].  The size of txvec can be smaller than the
///   dimension of the foam.
///
/// Return value:
///
/// - vector of cells, that fit txvec

std::vector<TMVA::PDEFoamCell*> TMVA::PDEFoam::FindCells(const std::vector<Float_t> &txvec) const
{
   // copy the coordinates from 'txvec' into a map
   std::map<Int_t, Float_t> txvec_map;
   for (UInt_t i=0; i<txvec.size(); ++i)
      txvec_map.insert(std::pair<Int_t, Float_t>(i, txvec.at(i)));

   // the cells found
   std::vector<PDEFoamCell*> cells(0);

   // loop over all target dimensions
   FindCells(txvec_map, fCells[0], cells);

   return cells;
}

////////////////////////////////////////////////////////////////////////////////
/// Find all cells, that contain the coordinates specified in txvec.
/// The key in 'txvec' is the dimension, and the corresponding value
/// is the coordinate.  Note, that not all coordinates have to be
/// specified in txvec.
///
/// Parameters:
///
/// - txvec - map of coordinates (transformed into foam coordinates
///   [0,1])
///
/// Return value:
///
/// - vector of cells, that fit txvec

std::vector<TMVA::PDEFoamCell*> TMVA::PDEFoam::FindCells(const std::map<Int_t, Float_t> &txvec) const
{
   // the cells found
   std::vector<PDEFoamCell*> cells(0);

   // loop over all target dimensions
   FindCells(txvec, fCells[0], cells);

   return cells;
}

////////////////////////////////////////////////////////////////////////////////
/// Draws 1-dimensional foam (= histogram)
///
/// Parameters:
///
/// - cell_value - the cell value to draw
///
/// - nbin - number of bins of result histogram
///
/// - kernel - a PDEFoam kernel.

TH1D* TMVA::PDEFoam::Draw1Dim( ECellValue cell_value, Int_t nbin, PDEFoamKernelBase *kernel )
{
   // avoid plotting of wrong dimensions
   if ( GetTotDim()!=1 )
      Log() << kFATAL << "<Draw1Dim>: function can only be used for 1-dimensional foams!"
            << Endl;

   TString hname("h_1dim");
   TH1D* h1=(TH1D*)gDirectory->Get(hname);
   if (h1) delete h1;
   h1= new TH1D(hname, "1-dimensional Foam", nbin, fXmin[0], fXmax[0]);

   if (!h1) Log() << kFATAL << "ERROR: Can not create histo" << hname << Endl;

   // loop over all bins
   for (Int_t ibinx=1; ibinx<=h1->GetNbinsX(); ++ibinx) {
      // get event vector corresponding to bin
      std::vector<Float_t> txvec;
      txvec.push_back( VarTransform(0, h1->GetBinCenter(ibinx)) );
      Float_t val = 0;
      if (kernel != NULL) {
         // get cell value using the kernel
         val = kernel->Estimate(this, txvec, cell_value);
      } else {
         val = GetCellValue(FindCell(txvec), cell_value);
      }
      // fill value to histogram
      h1->SetBinContent(ibinx, val + h1->GetBinContent(ibinx));
   }

   return h1;
}

////////////////////////////////////////////////////////////////////////////////
/// Project foam variable idim1 and variable idim2 to histogram.
///
/// Parameters:
///
/// - idim1, idim2 - dimensions to project to
///
/// - cell_value - the cell value to draw
///
/// - kernel - a PDEFoam kernel (optional).  If NULL is given, the
///            kernel is ignored and the pure cell values are
///            plotted.
///
/// - nbin - number of bins in x and y direction of result histogram
///          (optional, default is 50).
///
/// Returns:
/// a 2-dimensional histogram

TH2D* TMVA::PDEFoam::Project2( Int_t idim1, Int_t idim2, ECellValue cell_value, PDEFoamKernelBase *kernel, UInt_t nbin )
{
   // avoid plotting of wrong dimensions
   if ((idim1>=GetTotDim()) || (idim1<0) ||
       (idim2>=GetTotDim()) || (idim2<0) ||
       (idim1==idim2) )
      Log() << kFATAL << "<Project2>: wrong dimensions given: "
            << idim1 << ", " << idim2 << Endl;

   // root can not handle too many bins in one histogram --> catch this
   // Furthermore, to have more than 1000 bins in the histogram doesn't make
   // sense.
   if (nbin>1000){
      Log() << kWARNING << "Warning: number of bins too big: " << nbin
            << " Using 1000 bins for each dimension instead." << Endl;
      nbin = 1000;
   } else if (nbin<1) {
      Log() << kWARNING << "Wrong bin number: " << nbin
            << "; set nbin=50" << Endl;
      nbin = 50;
   }

   // create result histogram
   TString hname(Form("h_%d_vs_%d",idim1,idim2));

   // if histogram with this name already exists, delete it
   TH2D* h1=(TH2D*)gDirectory->Get(hname.Data());
   if (h1) delete h1;
   h1= new TH2D(hname.Data(), Form("var%d vs var%d",idim1,idim2), nbin, fXmin[idim1], fXmax[idim1], nbin, fXmin[idim2], fXmax[idim2]);

   if (!h1) Log() << kFATAL << "ERROR: Can not create histo" << hname << Endl;

   // ============== start projection algorithm ================
   // loop over all histogram bins (2-dim)
   for (Int_t xbin = 1; xbin <= h1->GetNbinsX(); ++xbin) {
      for (Int_t ybin = 1; ybin <= h1->GetNbinsY(); ++ybin) {
         // calculate the phase space point, which corresponds to this
         // bin combination
         std::map<Int_t, Float_t> txvec;
         txvec[idim1] = VarTransform(idim1, h1->GetXaxis()->GetBinCenter(xbin));
         txvec[idim2] = VarTransform(idim2, h1->GetYaxis()->GetBinCenter(ybin));

         // find the cells, which corresponds to this phase space
         // point
         std::vector<TMVA::PDEFoamCell*> cells = FindCells(txvec);

         // loop over cells and fill the histogram with the cell
         // values
         Float_t sum_cv = 0; // sum of the cell values
         for (std::vector<TMVA::PDEFoamCell*>::const_iterator it = cells.begin();
              it != cells.end(); ++it) {
            // get cell position and size
            PDEFoamVect cellPosi(GetTotDim()), cellSize(GetTotDim());
            (*it)->GetHcub(cellPosi,cellSize);
            // Create complete event vector from txvec.  The missing
            // coordinates of txvec are set to the cell center.
            std::vector<Float_t> tvec;
            for (Int_t i=0; i<GetTotDim(); ++i) {
               if ( i != idim1 && i != idim2 )
                  tvec.push_back(cellPosi[i] + 0.5*cellSize[i]);
               else
                  tvec.push_back(txvec[i]);
            }
            if (kernel != NULL) {
               // get the cell value using the kernel
               sum_cv += kernel->Estimate(this, tvec, cell_value);
            } else {
               sum_cv += GetCellValue(FindCell(tvec), cell_value);
            }
         }

         // fill the bin content
         h1->SetBinContent(xbin, ybin, sum_cv + h1->GetBinContent(xbin, ybin));
      }
   }

   return h1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the cell value of 'cell' corresponding to the given
/// option 'cv'.  This function should be overridden by the subclass
/// in order to specify which cell elements to return for a given
/// cell value 'cv'.  By default kValue returns cell element 0, and
/// kValueError returns cell element 1.

Float_t TMVA::PDEFoam::GetCellValue(const PDEFoamCell* cell, ECellValue cv)
{
   // calculate cell value (depending on the given option 'cv')
   switch (cv) {

   case kValue:
      return GetCellElement(cell, 0);

   case kValueError:
      return GetCellElement(cell, 1);

   case kValueDensity: {

      Double_t volume  = cell->GetVolume();
      if (volume > numeric_limits<double>::epsilon()) {
         return GetCellValue(cell, kValue)/volume;
      } else {
         if (volume<=0){
            cell->Print("1"); // debug output
            Log() << kWARNING << "<GetCellDensity(cell)>: ERROR: cell volume"
                  << " negative or zero!"
                  << " ==> return cell density 0!"
                  << " cell volume=" << volume
                  << " cell entries=" << GetCellValue(cell, kValue) << Endl;
         } else {
            Log() << kWARNING << "<GetCellDensity(cell)>: WARNING: cell volume"
                  << " close to zero!"
                  << " cell volume: " << volume << Endl;
         }
      }
   }
      return 0;

   case kMeanValue:
      return cell->GetIntg();

   case kRms:
      return cell->GetDriv();

   case kRmsOvMean:
      if (cell->GetIntg() != 0)
         return cell->GetDriv()/cell->GetIntg();
      else
         return 0;

   case kCellVolume:
      return cell->GetVolume();

   default:
      Log() << kFATAL << "<GetCellValue>: unknown cell value" << Endl;
      return 0;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns cell element i of cell 'cell'.  If the cell has no
/// elements or the index 'i' is out of range, than 0 is returned.

Double_t TMVA::PDEFoam::GetCellElement( const PDEFoamCell *cell, UInt_t i ) const
{
   // dynamic_cast doesn't seem to work here ?!
   TVectorD *vec = (TVectorD*)cell->GetElement();

   // if vec is not set or index out of range, return 0
   if (!vec || i >= (UInt_t) vec->GetNrows())
      return 0;

   return (*vec)(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Set cell element i of cell to value.  If the cell element i does
/// not exist, it is created.

void TMVA::PDEFoam::SetCellElement( PDEFoamCell *cell, UInt_t i, Double_t value )
{
   TVectorD *vec = NULL;

   // if no cell elements are set, create TVectorD with i+1 entries,
   // ranging from [0,i]
   if (cell->GetElement() == NULL) {
      vec = new TVectorD(i+1);
      vec->Zero();       // set all values to zero
      (*vec)(i) = value; // set element i to value
      cell->SetElement(vec);
   } else {
      // dynamic_cast doesn't seem to work here ?!
      vec = (TVectorD*)cell->GetElement();
      if (!vec)
         Log() << kFATAL << "<SetCellElement> ERROR: cell element is not a TVectorD*" << Endl;
      // check vector size and resize if necessary
      if (i >= (UInt_t) vec->GetNrows())
         vec->ResizeTo(0,i);
      // set element i to value
      (*vec)(i) = value;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden function of PDEFoam to avoid native foam output.
/// Draw TMVA-process bar instead.

void TMVA::PDEFoam::OutputGrow( Bool_t finished )
{
   if (finished) {
      Log() << kINFO << "Elapsed time: " + fTimer->GetElapsedTime()
            << "                                 " << Endl;
      return;
   }

   Int_t modulo = 1;

   if (fNCells        >= 100) modulo = Int_t(fNCells/100);
   if (fLastCe%modulo == 0)   fTimer->DrawProgressBar( fLastCe );
}

////////////////////////////////////////////////////////////////////////////////
/// Debugging tool which plots the cells of a 2-dimensional PDEFoam
/// as rectangles in C++ format readable for ROOT.
///
/// Parameters:
/// - filename - filename of output root macro
///
/// - opt - cell_value, rms, rms_ov_mean
///   If cell_value is set, the following values will be filled into
///   the result histogram:
///    - number of events - in case of classification with 2 separate
///                         foams or multi-target regression
///    - discriminator    - in case of classification with one
///                         unified foam
///    - target           - in case of mono-target regression
///   If none of {cell_value, rms, rms_ov_mean} is given, the cells
///   will not be filled.
///   If 'opt' contains the string 'cellnumber', the index of
///   each cell is draw in addition.
///
/// - createCanvas - whether to create a new canvas or not
///
/// - colors - whether to fill cells with colors or shades of grey
///
/// Example:
///
///   The following commands load a mono-target regression foam from
///   file 'foam.root' and create a ROOT macro 'output.C', which
///   draws all PDEFoam cells with little boxes.  The latter are
///   filled with colors according to the target value stored in the
///   cell.  Also the cell number is drawn.
///
///   TFile file("foam.root");
///   TMVA::PDEFoam *foam = (TMVA::PDEFoam*) gDirectory->Get("MonoTargetRegressionFoam");
///   foam->RootPlot2dim("output.C","cell_value,cellnumber");
///   gROOT->Macro("output.C");

void TMVA::PDEFoam::RootPlot2dim( const TString& filename, TString opt,
                                  Bool_t createCanvas, Bool_t colors )
{
   if (GetTotDim() != 2)
      Log() << kFATAL << "RootPlot2dim() can only be used with "
            << "two-dimensional foams!" << Endl;

   // select value to plot
   ECellValue cell_value = kValue;
   Bool_t plotcellnumber = kFALSE;
   Bool_t fillcells      = kTRUE;
   if (opt.Contains("cell_value")){
      cell_value = kValue;
   } else if (opt.Contains("rms_ov_mean")){
      cell_value = kRmsOvMean;
   } else if (opt.Contains("rms")){
      cell_value = kRms;
   } else {
      fillcells = kFALSE;
   }
   if (opt.Contains("cellnumber"))
      plotcellnumber = kTRUE;

   // open file (root macro)
   std::ofstream outfile(filename, std::ios::out);

   outfile<<"{" << std::endl;

   // declare boxes and set the fill styles
   if (!colors) { // define grayscale colors from light to dark,
      // starting from color index 1000
      outfile << "TColor *graycolors[100];" << std::endl;
      outfile << "for (Int_t i=0.; i<100; i++)" << std::endl;
      outfile << "  graycolors[i]=new TColor(1000+i, 1-(Float_t)i/100.,1-(Float_t)i/100.,1-(Float_t)i/100.);"<< std::endl;
   }
   if (createCanvas)
      outfile << "cMap = new TCanvas(\"" << fName << "\",\"Cell Map for "
              << fName << "\",600,600);" << std::endl;

   outfile<<"TBox*a=new TBox();"<<std::endl;
   outfile<<"a->SetFillStyle(0);"<<std::endl;  // big frame
   outfile<<"a->SetLineWidth(4);"<<std::endl;
   outfile<<"TBox *b1=new TBox();"<<std::endl;  // single cell
   outfile<<"TText*t=new TText();"<<std::endl;  // text for numbering
   if (fillcells) {
      outfile << (colors ? "gStyle->SetPalette(1, 0);" : "gStyle->SetPalette(0);")
              << std::endl;
      outfile <<"b1->SetFillStyle(1001);"<<std::endl;
      outfile<<"TBox *b2=new TBox();"<<std::endl;  // single cell
      outfile <<"b2->SetFillStyle(0);"<<std::endl;
   }
   else {
      outfile <<"b1->SetFillStyle(0);"<<std::endl;
   }

   if (fillcells)
      (colors ? gStyle->SetPalette(1, 0) : gStyle->SetPalette(0) );

   Float_t zmin = 1E8;  // minimal value (for color calculation)
   Float_t zmax = -1E8; // maximal value (for color calculation)

   // if cells shall be filled, calculate minimal and maximal plot
   // value --> store in zmin and zmax
   if (fillcells) {
      for (Long_t iCell=1; iCell<=fLastCe; iCell++) {
         if ( fCells[iCell]->GetStat() == 1) {
            Float_t value = GetCellValue(fCells[iCell], cell_value);
            if (value<zmin)
               zmin=value;
            if (value>zmax)
               zmax=value;
         }
      }
      outfile << "// observed minimum and maximum of distribution: " << std::endl;
      outfile << "// Float_t zmin = "<< zmin << ";" << std::endl;
      outfile << "// Float_t zmax = "<< zmax << ";" << std::endl;
   }

   outfile << "// used minimum and maximum of distribution (taking into account log scale if applicable): " << std::endl;
   outfile << "Float_t zmin = "<< zmin << ";" << std::endl;
   outfile << "Float_t zmax = "<< zmax << ";" << std::endl;

   Float_t x1,y1,x2,y2,x,y; // box and text coordinates
   Float_t offs  = 0.01;
   Float_t lpag  = 1-2*offs;
   Int_t ncolors  = colors ? gStyle->GetNumberOfColors() : 100;
   Float_t scale = (ncolors-1)/(zmax - zmin);
   PDEFoamVect cellPosi(GetTotDim()), cellSize(GetTotDim());

   // loop over cells and draw a box for every cell (and maybe the
   // cell number as well)
   outfile << "// =========== Rectangular cells  ==========="<< std::endl;
   for (Long_t iCell=1; iCell<=fLastCe; iCell++) {
      if ( fCells[iCell]->GetStat() == 1) {
         fCells[iCell]->GetHcub(cellPosi,cellSize);
         x1 = offs+lpag*(cellPosi[0]);
         y1 = offs+lpag*(cellPosi[1]);
         x2 = offs+lpag*(cellPosi[0]+cellSize[0]);
         y2 = offs+lpag*(cellPosi[1]+cellSize[1]);

         if (fillcells) {
            // get cell value
            Float_t value = GetCellValue(fCells[iCell], cell_value);

            // calculate fill color
            Int_t color;
            if (colors)
               color = gStyle->GetColorPalette(Int_t((value-zmin)*scale));
            else
               color = 1000+(Int_t((value-zmin)*scale));

            // set fill color of box b1
            outfile << "b1->SetFillColor(" << color << ");" << std::endl;
         }

         //     cell rectangle
         outfile<<"b1->DrawBox("<<x1<<","<<y1<<","<<x2<<","<<y2<<");"<<std::endl;
         if (fillcells)
            outfile<<"b2->DrawBox("<<x1<<","<<y1<<","<<x2<<","<<y2<<");"<<std::endl;

         //     cell number
         if (plotcellnumber) {
            outfile<<"t->SetTextColor(4);"<<std::endl;
            if(fLastCe<51)
               outfile<<"t->SetTextSize(0.025);"<<std::endl;  // text for numbering
            else if(fLastCe<251)
               outfile<<"t->SetTextSize(0.015);"<<std::endl;
            else
               outfile<<"t->SetTextSize(0.008);"<<std::endl;
            x = offs+lpag*(cellPosi[0]+0.5*cellSize[0]);
            y = offs+lpag*(cellPosi[1]+0.5*cellSize[1]);
            outfile<<"t->DrawText("<<x<<","<<y<<","<<"\""<<iCell<<"\""<<");"<<std::endl;
         }
      }
   }
   outfile<<"// ============== End Rectangles ==========="<< std::endl;

   outfile << "}" << std::endl;
   outfile.flush();
   outfile.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert event to internal foam's density estimator
/// PDEFoamDensityBase.

void TMVA::PDEFoam::FillBinarySearchTree( const Event* ev )
{
   GetDistr()->FillBinarySearchTree(ev);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the foam's density estimator, which contains the binary
/// search tree.

void TMVA::PDEFoam::DeleteBinarySearchTree()
{
   if(fDistr) delete fDistr;
   fDistr = NULL;
}
