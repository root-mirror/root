// Author: Kerstin Tackmann, Andreas Hoecker, Heiko Lacker

/**********************************************************************************
 *                                                                                *
 * Project: TSVDUnfold - data unfolding based on Singular Value Decomposition     *
 * Package: ROOT                                                                  *
 * Class  : TSVDUnfold                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Single class implementation of SVD data unfolding based on:               *
 *          A. Hoecker, V. Kartvelishvili,                                        *
 *          "SVD approach to data unfolding"                                      *
 *          NIM A372, 469 (1996) [hep-ph/9509307]                                 *
 *                                                                                *
 * Authors:                                                                       *
 *      Kerstin Tackmann <Kerstin.Tackmann@cern.ch>   - CERN, Switzerland         *
 *      Andreas Hoecker  <Andreas.Hoecker@cern.ch>    - CERN, Switzerland         *
 *      Heiko Lacker     <lacker@physik.hu-berlin.de> - Humboldt U, Germany       *
 *                                                                                *
 * Copyright (c) 2010:                                                            *
 *      CERN, Switzerland                                                         *
 *      Humboldt University, Germany                                              *
 *                                                                                *
 **********************************************************************************/

////////////////////////////////////////////////////////////////////////////////

/** \class TSVDUnfold
   \ingroup Hist
 SVD Approach to Data Unfolding
 <p>
 Reference: <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a>
 <p>
 TSVDUnfold implements the singular value decomposition based unfolding method (see reference). Currently, the unfolding of one-dimensional histograms is supported, with the same number of bins for the measured and the unfolded spectrum.
 <p>
 The unfolding procedure is based on singular value decomposition of the response matrix. The regularisation of the unfolding is implemented via a discrete minimum-curvature condition.
 <p>
 Monte Carlo inputs:
 <ul>
 <li><tt>xini</tt>: true underlying spectrum (TH1D, n bins)
 <li><tt>bini</tt>: reconstructed spectrum (TH1D, n bins)
 <li><tt>Adet</tt>: response matrix (TH2D, nxn bins)
 </ul>
 Consider the unfolding of a measured spectrum <tt>bdat</tt> with covariance matrix <tt>Bcov</tt> (if not passed explicitly, a diagonal covariance will be built given the errors of <tt>bdat</tt>). The corresponding spectrum in the Monte Carlo is given by <tt>bini</tt>, with the true underlying spectrum given by <tt>xini</tt>. The detector response is described by <tt>Adet</tt>, with <tt>Adet</tt> filled with events (not probabilities) with the true observable on the y-axis and the reconstructed observable on the x-axis.
 <p>
 The measured distribution can be unfolded for any combination of resolution, efficiency and acceptance effects, provided an appropriate definition of <tt>xini</tt> and <tt>Adet</tt>.<br><br>
 <p>
 The unfolding can be performed by
 \code{.cpp}
 TSVDUnfold *tsvdunf = new TSVDUnfold( bdat, Bcov, bini, xini, Adet );
 TH1D* unfresult = tsvdunf->Unfold( kreg );
 \endcode
 where <tt>kreg</tt> determines the regularisation of the unfolding. In general, overregularisation (too small <tt>kreg</tt>) will bias the unfolded spectrum towards the Monte Carlo input, while underregularisation (too large <tt>kreg</tt>) will lead to large fluctuations in the unfolded spectrum. The optimal regularisation can be determined following guidelines in <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a> using the distribution of the <tt>|d_i|</tt> that can be obtained by <tt>tsvdunf->GetD()</tt> and/or using pseudo-experiments.
 <p>
 Covariance matrices on the measured spectrum (for either the total uncertainties or individual sources of uncertainties) can be propagated to covariance matrices using the <tt>GetUnfoldCovMatrix</tt> method, which uses pseudo experiments for the propagation. In addition, <tt>GetAdetCovMatrix</tt> allows for the propagation of the statistical uncertainties on the response matrix using pseudo experiments. The covariance matrix corresponding to <tt>Bcov</tt> is also computed as described in <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a> and can be obtained from <tt>tsvdunf->GetXtau()</tt> and its (regularisation independent) inverse from  <tt>tsvdunf->GetXinv()</tt>. The distribution of singular values can be retrieved using <tt>tsvdunf->GetSV()</tt>.
 <p>
 See also the tutorial for a toy example.
*/

#include <iostream>

#include "TSVDUnfold.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDecompSVD.h"
#include "TRandom3.h"
#include "TMath.h"

ClassImp(TSVDUnfold);

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Alternative constructor
/// User provides data and MC test spectra, as well as detector response matrix, diagonal covariance matrix of measured spectrum built from the uncertainties on measured spectrum

TSVDUnfold::TSVDUnfold( const TH1D *bdat, const TH1D *bini, const TH1D *xini, const TH2D *Adet )
: TObject     (),
fNdim       (0),
fDdim       (2),
fNormalize  (kFALSE),
fKReg       (-1),
fDHist      (NULL),
fSVHist     (NULL),
fXtau       (NULL),
fXinv       (NULL),
fBdat       (bdat),
fBini       (bini),
fXini       (xini),
fAdet       (Adet),
fToyhisto   (NULL),
fToymat     (NULL),
fToyMode    (kFALSE),
fMatToyMode (kFALSE)
{
   if (bdat->GetNbinsX() != bini->GetNbinsX() ||
       bdat->GetNbinsX() != xini->GetNbinsX() ||
       bdat->GetNbinsX() != Adet->GetNbinsX() ||
       bdat->GetNbinsX() != Adet->GetNbinsY()) {
      TString msg = "All histograms must have equal dimension.\n";
      msg += Form( "  Found: dim(bdat)=%i\n",    bdat->GetNbinsX() );
      msg += Form( "  Found: dim(bini)=%i\n",    bini->GetNbinsX() );
      msg += Form( "  Found: dim(xini)=%i\n",    xini->GetNbinsX() );
      msg += Form( "  Found: dim(Adet)=%i,%i\n", Adet->GetNbinsX(), Adet->GetNbinsY() );
      msg += "Please start again!";

      Fatal( "Init", msg, "%s" );
   }

   fBcov = (TH2D*)fAdet->Clone("bcov");

   for(int i=1; i<=fBdat->GetNbinsX(); i++){
      fBcov->SetBinContent(i, i, fBdat->GetBinError(i)*fBdat->GetBinError(i));
      for(int j=1; j<=fBdat->GetNbinsX(); j++){
         if(i==j) continue;
         fBcov->SetBinContent(i,j,0.);
      }
   }
   // Get the input histos
   fNdim = bdat->GetNbinsX();
   fDdim = 2; // This is the derivative used to compute the curvature matrix
}


////////////////////////////////////////////////////////////////////////////////
/// Default constructor
/// Initialisation of TSVDUnfold
/// User provides data and MC test spectra, as well as detector response matrix and the covariance matrix of the measured distribution

TSVDUnfold::TSVDUnfold( const TH1D *bdat, TH2D* Bcov, const TH1D *bini, const TH1D *xini, const TH2D *Adet )
: TObject     (),
fNdim       (0),
fDdim       (2),
fNormalize  (kFALSE),
fKReg       (-1),
fDHist      (NULL),
fSVHist     (NULL),
fXtau       (NULL),
fXinv       (NULL),
fBdat       (bdat),
fBcov       (Bcov),
fBini       (bini),
fXini       (xini),
fAdet       (Adet),
fToyhisto   (NULL),
fToymat     (NULL),
fToyMode    (kFALSE),
fMatToyMode (kFALSE)
{
   if (bdat->GetNbinsX() != bini->GetNbinsX() ||
       bdat->GetNbinsX() != xini->GetNbinsX() ||
       bdat->GetNbinsX() != Bcov->GetNbinsX() ||
       bdat->GetNbinsX() != Bcov->GetNbinsY() ||
       bdat->GetNbinsX() != Adet->GetNbinsX() ||
       bdat->GetNbinsX() != Adet->GetNbinsY()) {
      TString msg = "All histograms must have equal dimension.\n";
      msg += Form( "  Found: dim(bdat)=%i\n",    bdat->GetNbinsX() );
      msg += Form( "  Found: dim(Bcov)=%i,%i\n",    Bcov->GetNbinsX(), Bcov->GetNbinsY() );
      msg += Form( "  Found: dim(bini)=%i\n",    bini->GetNbinsX() );
      msg += Form( "  Found: dim(xini)=%i\n",    xini->GetNbinsX() );
      msg += Form( "  Found: dim(Adet)=%i,%i\n", Adet->GetNbinsX(), Adet->GetNbinsY() );
      msg += "Please start again!";

      Fatal( "Init", msg, "%s" );
   }

   // Get the input histos
   fNdim = bdat->GetNbinsX();
   fDdim = 2; // This is the derivative used to compute the curvature matrix
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TSVDUnfold::TSVDUnfold( const TSVDUnfold& other )
: TObject     ( other ),
fNdim       (other.fNdim),
fDdim       (other.fDdim),
fNormalize  (other.fNormalize),
fKReg       (other.fKReg),
fDHist      (other.fDHist),
fSVHist     (other.fSVHist),
fXtau       (other.fXtau),
fXinv       (other.fXinv),
fBdat       (other.fBdat),
fBcov       (other.fBcov),
fBini       (other.fBini),
fXini       (other.fXini),
fAdet       (other.fAdet),
fToyhisto   (other.fToyhisto),
fToymat     (other.fToymat),
fToyMode    (other.fToyMode),
fMatToyMode (other.fMatToyMode)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TSVDUnfold::~TSVDUnfold()
{
   if(fToyhisto){
      delete fToyhisto;
      fToyhisto = 0;
   }

   if(fToymat){
      delete fToymat;
      fToymat = 0;
   }

   if(fDHist){
      delete fDHist;
      fDHist = 0;
   }

   if(fSVHist){
      delete fSVHist;
      fSVHist = 0;
   }

   if(fXtau){
      delete fXtau;
      fXtau = 0;
   }

   if(fXinv){
      delete fXinv;
      fXinv = 0;
   }

   if(fBcov){
      delete fBcov;
      fBcov = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the unfolding with regularisation parameter kreg

TH1D* TSVDUnfold::Unfold( Int_t kreg )
{
   fKReg = kreg;

   // Make the histos
   if (!fToyMode && !fMatToyMode) InitHistos( );

   // Create vectors and matrices
   TVectorD vb(fNdim), vbini(fNdim), vxini(fNdim), vberr(fNdim);
   TMatrixD mB(fNdim, fNdim), mA(fNdim, fNdim), mCurv(fNdim, fNdim), mC(fNdim, fNdim);

   Double_t eps = 1e-12;
   Double_t sreg;

   // Copy histogams entries into vector
   if (fToyMode) { H2V( fToyhisto, vb ); H2Verr( fToyhisto, vberr ); }
   else          { H2V( fBdat,     vb ); H2Verr( fBdat,     vberr ); }

   H2M( fBcov, mB);
   H2V( fBini, vbini );
   H2V( fXini, vxini );
   if (fMatToyMode) H2M( fToymat, mA );
   else        H2M( fAdet,   mA );

   // Fill and invert the second derivative matrix
   FillCurvatureMatrix( mCurv, mC );

   // Inversion of mC by help of SVD
   TDecompSVD CSVD(mC);
   TMatrixD CUort = CSVD.GetU();
   TMatrixD CVort = CSVD.GetV();
   TVectorD CSV   = CSVD.GetSig();

   TMatrixD CSVM(fNdim, fNdim);
   for (Int_t i=0; i<fNdim; i++) CSVM(i,i) = 1/CSV(i);

   CUort.Transpose( CUort );
   TMatrixD mCinv = (CVort*CSVM)*CUort;

   //Rescale using the data covariance matrix
   TDecompSVD BSVD( mB );
   TMatrixD QT = BSVD.GetU();
   QT.Transpose(QT);
   TVectorD B2SV = BSVD.GetSig();
   TVectorD BSV(B2SV);

   for(int i=0; i<fNdim; i++){
      BSV(i) = TMath::Sqrt(B2SV(i));
   }
   TMatrixD mAtmp(fNdim,fNdim);
   TVectorD vbtmp(fNdim);
   mAtmp *= 0;
   vbtmp *= 0;
   for(int i=0; i<fNdim; i++){
      for(int j=0; j<fNdim; j++){
         if(BSV(i)){
            vbtmp(i) += QT(i,j)*vb(j)/BSV(i);
         }
         for(int m=0; m<fNdim; m++){
            if(BSV(i)){
               mAtmp(i,j) += QT(i,m)*mA(m,j)/BSV(i);
            }
         }
      }
   }
   mA = mAtmp;
   vb = vbtmp;

   // Singular value decomposition and matrix operations
   TDecompSVD ASVD( mA*mCinv );
   TMatrixD Uort = ASVD.GetU();
   TMatrixD Vort = ASVD.GetV();
   TVectorD ASV  = ASVD.GetSig();

   if (!fToyMode && !fMatToyMode) {
      V2H(ASV, *fSVHist);
   }

   TMatrixD Vreg = mCinv*Vort;

   Uort.Transpose(Uort);
   TVectorD vd    = Uort*vb;

   if (!fToyMode && !fMatToyMode) {
      V2H(vd, *fDHist);
   }

   // Damping coefficient
   Int_t k = GetKReg()-1;

   TVectorD vx(fNdim); // Return variable

   // Damping factors
   TVectorD vdz(fNdim);
   TMatrixD Z(fNdim, fNdim);
   for (Int_t i=0; i<fNdim; i++) {
      if (ASV(i)<ASV(0)*eps) sreg = ASV(0)*eps;
      else                   sreg = ASV(i);
      vdz(i) = sreg/(sreg*sreg + ASV(k)*ASV(k));
      Z(i,i) = vdz(i)*vdz(i);
   }
   TVectorD vz = CompProd( vd, vdz );

   TMatrixD VortT(Vort);
   VortT.Transpose(VortT);
   TMatrixD W = mCinv*Vort*Z*VortT*mCinv;

   TMatrixD Xtau(fNdim, fNdim);
   TMatrixD Xinv(fNdim, fNdim);
   Xtau *= 0;
   Xinv *= 0;
   for (Int_t i=0; i<fNdim; i++) {
      for (Int_t j=0; j<fNdim; j++) {
         Xtau(i,j) =  vxini(i) * vxini(j) * W(i,j);

         double a=0;
         for (Int_t m=0; m<fNdim; m++) {
            a += mA(m,i)*mA(m,j);
         }
         if(vxini(i) && vxini(j))
            Xinv(i,j) = a/vxini(i)/vxini(j);
      }
   }

   // Compute the weights
   TVectorD vw = Vreg*vz;

   // Rescale by xini
   vx = CompProd( vw, vxini );

   if(fNormalize){ // Scale result to unit area
      Double_t scale = vx.Sum();
      if (scale > 0){
         vx *= 1.0/scale;
         Xtau *= 1./scale/scale;
         Xinv *= scale*scale;
      }
   }

   if (!fToyMode && !fMatToyMode) {
      M2H(Xtau, *fXtau);
      M2H(Xinv, *fXinv);
   }

   // Get Curvature and also chi2 in case of MC unfolding
   if (!fToyMode && !fMatToyMode) {
      Info( "Unfold", "Unfolding param: %i",k+1 );
      Info( "Unfold", "Curvature of weight distribution: %f", GetCurvature( vw, mCurv ) );
   }

   TH1D* h = (TH1D*)fBdat->Clone("unfoldingresult");
   for(int i=1; i<=fNdim; i++){
      h->SetBinContent(i,0.);
      h->SetBinError(i,0.);
   }
   V2H( vx, *h );

   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine for given input error matrix covariance matrix of unfolded
/// spectrum from toy simulation given the passed covariance matrix on measured spectrum
/// "cov"    - covariance matrix on the measured spectrum, to be propagated
/// "ntoys"  - number of pseudo experiments used for the propagation
/// "seed"   - seed for pseudo experiments
/// Note that this covariance matrix will contain effects of forced normalisation if spectrum is normalised to unit area.

TH2D* TSVDUnfold::GetUnfoldCovMatrix( const TH2D* cov, Int_t ntoys, Int_t seed )
{
   fToyMode = true;
   TH1D* unfres = 0;
   TH2D* unfcov = (TH2D*)fAdet->Clone("unfcovmat");
   unfcov->SetTitle("Toy covariance matrix");
   for(int i=1; i<=fNdim; i++)
      for(int j=1; j<=fNdim; j++)
         unfcov->SetBinContent(i,j,0.);

   // Code for generation of toys (taken from RooResult and modified)
   // Calculate the elements of the upper-triangular matrix L that
   // gives Lt*L = C, where Lt is the transpose of L (the "square-root method")
   TMatrixD L(fNdim,fNdim); L *= 0;

   for (Int_t iPar= 0; iPar < fNdim; iPar++) {

      // Calculate the diagonal term first
      L(iPar,iPar) = cov->GetBinContent(iPar+1,iPar+1);
      for (Int_t k=0; k<iPar; k++) L(iPar,iPar) -= TMath::Power( L(k,iPar), 2 );
      if (L(iPar,iPar) > 0.0) L(iPar,iPar) = TMath::Sqrt(L(iPar,iPar));
      else                    L(iPar,iPar) = 0.0;

      // ...then the off-diagonal terms
      for (Int_t jPar=iPar+1; jPar<fNdim; jPar++) {
         L(iPar,jPar) = cov->GetBinContent(iPar+1,jPar+1);
         for (Int_t k=0; k<iPar; k++) L(iPar,jPar) -= L(k,iPar)*L(k,jPar);
         if (L(iPar,iPar)!=0.) L(iPar,jPar) /= L(iPar,iPar);
         else                  L(iPar,jPar) = 0;
      }
   }

   // Remember it
   TMatrixD *Lt = new TMatrixD(TMatrixD::kTransposed,L);
   TRandom3 random(seed);

   fToyhisto = (TH1D*)fBdat->Clone("toyhisto");
   TH1D *toymean = (TH1D*)fBdat->Clone("toymean");
   for (Int_t j=1; j<=fNdim; j++) toymean->SetBinContent(j,0.);

   // Get the mean of the toys first
   for (int i=1; i<=ntoys; i++) {

      // create a vector of unit Gaussian variables
      TVectorD g(fNdim);
      for (Int_t k= 0; k < fNdim; k++) g(k) = random.Gaus(0.,1.);

      // Multiply this vector by Lt to introduce the appropriate correlations
      g *= (*Lt);

      // Add the mean value offsets and store the results
      for (int j=1; j<=fNdim; j++) {
         fToyhisto->SetBinContent(j,fBdat->GetBinContent(j)+g(j-1));
         fToyhisto->SetBinError(j,fBdat->GetBinError(j));
      }

      unfres = Unfold(GetKReg());

      for (Int_t j=1; j<=fNdim; j++) {
         toymean->SetBinContent(j, toymean->GetBinContent(j) + unfres->GetBinContent(j)/ntoys);
      }
      delete unfres;
      unfres = 0;
   }

   // Reset the random seed
   random.SetSeed(seed);

   //Now the toys for the covariance matrix
   for (int i=1; i<=ntoys; i++) {

      // Create a vector of unit Gaussian variables
      TVectorD g(fNdim);
      for (Int_t k= 0; k < fNdim; k++) g(k) = random.Gaus(0.,1.);

      // Multiply this vector by Lt to introduce the appropriate correlations
      g *= (*Lt);

      // Add the mean value offsets and store the results
      for (int j=1; j<=fNdim; j++) {
         fToyhisto->SetBinContent( j, fBdat->GetBinContent(j)+g(j-1) );
         fToyhisto->SetBinError  ( j, fBdat->GetBinError(j) );
      }
      unfres = Unfold(GetKReg());

      for (Int_t j=1; j<=fNdim; j++) {
         for (Int_t k=1; k<=fNdim; k++) {
            unfcov->SetBinContent(j,k,unfcov->GetBinContent(j,k) + ( (unfres->GetBinContent(j) - toymean->GetBinContent(j))* (unfres->GetBinContent(k) - toymean->GetBinContent(k))/(ntoys-1)) );
         }
      }
      delete unfres;
      unfres = 0;
   }
   delete Lt;
   delete toymean;
   fToyMode = kFALSE;

   return unfcov;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine covariance matrix of unfolded spectrum from finite statistics in
/// response matrix using pseudo experiments
/// "ntoys"  - number of pseudo experiments used for the propagation
/// "seed"   - seed for pseudo experiments

TH2D* TSVDUnfold::GetAdetCovMatrix( Int_t ntoys, Int_t seed )
{
   fMatToyMode = true;
   TH1D* unfres = 0;
   TH2D* unfcov = (TH2D*)fAdet->Clone("unfcovmat");
   unfcov->SetTitle("Toy covariance matrix");
   for(int i=1; i<=fNdim; i++)
      for(int j=1; j<=fNdim; j++)
         unfcov->SetBinContent(i,j,0.);

   //Now the toys for the detector response matrix
   TRandom3 random(seed);

   fToymat = (TH2D*)fAdet->Clone("toymat");
   TH1D *toymean = (TH1D*)fXini->Clone("toymean");
   for (Int_t j=1; j<=fNdim; j++) toymean->SetBinContent(j,0.);

   for (int i=1; i<=ntoys; i++) {
      for (Int_t k=1; k<=fNdim; k++) {
         for (Int_t m=1; m<=fNdim; m++) {
            if (fAdet->GetBinContent(k,m)) {
               fToymat->SetBinContent(k, m, random.Poisson(fAdet->GetBinContent(k,m)));
            }
         }
      }

      unfres = Unfold(GetKReg());

      for (Int_t j=1; j<=fNdim; j++) {
         toymean->SetBinContent(j, toymean->GetBinContent(j) + unfres->GetBinContent(j)/ntoys);
      }
      delete unfres;
      unfres = 0;
   }

   // Reset the random seed
   random.SetSeed(seed);

   for (int i=1; i<=ntoys; i++) {
      for (Int_t k=1; k<=fNdim; k++) {
         for (Int_t m=1; m<=fNdim; m++) {
            if (fAdet->GetBinContent(k,m))
               fToymat->SetBinContent(k, m, random.Poisson(fAdet->GetBinContent(k,m)));
         }
      }

      unfres = Unfold(GetKReg());

      for (Int_t j=1; j<=fNdim; j++) {
         for (Int_t k=1; k<=fNdim; k++) {
            unfcov->SetBinContent(j,k,unfcov->GetBinContent(j,k) + ( (unfres->GetBinContent(j) - toymean->GetBinContent(j))*(unfres->GetBinContent(k) - toymean->GetBinContent(k))/(ntoys-1)) );
         }
      }
      delete unfres;
      unfres = 0;
   }
   delete toymean;
   fMatToyMode = kFALSE;

   return unfcov;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns d vector (for choosing appropriate regularisation)

TH1D* TSVDUnfold::GetD() const
{
   for (int i=1; i<=fDHist->GetNbinsX(); i++) {
      if (fDHist->GetBinContent(i)<0.) fDHist->SetBinContent(i, TMath::Abs(fDHist->GetBinContent(i)));
   }
   return fDHist;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns singular values vector

TH1D* TSVDUnfold::GetSV() const
{
   return fSVHist;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the computed regularized covariance matrix corresponding to total uncertainties on measured spectrum as passed in the constructor.
/// Note that this covariance matrix will not contain the effects of forced normalization if spectrum is normalized to unit area.

TH2D* TSVDUnfold::GetXtau() const
{
   return fXtau;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the computed inverse of the covariance matrix

TH2D* TSVDUnfold::GetXinv() const
{
   return fXinv;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the covariance matrix

TH2D* TSVDUnfold::GetBCov() const
{
   return fBcov;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill 1D histogram into vector

void TSVDUnfold::H2V( const TH1D* histo, TVectorD& vec )
{
   for (Int_t i=0; i<histo->GetNbinsX(); i++) vec(i) = histo->GetBinContent(i+1);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill 1D histogram errors into vector

void TSVDUnfold::H2Verr( const TH1D* histo, TVectorD& vec )
{
   for (Int_t i=0; i<histo->GetNbinsX(); i++) vec(i) = histo->GetBinError(i+1);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector into 1D histogram

void TSVDUnfold::V2H( const TVectorD& vec, TH1D& histo )
{
   for(Int_t i=0; i<vec.GetNrows(); i++) histo.SetBinContent(i+1, vec(i));
}

////////////////////////////////////////////////////////////////////////////////
/// Fill 2D histogram into matrix

void TSVDUnfold::H2M( const TH2D* histo, TMatrixD& mat )
{
   for (Int_t j=0; j<histo->GetNbinsX(); j++) {
      for (Int_t i=0; i<histo->GetNbinsY(); i++) {
         mat(i,j) = histo->GetBinContent(i+1,j+1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill 2D histogram into matrix

void TSVDUnfold::M2H( const TMatrixD& mat, TH2D& histo )
{
   for (Int_t j=0; j<mat.GetNcols(); j++) {
      for (Int_t i=0; i<mat.GetNrows(); i++) {
         histo.SetBinContent(i+1,j+1, mat(i,j));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Divide entries of two vectors

TVectorD TSVDUnfold::VecDiv( const TVectorD& vec1, const TVectorD& vec2, Int_t zero )
{
   TVectorD quot(vec1.GetNrows());
   for (Int_t i=0; i<vec1.GetNrows(); i++) {
      if (vec2(i) != 0) quot(i) = vec1(i) / vec2(i);
      else {
         if   (zero) quot(i) = 0;
         else        quot(i) = vec1(i);
      }
   }
   return quot;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide matrix entries by vector

TMatrixD TSVDUnfold::MatDivVec( const TMatrixD& mat, const TVectorD& vec, Int_t zero )
{
   TMatrixD quotmat(mat.GetNrows(), mat.GetNcols());
   for (Int_t i=0; i<mat.GetNrows(); i++) {
      for (Int_t j=0; j<mat.GetNcols(); j++) {
         if (vec(i) != 0) quotmat(i,j) = mat(i,j) / vec(i);
         else {
            if   (zero) quotmat(i,j) = 0;
            else        quotmat(i,j) = mat(i,j);
         }
      }
   }
   return quotmat;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply entries of two vectors

TVectorD TSVDUnfold::CompProd( const TVectorD& vec1, const TVectorD& vec2 )
{
   TVectorD res(vec1.GetNrows());
   for (Int_t i=0; i<vec1.GetNrows(); i++) res(i) = vec1(i) * vec2(i);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute curvature of vector

Double_t TSVDUnfold::GetCurvature(const TVectorD& vec, const TMatrixD& curv)
{
   return vec*(curv*vec);
}

////////////////////////////////////////////////////////////////////////////////

void TSVDUnfold::FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const
{
   Double_t eps = 0.00001;

   Int_t ndim = tCurv.GetNrows();

   // Init
   tCurv *= 0;
   tC    *= 0;

   if (fDdim == 0) for (Int_t i=0; i<ndim; i++) tC(i,i) = 1;
   else if (fDdim == 1) {
      for (Int_t i=0; i<ndim; i++) {
         if (i < ndim-1) tC(i,i+1) = 1.0;
         tC(i,i) = 1.0;
      }
   }
   else if (fDdim == 2) {
      for (Int_t i=0; i<ndim; i++) {
         if (i > 0)      tC(i,i-1) = 1.0;
         if (i < ndim-1) tC(i,i+1) = 1.0;
         tC(i,i) = -2.0;
      }
      tC(0,0) = -1.0;
      tC(ndim-1,ndim-1) = -1.0;
   }
   else if (fDdim == 3) {
      for (Int_t i=1; i<ndim-2; i++) {
         tC(i,i-1) =  1.0;
         tC(i,i)   = -3.0;
         tC(i,i+1) =  3.0;
         tC(i,i+2) = -1.0;
      }
   }
   else if (fDdim==4) {
      for (Int_t i=0; i<ndim; i++) {
         if (i > 0)      tC(i,i-1) = -4.0;
         if (i < ndim-1) tC(i,i+1) = -4.0;
         if (i > 1)      tC(i,i-2) =  1.0;
         if (i < ndim-2) tC(i,i+2) =  1.0;
         tC(i,i) = 6.0;
      }
      tC(0,0) = 2.0;
      tC(ndim-1,ndim-1) = 2.0;
      tC(0,1) = -3.0;
      tC(ndim-2,ndim-1) = -3.0;
      tC(1,0) = -3.0;
      tC(ndim-1,ndim-2) = -3.0;
      tC(1,1) =  6.0;
      tC(ndim-2,ndim-2) =  6.0;
   }
   else if (fDdim == 5) {
      for (Int_t i=2; i < ndim-3; i++) {
         tC(i,i-2) = 1.0;
         tC(i,i-1) = -5.0;
         tC(i,i)   = 10.0;
         tC(i,i+1) = -10.0;
         tC(i,i+2) = 5.0;
         tC(i,i+3) = -1.0;
      }
   }
   else if (fDdim == 6) {
      for (Int_t i = 3; i < ndim - 3; i++) {
         tC(i,i-3) = 1.0;
         tC(i,i-2) = -6.0;
         tC(i,i-1) = 15.0;
         tC(i,i)   = -20.0;
         tC(i,i+1) = 15.0;
         tC(i,i+2) = -6.0;
         tC(i,i+3) = 1.0;
      }
   }

   // Add epsilon to avoid singularities
   for (Int_t i=0; i<ndim; i++) tC(i,i) = tC(i,i) + eps;

   //Get curvature matrix
   for (Int_t i=0; i<ndim; i++) {
      for (Int_t j=0; j<ndim; j++) {
         for (Int_t k=0; k<ndim; k++) {
            tCurv(i,j) = tCurv(i,j) + tC(k,i)*tC(k,j);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TSVDUnfold::InitHistos( )
{
   fDHist = new TH1D( "dd", "d vector after orthogonal transformation", fNdim, 0, fNdim );
   fDHist->Sumw2();

   fSVHist = new TH1D( "sv", "Singular values of AC^-1", fNdim, 0, fNdim );
   fSVHist->Sumw2();

   fXtau = (TH2D*)fAdet->Clone("Xtau");
   fXtau->SetTitle("Regularized covariance matrix");
   fXtau->Sumw2();

   fXinv = (TH2D*)fAdet->Clone("Xinv");
   fXinv->SetTitle("Inverse covariance matrix");
   fXinv->Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// naive regularised inversion cuts off small elements

void TSVDUnfold::RegularisedSymMatInvert( TMatrixDSym& mat, Double_t eps )
{
   // init reduced matrix
   const UInt_t n = mat.GetNrows();
   UInt_t nn = 0;

   UInt_t *ipos = new UInt_t[n];
   //   UInt_t ipos[n];

   // find max diagonal entries
   Double_t ymax = 0;
   for (UInt_t i=0; i<n; i++) if (TMath::Abs(mat[i][i]) > ymax) ymax = TMath::Abs(mat[i][i]);

   for (UInt_t i=0; i<n; i++) {

      // save position of accepted entries
      if (TMath::Abs(mat[i][i])/ymax > eps) ipos[nn++] = i;
   }

   // effective matrix
   TMatrixDSym matwork( nn );
   for (UInt_t in=0; in<nn; in++) for (UInt_t jn=0; jn<nn; jn++) matwork(in,jn) = 0;

   // fill non-zero effective working matrix
   for (UInt_t in=0; in<nn; in++) {

      matwork[in][in] = mat[ipos[in]][ipos[in]];
      for (UInt_t jn=in+1; jn<nn; jn++) {
         matwork[in][jn] = mat[ipos[in]][ipos[jn]];
         matwork[jn][in] = matwork[in][jn];
      }
   }

   // invert
   matwork.Invert();

   // reinitialise old matrix
   for (UInt_t i=0; i<n; i++) for (UInt_t j=0; j<n; j++) mat[i][j] = 0;

   // refill inverted matrix in old one
   for (UInt_t in=0; in<nn; in++) {
      mat[ipos[in]][ipos[in]] = matwork[in][in];
      for (UInt_t jn=in+1; jn<nn; jn++) {
         mat[ipos[in]][ipos[jn]] = matwork[in][jn];
         mat[ipos[jn]][ipos[in]] = mat[ipos[in]][ipos[jn]];
      }
   }
   delete []  ipos;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper routine to compute chi-squared between distributions using the computed inverse of the covariance matrix for the unfolded spectrum as given in paper.

Double_t TSVDUnfold::ComputeChiSquared( const TH1D& truspec, const TH1D& unfspec)
{
   UInt_t n = truspec.GetNbinsX();

   // compute chi2
   Double_t chi2 = 0;
   for (UInt_t i=0; i<n; i++) {
      for (UInt_t j=0; j<n; j++) {
         chi2 += ( (truspec.GetBinContent( i+1 )-unfspec.GetBinContent( i+1 )) *
                  (truspec.GetBinContent( j+1 )-unfspec.GetBinContent( j+1 )) * fXinv->GetBinContent(i+1,j+1) );
      }
   }

   return chi2;
}

