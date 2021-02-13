// Author: Stefan Schmitt
// DESY, 13/10/08

//  Version 17.6,  updated doxygen-style comments, add one argument for scanLCurve
//
//  History:
//    Version 17.5, fix memory leak and other bugs
//    Version 17.4, in parallel to changes in TUnfoldBinning
//    Version 17.3, in parallel to changes in TUnfoldBinning
//    Version 17.2, in parallel to changes in TUnfoldBinning
//    Version 17.1, bug fixes in GetFoldedOutput, GetOutput
//    Version 17.0, error matrix with SetInput, store fL not fLSquared
//    Version 16.2, in parallel to bug-fix in TUnfoldSys
//    Version 16.1, in parallel to bug-fix in TUnfold.C
//    Version 16.0, some cleanup, more getter functions, query version number
//    Version 15, simplified L-curve scan, new tau definition, new eror calc.
//    Version 14, with changes in TUnfoldSys.cxx
//    Version 13, new methods for derived classes
//    Version 12, with support for preconditioned matrix inversion
//    Version 11, regularisation methods have return values
//    Version 10, with bug-fix in TUnfold.cxx
//    Version 9, implements method for optimized inversion of sparse matrix
//    Version 8, replace all TMatrixSparse matrix operations by private code
//    Version 7, fix problem with TMatrixDSparse,TMatrixD multiplication
//    Version 6, completely remove definition of class XY
//    Version 5, move definition of class XY from TUnfold.C to this file
//    Version 4, with bug-fix in TUnfold.C
//    Version 3, with bug-fix in TUnfold.C
//    Version 2, with changed ScanLcurve() arguments
//    Version 1, added ScanLcurve() method
//    Version 0, stable version of basic unfolding algorithm


#ifndef ROOT_TUnfold
#define ROOT_TUnfold

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  TUnfold provides functionality to correct data                      //
//   for migration effects.                                             //
//                                                                      //
//  Citation: S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]        //
//                                                                      //
//                                                                      //
//  TUnfold solves the inverse problem                                  //
//                                                                      //
//                   T   -1            2          T                     //
//    chi**2 = (y-Ax) Vyy  (y-Ax) + tau  (L(x-x0)) L(x-x0)              //
//                                                                      //
//  Monte Carlo input                                                   //
//    y: vector of measured quantities  (dimension ny)                  //
//    Vyy: covariance matrix for y (dimension ny x ny)                  //
//    A: migration matrix               (dimension ny x nx)             //
//    x: unknown underlying distribution (dimension nx)                 //
//  Regularisation                                                      //
//    tau: parameter, defining the regularisation strength              //
//    L: matrix of regularisation conditions (dimension nl x nx)        //
//    x0: underlying distribution bias                                  //
//                                                                      //
//  where chi**2 is minimized as a function of x                        //
//                                                                      //
//  The algorithm is based on "standard" matrix inversion, with the     //
//  known limitations in numerical accuracy and computing cost for      //
//  matrices with large dimensions.                                     //
//                                                                      //
//  Thus the algorithm should not used for large dimensions of x and y  //
//    dim(x) should not exceed O(100)                                   //    
//    dim(y) should not exceed O(500)                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/*
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
*/

#include <TH1D.h>
#include <TH2D.h>
#include <TObject.h>
#include <TArrayI.h>
#include <TSpline.h>
#include <TMatrixDSparse.h>
#include <TMatrixD.h>
#include <TString.h>

#define TUnfold_VERSION "V17.6"
#define TUnfold_CLASS_VERSION 17


class TUnfold : public TObject {
 private:
   void InitTUnfold(void);     // initialize all data members
 public:

   /// type of extra constraint
   enum EConstraint {

      /// use no extra constraint
      kEConstraintNone =0,

      /// enforce preservation of the area
      kEConstraintArea =1
   };

   /// choice of regularisation scheme
   enum ERegMode {

      /// no regularisation, or defined later by RegularizeXXX() methods
      kRegModeNone = 0,

      /// regularise the amplitude of the output distribution
      kRegModeSize = 1,

      /// regularize the 1st derivative of the output distribution
      kRegModeDerivative = 2,

      /// regularize the 2nd derivative of the output distribution
      kRegModeCurvature = 3,


      /// mixed regularisation pattern
      kRegModeMixed = 4
   };

   /// arrangement of axes for the response matrix (TH2 histogram)
   enum EHistMap {

      /// truth level on x-axis of the response matrix
      kHistMapOutputHoriz = 0,

      /// truth level on y-axis of the response matrix
      kHistMapOutputVert = 1
   };

 protected:
   /// response matrix A
   TMatrixDSparse * fA;
   /// regularisation conditions L
   TMatrixDSparse *fL;
   /// input (measured) data y
   TMatrixD *fY;
   /// covariance matrix Vyy corresponding to y
   TMatrixDSparse *fVyy;
   /// scale factor for the bias
   Double_t fBiasScale;
   /// bias vector x0
   TMatrixD *fX0;
   /// regularisation parameter tau squared
   Double_t fTauSquared;
   /// mapping of matrix indices to histogram bins
   TArrayI fXToHist;
   /// mapping of histogram bins to matrix indices
   TArrayI fHistToX;
   /// truth vector calculated from the non-normalized response matrix
   TArrayD fSumOverY;
   /// type of constraint to use for the unfolding
   EConstraint fConstraint;
   /// type of regularisation
   ERegMode fRegMode;
 private:
   /// number of input bins which are dropped because they have error=0
   Int_t fIgnoredBins;
   /// machine accuracy used to determine matrix rank after eigenvalue analysis
   Double_t fEpsMatrix;
   /// unfolding result x
   TMatrixD *fX;
   /// covariance matrix Vxx
   TMatrixDSparse *fVxx;
   /// inverse of covariance matrix Vxx<sup>-1</sup>
   TMatrixDSparse *fVxxInv;
   /// inverse of the input covariance matrix Vyy<sup>-1</sup>
   TMatrixDSparse *fVyyInv;
   /// result x folded back A*x
   TMatrixDSparse *fAx;
   /// chi**2 contribution from (y-Ax)Vyy<sup>-1</sup>(y-Ax)
   Double_t fChi2A;
   /// chi**2 contribution from (x-s*x0)<sup>T</sup>L<sup>T</sup>L(x-s*x0)
   Double_t fLXsquared;
   /// maximum global correlation coefficient
   Double_t fRhoMax;
   /// average global correlation coefficient
   Double_t fRhoAvg;
   /// number of degrees of freedom
   Int_t fNdf;
   /// matrix contribution to the of derivative dx_k/dA_ij
   TMatrixDSparse *fDXDAM[2];
   /// vector contribution to the of derivative dx_k/dA_ij
   TMatrixDSparse *fDXDAZ[2];
   /// derivative of the result wrt tau squared
   TMatrixDSparse *fDXDtauSquared;
   /// derivative of the result wrt dx/dy
   TMatrixDSparse *fDXDY;
   /// matrix E^(-1)
   TMatrixDSparse *fEinv;
   /// matrix E
   TMatrixDSparse *fE;
 protected:
   // Int_t IsNotSymmetric(TMatrixDSparse const &m) const;
   virtual Double_t DoUnfold(void);     // the unfolding algorithm
   virtual void ClearResults(void);     // clear all results
   void ClearHistogram(TH1 *h,Double_t x=0.) const;
   virtual TString GetOutputBinName(Int_t iBinX) const; // name a bin
   TMatrixDSparse *MultiplyMSparseM(const TMatrixDSparse *a,const TMatrixD *b) const; // multiply sparse and non-sparse matrix
   TMatrixDSparse *MultiplyMSparseMSparse(const TMatrixDSparse *a,const TMatrixDSparse *b) const; // multiply sparse and sparse matrix
   TMatrixDSparse *MultiplyMSparseTranspMSparse(const TMatrixDSparse *a,const TMatrixDSparse *b) const; // multiply transposed sparse and sparse matrix
   TMatrixDSparse *MultiplyMSparseMSparseTranspVector
      (const TMatrixDSparse *m1,const TMatrixDSparse *m2,
       const TMatrixTBase<Double_t> *v) const; // calculate M_ij = sum_k [m1_ik*m2_jk*v[k] ]. the pointer v may be zero (means no scaling).
   TMatrixDSparse *InvertMSparseSymmPos(const TMatrixDSparse *A,Int_t *rank) const; // invert symmetric (semi-)positive sparse matrix
   void AddMSparse(TMatrixDSparse *dest,Double_t f,const TMatrixDSparse *src) const; // replacement for dest += f*src
   TMatrixDSparse *CreateSparseMatrix(Int_t nrow,Int_t ncol,Int_t nele,Int_t *row,Int_t *col,Double_t *data) const; // create a TMatrixDSparse from an array
   /// returns internal number of output (truth) matrix rows
   inline Int_t GetNx(void) const {
      return fA->GetNcols();
   }
   /// converts truth histogram bin number to matrix row
   inline Int_t GetRowFromBin(int ix) const { return fHistToX[ix]; }
   /// converts matrix row to truth histogram bin number
   inline Int_t GetBinFromRow(int ix) const { return fXToHist[ix]; }
   /// returns the number of measurement bins
   inline Int_t GetNy(void) const {
      return fA->GetNrows();
   }
   /// vector of the unfolding result
   inline const TMatrixD *GetX(void) const { return fX; }
   /// covariance matrix of the result
   inline const TMatrixDSparse *GetVxx(void) const { return fVxx; }
   /// inverse of covariance matrix of the result
   inline const TMatrixDSparse *GetVxxInv(void) const { return fVxxInv; }
   /// vector of folded-back result
   inline const TMatrixDSparse *GetAx(void) const { return fAx; }
   /// matrix of derivatives dx/dy
   inline const TMatrixDSparse *GetDXDY(void) const { return fDXDY; }
   /// matrix contributions of the derivative dx/dA
   inline const TMatrixDSparse *GetDXDAM(int i) const { return fDXDAM[i]; }
   /// vector contributions of the derivative dx/dA
   inline const TMatrixDSparse *GetDXDAZ(int i) const { return fDXDAZ[i]; }
   /// matrix E<sup>-1</sup>, using internal bin counting
   inline const TMatrixDSparse *GetEinv(void) const { return fEinv; }
   /// matrix E, using internal bin counting
   inline const TMatrixDSparse *GetE(void) const { return fE; }
   /// inverse of covariance matrix of the data y
   inline const TMatrixDSparse *GetVyyInv(void) const { return fVyyInv; }

   void ErrorMatrixToHist(TH2 *ematrix,const TMatrixDSparse *emat,const Int_t *binMap,Bool_t doClear) const; // return an error matrix as histogram
   Double_t GetRhoIFromMatrix(TH1 *rhoi,const TMatrixDSparse *eOrig,const Int_t *binMap,TH2 *invEmat) const; // return global correlation coefficients
   /// vector of derivative dx/dtauSquared, using internal bin counting
   inline const TMatrixDSparse *GetDXDtauSquared(void) const { return fDXDtauSquared; }
   /// delete matrix and invalidate pointer
   static void DeleteMatrix(TMatrixD **m);
   /// delete sparse matrix and invalidate pointer
   static void DeleteMatrix(TMatrixDSparse **m);

   Bool_t AddRegularisationCondition(Int_t i0,Double_t f0,Int_t i1=-1,Double_t f1=0.,Int_t i2=-1,Double_t f2=0.); // add regularisation condition for a triplet of output bins
   Bool_t AddRegularisationCondition(Int_t nEle,const Int_t *indices,const Double_t *rowData); // add a regularisation condition
public:
   static const char*GetTUnfoldVersion(void);
   // Set up response matrix and regularisation scheme
   TUnfold(const TH2 *hist_A, EHistMap histmap,
           ERegMode regmode = kRegModeSize,
           EConstraint constraint=kEConstraintArea);
   // for root streamer and derived classes
   TUnfold(void);
   virtual ~TUnfold(void);
   // define input distribution
   virtual Int_t SetInput(const TH1 *hist_y, Double_t scaleBias=0.0,Double_t oneOverZeroError=0.0,const TH2 *hist_vyy=0,const TH2 *hist_vyy_inv=0);
   // Unfold with given choice of tau and input 
   virtual Double_t DoUnfold(Double_t tau);
   Double_t DoUnfold(Double_t tau,const TH1 *hist_y, Double_t scaleBias=0.0);
   // scan the L curve using successive calls to DoUnfold(Double_t) at various tau
   virtual Int_t ScanLcurve(Int_t nPoint,Double_t tauMin,
                            Double_t tauMax,TGraph **lCurve,
			    TSpline **logTauX=0,TSpline **logTauY=0,
                            TSpline **logTauCurvature=0);

   // access unfolding results
   Double_t GetTau(void) const;
   void GetOutput(TH1 *output,const Int_t *binMap=0) const;
   void GetEmatrix(TH2 *ematrix,const Int_t *binMap=0) const;
   void GetRhoIJ(TH2 *rhoij,const Int_t *binMap=0) const;
   Double_t GetRhoI(TH1 *rhoi,const Int_t *binMap=0,TH2 *invEmat=0) const;
   void GetFoldedOutput(TH1 *folded,const Int_t *binMap=0) const;

   // access input parameters
   void GetProbabilityMatrix(TH2 *A,EHistMap histmap) const;
   void GetNormalisationVector(TH1 *s,const Int_t *binMap=0) const; // get the vector of normalisation factors, equivalent to the initial bias vector
   void GetInput(TH1 *inputData,const Int_t *binMap=0) const; // get input data
   void GetInputInverseEmatrix(TH2 *ematrix);   // get input data inverse of error matrix
   void GetBias(TH1 *bias,const Int_t *binMap=0) const; // get bias (includind biasScale)
   Int_t GetNr(void) const; // number of regularisation conditions
   void GetL(TH2 *l) const; // get matrix of regularisation conditions
   void GetLsquared(TH2 *lsquared) const;

   // access various properties of the result
   /// get maximum global correlation determined in recent unfolding
   inline Double_t GetRhoMax(void) const { return fRhoMax; }
   /// get average global correlation determined in recent unfolding
   inline Double_t GetRhoAvg(void) const { return fRhoAvg; }
   /// get &chi;<sup>2</sup><sub>A</sub> contribution determined in recent unfolding
   inline Double_t GetChi2A(void) const { return fChi2A; }

   Double_t GetChi2L(void) const; // get &chi;<sup>2</sup><sub>L</sub> contribution determined in recent unfolding
   virtual Double_t GetLcurveX(void) const;        // get value on x axis of L curve
   virtual Double_t GetLcurveY(void) const;        // get value on y axis of L curve
   /// get number of degrees of freedom determined in recent unfolding
   ///
   /// This returns the number of valid measurements minus the number
   /// of unfolded truth bins. If the area constraint is active, one
   /// further degree of freedom is subtracted
   inline Int_t GetNdf(void) const { return fNdf; }
   Int_t GetNpar(void) const;  // get number of parameters

   // advanced features
   void SetBias(const TH1 *bias);       // set alternative bias
   void SetConstraint(EConstraint constraint); // set type of constraint for the next unfolding
   Int_t RegularizeSize(int bin, Double_t scale = 1.0);   // regularise the size of one output bin
   Int_t RegularizeDerivative(int left_bin, int right_bin, Double_t scale = 1.0); // regularize difference of two output bins (1st derivative)
   Int_t RegularizeCurvature(int left_bin, int center_bin, int right_bin, Double_t scale_left = 1.0, Double_t scale_right = 1.0);  // regularize curvature of three output bins (2nd derivative)
   Int_t RegularizeBins(int start, int step, int nbin, ERegMode regmode);        // regularize a 1-dimensional curve
   Int_t RegularizeBins2D(int start_bin, int step1, int nbin1, int step2, int nbin2, ERegMode regmode);  // regularize a 2-dimensional grid
   /// get numerical accuracy for Eigenvalue analysis when inverting
   /// matrices with rank problems
   inline Double_t GetEpsMatrix(void) const { return  fEpsMatrix; } 
   /// set numerical accuracy for Eigenvalue analysis when inverting
   /// matrices with rank problems
   void SetEpsMatrix(Double_t eps); // set accuracy for eigenvalue analysis

   ClassDef(TUnfold, TUnfold_CLASS_VERSION) //Unfolding with support for L-curve analysis
};

#endif
