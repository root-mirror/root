// @(#)root/hist:$Id$
// Author: Christian Holm Christensen    1/8/2000

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPrincipal
#define ROOT_TPrincipal

#include "TNamed.h"
#include "TVectorD.h"
#include "TMatrixD.h"

class TList;

class TPrincipal : public TNamed {

protected:
   Int_t       fNumberOfDataPoints;   // Number of data points
   Int_t       fNumberOfVariables;    // Number of variables

   TVectorD    fMeanValues;           // Mean value over all data points
   TVectorD    fSigmas;               // vector of sigmas
   TMatrixD    fCovarianceMatrix;     // Covariance matrix

   TMatrixD    fEigenVectors;         // Eigenvector matrix of trans
   TVectorD    fEigenValues;          // Eigenvalue vector of trans

   TVectorD    fOffDiagonal;          // elements of the tridiagonal

   TVectorD    fUserData;             // Vector of original data points

   Double_t    fTrace;                // Trace of covarience matrix

   TList      *fHistograms;           // List of histograms

   Bool_t      fIsNormalised;         // Normalize matrix?
   Bool_t      fStoreData;            // Should we store input data?

   TPrincipal(const TPrincipal&);
   TPrincipal& operator=(const TPrincipal&);

   void        MakeNormalised();
   void        MakeRealCode(const char *filename, const char *prefix, Option_t *option="");

public:
   TPrincipal();
   virtual ~TPrincipal();
   TPrincipal(Int_t nVariables, Option_t *opt="ND");

   virtual void       AddRow(const Double_t *x);
   virtual void       Browse(TBrowser *b);
   virtual void       Clear(Option_t *option="");
   const TMatrixD    *GetCovarianceMatrix() const {return &fCovarianceMatrix;}
   const TVectorD    *GetEigenValues() const      {return &fEigenValues;}
   const TMatrixD    *GetEigenVectors() const     {return &fEigenVectors;}
   TList             *GetHistograms() const {return fHistograms;}
   const TVectorD    *GetMeanValues() const       {return &fMeanValues;}
   const Double_t    *GetRow(Int_t row);
   const TVectorD    *GetSigmas() const           {return &fSigmas;}
   const TVectorD    *GetUserData() const         {return &fUserData;}
   Bool_t             IsFolder() const { return kTRUE;}
   virtual void       MakeCode(const char *filename ="pca", Option_t *option="");  // *MENU*
   virtual void       MakeHistograms(const char *name = "pca", Option_t *option="epsdx"); // *MENU*
   virtual void       MakeMethods(const char *classname = "PCA", Option_t *option=""); // *MENU*
   virtual void       MakePrincipals();            // *MENU*
   virtual void       P2X(const Double_t *p, Double_t *x, Int_t nTest);
   virtual void       Print(Option_t *opt="MSE") const;         // *MENU*
   virtual void       SumOfSquareResiduals(const Double_t *x, Double_t *s);
   void               Test(Option_t *option="");       // *MENU*
   virtual void       X2P(const Double_t *x, Double_t *p);

   ClassDef(TPrincipal,2) // Principal Components Analysis
}
;

#endif
