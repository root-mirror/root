// @(#)root/hist:$Id$
// Authors: Lorenzo Moneta, Aurélie Flandi  27/08/14

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  ROOT  Team, CERN/PH-SFT                        *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_TF1Convolution__
#define ROOT_TF1Convolution__

#include "TF1AbsComposition.h"
#include <memory>
#include <vector>
#include "TF1.h"
#include "TGraph.h"

class TF1Convolution : public TF1AbsComposition {
   std::unique_ptr<TF1> fFunction1;    ///< First function to be convolved
   std::unique_ptr<TF1> fFunction2;    ///< Second function to be convolved
   std::unique_ptr<TGraph> fGraphConv; ///<! Graph of the convolution

   std::vector < Double_t >   fParams1;
   std::vector < Double_t >   fParams2;

   std::vector< TString >    fParNames;  ///< Parameters' names

   Double_t fXmin;                       ///< Minimal bound of the range of the convolution
   Double_t fXmax;                       ///< Maximal bound of the range of the convolution
   Int_t    fNofParams1;
   Int_t    fNofParams2;
   Int_t    fCstIndex;                   ///< Index of the constant parameter f the first function
   Int_t    fNofPoints;                  ///< Number of point for FFT array
   Bool_t   fFlagFFT;                    ///< Choose FFT or numerical convolution
   Bool_t fFlagGraph = false;            ///<! Tells if the graph is already done or not

   Double_t EvalNumConv(Double_t t);
   Double_t EvalFFTConv(Double_t t);
   void     InitializeDataMembers(TF1* function1, TF1* function2, Bool_t useFFT);
   void     MakeFFTConv();

public:
   TF1Convolution();
   TF1Convolution(TF1 *function1, TF1 *function2, Bool_t useFFT = true);
   TF1Convolution(TF1 *function1, TF1 *function2, Double_t xmin, Double_t xmax, Bool_t useFFT = true);
   TF1Convolution(TString formula, Double_t xmin = 1., Double_t xmax = 0., Bool_t useFFT = true);
   TF1Convolution(TString formula1, TString formula2, Double_t xmin = 1., Double_t xmax = 0., Bool_t useFFT = true);

   // Copy constructor
   TF1Convolution(const TF1Convolution &conv);

   TF1Convolution &operator=(const TF1Convolution &rhs);
   virtual ~TF1Convolution() {}

   void SetParameters(const Double_t *params);
   void SetParameters(Double_t p0, Double_t p1, Double_t p2 = 0., Double_t p3 = 0., Double_t p4 = 0., Double_t p5 = 0.,
                      Double_t p6 = 0., Double_t p7 = 0.);
   void SetRange(Double_t a, Double_t b);
   void SetExtraRange(Double_t percentage);
   void SetNofPointsFFT(Int_t n);
   void SetNumConv(Bool_t flag = true) { fFlagFFT = !flag; }

   Int_t GetNpar() const { return (fNofParams1 + fNofParams2); }
   Double_t GetXmin() const { return fXmin; }
   Double_t GetXmax() const { return fXmax; }
   const char *GetParName(Int_t ipar) const { return fParNames.at(ipar).Data(); }
   void GetRange(Double_t &a, Double_t &b) const;

   void Update();

   Double_t operator()(const Double_t *x, const Double_t *p);

   void Copy(TObject &obj) const;

   ClassDef(TF1Convolution, 1);
};


#endif
