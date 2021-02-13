// @(#)root/hist:$Id$
// Authors: Lorenzo Moneta, Aurélie Flandi  27/08/14

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  ROOT  Team, CERN/PH-SFT                        *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_TF1NormSum__
#define ROOT_TF1NormSum__

#include "TF1AbsComposition.h"
#include <vector>
#include <memory>
#include "TF1.h"

class TF1NormSum : public TF1AbsComposition {

protected:
   unsigned int fNOfFunctions;                         ///< Number of functions to add
   Double_t     fScale;                                ///< Fixed Scale parameter to normalize function (e.g. bin width)
   Double_t fXmin;                                     /// Minimal bound of range of NormSum
   Double_t fXmax;                                     /// Maximal bound of range of NormSum
   std::vector<std::unique_ptr<TF1>> fFunctions;       ///< Vector of size fNOfFunctions containing TF1 functions
   std::vector < Double_t  > fCoeffs;                  ///< Vector of size afNOfFunctions containing coefficients in front of each function
   std::vector < Int_t     > fCstIndexes;              ///< Vector with size of fNOfFunctions containing the index of the constant parameter/ function (the removed ones)
   std::vector< TString >    fParNames;                ///< Parameter names

   void InitializeDataMembers(const std::vector<TF1 *> &functions, const std::vector<Double_t> &coeffs,
                              Double_t scale); // acts as a constructor

public:

   TF1NormSum();
   TF1NormSum(const std::vector <TF1*>&functions, const std::vector <Double_t> &coeffs, Double_t scale = 1.);
   TF1NormSum(TF1* function1, TF1* function2, Double_t coeff1 = 1., Double_t coeff2 = 1., Double_t scale = 1.);
   TF1NormSum(TF1* function1, TF1* function2, TF1*function3, Double_t coeff1 = 1., Double_t coeff2 = 1., Double_t coeff3 = 1., Double_t scale = 1.);
   TF1NormSum(const TString &formula, Double_t xmin, Double_t xmax);

   // Copy constructor
   TF1NormSum(const TF1NormSum &nsum);

   TF1NormSum &operator=(const TF1NormSum &rhs);

   virtual ~TF1NormSum() {}

   double operator()(const Double_t *x, const Double_t *p);

   std::vector<double> GetParameters() const;

   void        SetScale(Double_t scale) { fScale = scale; };

   void SetParameters(const Double_t *params);

   void        SetParameters(Double_t p0, Double_t p1, Double_t p2=0., Double_t p3=0., Double_t p4=0.,
                                   Double_t p5=0., Double_t p6=0., Double_t p7=0., Double_t p8=0., Double_t p9=0., Double_t p10=0.);

   void SetRange(Double_t a, Double_t b);

   Int_t       GetNpar() const;

   Double_t    GetScale() const { return fScale; }

   const char* GetParName(Int_t ipar) const { return fParNames.at(ipar).Data(); }

   Double_t GetXmin() const { return fXmin; }

   Double_t GetXmax() const { return fXmax; }

   void GetRange(Double_t &a, Double_t &b) const;

   void Update();

   void Copy(TObject &obj) const;

   ClassDef(TF1NormSum, 1);
};
#endif /* defined(ROOT_TF1NormSum__) */
