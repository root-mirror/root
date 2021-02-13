// @(#)root/minuit:$Id$
// Author: Anna Kreshuk 04/03/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLinearFitter
#define ROOT_TLinearFitter

//////////////////////////////////////////////////////////////////////////
//
// The Linear Fitter - fitting functions that are LINEAR IN PARAMETERS
//
// Linear fitter is used to fit a set of data points with a linear
// combination of specified functions. Note, that "linear" in the name
// stands only for the model dependency on parameters, the specified
// functions can be nonlinear.
// The general form of this kind of model is
//
//          y(x) = a[0] + a[1]*f[1](x)+...a[n]*f[n](x)
//
// Functions f are fixed functions of x. For example, fitting with a
// polynomial is linear fitting in this sense.
//
//                         The fitting method
//
// The fit is performed using the Normal Equations method with Cholesky
// decomposition.
//
//                         Why should it be used?
//
// The linear fitter is considerably faster than general non-linear
// fitters and doesn't require to set the initial values of parameters.
//
//                          Using the fitter:
//
// 1.Adding the data points:
//  1.1 To store or not to store the input data?
//      - There are 2 options in the constructor - to store or not
//        store the input data. The advantages of storing the data
//        are that you'll be able to reset the fitting model without
//        adding all the points again, and that for very large sets
//        of points the chisquare is calculated more precisely.
//        The obvious disadvantage is the amount of memory used to
//        keep all the points.
//      - Before you start adding the points, you can change the
//        store/not store option by StoreData() method.
//  1.2 The data can be added:
//      - simply point by point - AddPoint() method
//      - an array of points at once:
//        If the data is already stored in some arrays, this data
//        can be assigned to the linear fitter without physically
//        coping bytes, thanks to the Use() method of
//        TVector and TMatrix classes - AssignData() method
//
// 2.Setting the formula
//  2.1 The linear formula syntax:
//      -Additive parts are separated by 2 plus signes "++"
//       --for example "1 ++ x" - for fitting a straight line
//      -All standard functions, undrestood by TFormula, can be used
//       as additive parts
//       --TMath functions can be used too
//      -Functions, used as additive parts, shouldn't have any parameters,
//       even if those parameters are set.
//       --for example, if normalizing a sum of a gaus(0, 1) and a
//         gaus(0, 2), don't use the built-in "gaus" of TFormula,
//         because it has parameters, take TMath::Gaus(x, 0, 1) instead.
//      -Polynomials can be used like "pol3", .."polN"
//      -If fitting a more than 3-dimensional formula, variables should
//       be numbered as follows:
//       -- x0, x1, x2... For example, to fit  "1 ++ x0 ++ x1 ++ x2 ++ x3*x3"
//  2.2 Setting the formula:
//    2.2.1 If fitting a 1-2-3-dimensional formula, one can create a
//          TF123 based on a linear expression and pass this function
//          to the fitter:
//          --Example:
//            TLinearFitter *lf = new TLinearFitter();
//            TF2 *f2 = new TF2("f2", "x ++ y ++ x*x*y*y", -2, 2, -2, 2);
//            lf->SetFormula(f2);
//          --The results of the fit are then stored in the function,
//            just like when the TH1::Fit or TGraph::Fit is used
//          --A linear function of this kind is by no means different
//            from any other function, it can be drawn, evaluated, etc.
//    2.2.2 There is no need to create the function if you don't want to,
//          the formula can be set by expression:
//          --Example:
//            // 2 is the number of dimensions
//            TLinearFitter *lf = new TLinearFitter(2);
//            lf->SetFormula("x ++ y ++ x*x*y*y");
//          --That's the only way to go, if you want to fit in more
//            than 3 dimensions
//    2.2.3 The fastest functions to compute are polynomials and hyperplanes.
//          --Polynomials are set the usual way: "pol1", "pol2",...
//          --Hyperplanes are set by expression "hyp3", "hyp4", ...
//          ---The "hypN" expressions only work when the linear fitter
//             is used directly, not through TH1::Fit or TGraph::Fit.
//             To fit a graph or a histogram with a hyperplane, define
//             the function as "1++x++y".
//          ---A constant term is assumed for a hyperplane, when using
//             the "hypN" expression, so "hyp3" is in fact fitting with
//             "1++x++y++z" function.
//          --Fitting hyperplanes is much faster than fitting other
//            expressions so if performance is vital, calculate the
//            function values beforehand and give them to the fitter
//            as variables
//          --Example:
//            You want to fit "sin(x)|cos(2*x)" very fast. Calculate
//            sin(x) and cos(2*x) beforehand and store them in array *data.
//            Then:
//            TLinearFitter *lf=new TLinearFitter(2, "hyp2");
//            lf->AssignData(npoint, 2, data, y);
//
//  2.3 Resetting the formula
//    2.3.1 If the input data is stored (or added via AssignData() function),
//          the fitting formula can be reset without re-adding all the points.
//          --Example:
//            TLinearFitter *lf=new TLinearFitter("1++x++x*x");
//            lf->AssignData(n, 1, x, y, e);
//            lf->Eval()
//            //looking at the parameter significance, you see,
//            // that maybe the fit will improve, if you take out
//            // the constant term
//            lf->SetFormula("x++x*x");
//            lf->Eval();
//            ...
//    2.3.2 If the input data is not stored, the fitter will have to be
//          cleared and the data will have to be added again to try a
//          different formula.
//
// 3.Accessing the fit results
//  3.1 There are methods in the fitter to access all relevant information:
//      --GetParameters, GetCovarianceMatrix, etc
//      --the t-values of parameters and their significance can be reached by
//        GetParTValue() and GetParSignificance() methods
//  3.2 If fitting with a pre-defined TF123, the fit results are also
//      written into this function.
//
//////////////////////////////////////////////////////////////////////////

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TFormula.h"
#include "TVirtualFitter.h"

#include <map>

class TLinearFitter: public TVirtualFitter {

private:
   TVectorD     fParams;         //vector of parameters
   TMatrixDSym  fParCovar;       //matrix of parameters' covariances
   TVectorD     fTValues;        //T-Values of parameters
   TVectorD     fParSign;        //significance levels of parameters
   TMatrixDSym  fDesign;         //matrix AtA
   TMatrixDSym  fDesignTemp;     //! temporary matrix, used for num.stability
   TMatrixDSym  fDesignTemp2;    //!
   TMatrixDSym  fDesignTemp3;    //!

   TVectorD     fAtb;            //vector Atb
   TVectorD     fAtbTemp;        //! temporary vector, used for num.stability
   TVectorD     fAtbTemp2;       //!
   TVectorD     fAtbTemp3;       //!

   static std::map<TString,TFormula*>   fgFormulaMap;  //! map of basis functions and formula
   TObjArray    fFunctions;      //array of basis functions
   TVectorD     fY;              //the values being fit
   Double_t     fY2;             //sum of square of y, used for chisquare
   Double_t     fY2Temp;         //! temporary variable used for num.stability
   TMatrixD     fX;              //values of x
   TVectorD     fE;              //the errors if they are known
   TFormula     *fInputFunction; //the function being fit
   Double_t     fVal[1000];      //! temporary

   Int_t        fNpoints;        //number of points
   Int_t        fNfunctions;     //number of basis functions
   Int_t        fFormulaSize;    //length of the formula
   Int_t        fNdim;           //number of dimensions in the formula
   Int_t        fNfixed;         //number of fixed parameters
   Int_t        fSpecial;        //=100+n if fitting a polynomial of deg.n
                                 //=200+n if fitting an n-dimensional hyperplane
   char         *fFormula;       //the formula
   Bool_t       fIsSet;          //Has the formula been set?
   Bool_t       fStoreData;      //Is the data stored?
   Double_t     fChisquare;      //Chisquare of the fit

   Int_t        fH;              //number of good points in robust fit
   Bool_t       fRobust;         //true when performing a robust fit
   TBits        fFitsample;      //indices of points, used in the robust fit

   Bool_t       *fFixedParams;   //[fNfixed] array of fixed/released params
   

   void  AddToDesign(Double_t *x, Double_t y, Double_t e);
   void  ComputeTValues();
   Int_t GraphLinearFitter(Double_t h);
   Int_t Graph2DLinearFitter(Double_t h);
   Int_t HistLinearFitter();
   Int_t MultiGraphLinearFitter(Double_t h);

   //robust fitting functions:
   Int_t     Partition(Int_t nmini, Int_t *indsubdat);
   void      RDraw(Int_t *subdat, Int_t *indsubdat);
   void      CreateSubset(Int_t ntotal, Int_t h, Int_t *index);
   Double_t  CStep(Int_t step, Int_t h, Double_t *residuals, Int_t *index, Int_t *subdat, Int_t start, Int_t end);
   Bool_t    Linf();

public:
   TLinearFitter();
   TLinearFitter(Int_t ndim, const char *formula, Option_t *opt="D");
   TLinearFitter(Int_t ndim);
   TLinearFitter(TFormula *function, Option_t *opt="D");
   TLinearFitter(const TLinearFitter& tlf);
   virtual ~TLinearFitter();

   TLinearFitter& operator=(const TLinearFitter& tlf);
   virtual void       Add(TLinearFitter *tlf);
   virtual void       AddPoint(Double_t *x, Double_t y, Double_t e=1);
   virtual void       AddTempMatrices();
   virtual void       AssignData(Int_t npoints, Int_t xncols, Double_t *x, Double_t *y, Double_t *e=0);

   virtual void       Clear(Option_t *option="");
   virtual void       ClearPoints();
   virtual void       Chisquare();
   virtual Int_t      Eval();
   virtual Int_t      EvalRobust(Double_t h=-1);
   virtual Int_t      ExecuteCommand(const char *command, Double_t *args, Int_t nargs);
   virtual void       FixParameter(Int_t ipar);
   virtual void       FixParameter(Int_t ipar, Double_t parvalue);
   virtual void       GetAtbVector(TVectorD &v);
   virtual Double_t   GetChisquare();
   virtual void       GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl=0.95);
   virtual void       GetConfidenceIntervals(TObject *obj, Double_t cl=0.95);
   virtual Double_t*  GetCovarianceMatrix() const;
   virtual void       GetCovarianceMatrix(TMatrixD &matr);
   virtual Double_t   GetCovarianceMatrixElement(Int_t i, Int_t j) const {return fParCovar(i, j);}
   virtual void       GetDesignMatrix(TMatrixD &matr);
   virtual void       GetErrors(TVectorD &vpar);
   virtual Int_t      GetNumberTotalParameters() const {return fNfunctions;}
   virtual Int_t      GetNumberFreeParameters() const {return fNfunctions-fNfixed;}
   virtual Int_t      GetNpoints() { return fNpoints; }
   virtual void       GetParameters(TVectorD &vpar);
   virtual Double_t   GetParameter(Int_t ipar) const {return fParams(ipar);}
   virtual Int_t      GetParameter(Int_t ipar,char* name,Double_t& value,Double_t& /*verr*/,Double_t& /*vlow*/, Double_t& /*vhigh*/) const;
   virtual const char *GetParName(Int_t ipar) const;
   virtual Double_t   GetParError(Int_t ipar) const;
   virtual Double_t   GetParTValue(Int_t ipar);
   virtual Double_t   GetParSignificance(Int_t ipar);
   virtual void       GetFitSample(TBits& bits);
   virtual Double_t   GetY2() const {return fY2;}
   virtual Bool_t     IsFixed(Int_t ipar) const {return fFixedParams[ipar];}
   virtual Int_t      Merge(TCollection *list);
   virtual void       PrintResults(Int_t level, Double_t amin=0) const;
   virtual void       ReleaseParameter(Int_t ipar);
   virtual void       SetBasisFunctions(TObjArray * functions);
   virtual void       SetDim(Int_t n);
   virtual void       SetFormula(const char* formula);
   virtual void       SetFormula(TFormula *function);
   virtual void       StoreData(Bool_t store) {fStoreData=store;}

   virtual Bool_t     UpdateMatrix();

   //dummy functions for TVirtualFitter:
   virtual Double_t  Chisquare(Int_t /*npar*/, Double_t * /*params*/) const {return 0;}
   virtual Int_t     GetErrors(Int_t /*ipar*/,Double_t & /*eplus*/, Double_t & /*eminus*/, Double_t & /*eparab*/, Double_t & /*globcc*/) const {return 0;}

   virtual Int_t     GetStats(Double_t& /*amin*/, Double_t& /*edm*/, Double_t& /*errdef*/, Int_t& /*nvpar*/, Int_t& /*nparx*/) const {return 0;}
   virtual Double_t  GetSumLog(Int_t /*i*/) {return 0;}
   virtual void      SetFitMethod(const char * /*name*/) {;}
   virtual Int_t     SetParameter(Int_t /*ipar*/,const char * /*parname*/,Double_t /*value*/,Double_t /*verr*/,Double_t /*vlow*/, Double_t /*vhigh*/) {return 0;}

   ClassDef(TLinearFitter, 2) //fit a set of data points with a linear combination of functions
};

#endif
