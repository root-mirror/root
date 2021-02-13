// @(#)root/minuit:$Id$
// Author: Anna Kreshuk 04/03/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TLinearFitter.h"
#include "TMath.h"
#include "TDecompChol.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TRandom2.h"
#include "TObjString.h"
#include "TF2.h"
#include "TH1.h"
#include "TList.h"
#include "TClass.h"
#include "TBuffer.h"
#include "TROOT.h"
#include "strlcpy.h"
#include "snprintf.h"

ClassImp(TLinearFitter);


std::map<TString,TFormula*> TLinearFitter::fgFormulaMap;



/**

\class TLinearFitter

\ingroup MinuitOld

The Linear Fitter - For fitting functions that are LINEAR IN PARAMETERS

## The Linear Fitter

Linear fitter is used to fit a set of data points with a linear
combination of specified functions. Note, that "linear" in the name
stands only for the model dependency on parameters, the specified
functions can be nonlinear.
The general form of this kind of model is
~~~~
         y(x) = a[0] + a[1]*f[1](x)+...a[n]*f[n](x)
~~~~

Functions f are fixed functions of x. For example, fitting with a
polynomial is linear fitting in this sense.

### Introduction

####   The fitting method

The fit is performed using the Normal Equations method with Cholesky
decomposition.

####                        Why should it be used?

The linear fitter is considerably faster than general non-linear
fitters and doesn't require to set the initial values of parameters.

###  Using the fitter:

### 1.Adding the data points:

#### 1.1 To store or not to store the input data?
     - There are 2 options in the constructor - to store or not
       store the input data. The advantages of storing the data
       are that you'll be able to reset the fitting model without
       adding all the points again, and that for very large sets
       of points the chisquare is calculated more precisely.
       The obvious disadvantage is the amount of memory used to
       keep all the points.
     - Before you start adding the points, you can change the
       store/not store option by StoreData() method.

#### 1.2 The data can be added:
     - simply point by point - AddPoint() method
     - an array of points at once:
       If the data is already stored in some arrays, this data
       can be assigned to the linear fitter without physically
       coping bytes, thanks to the Use() method of
       TVector and TMatrix classes - AssignData() method

### 2.Setting the formula

#### 2.1 The linear formula syntax:
     -Additive parts are separated by 2 plus signes "++"
      --for example "1 ++ x" - for fitting a straight line
     -All standard functions, undrestood by TFormula, can be used
      as additive parts
      --TMath functions can be used too
     -Functions, used as additive parts, shouldn't have any parameters,
      even if those parameters are set.
      --for example, if normalizing a sum of a gaus(0, 1) and a
        gaus(0, 2), don't use the built-in "gaus" of TFormula,
        because it has parameters, take TMath::Gaus(x, 0, 1) instead.
     -Polynomials can be used like "pol3", .."polN"
     -If fitting a more than 3-dimensional formula, variables should
      be numbered as follows:
      -- x[0], x[1], x[2]... For example, to fit  "1 ++ x[0] ++ x[1] ++ x[2] ++ x[3]*x[3]"

#### 2.2 Setting the formula:

#####   2.2.1 If fitting a 1-2-3-dimensional formula, one can create a
         TF123 based on a linear expression and pass this function
         to the fitter:
         --Example:
~~~~
           TLinearFitter *lf = new TLinearFitter();
           TF2 *f2 = new TF2("f2", "x ++ y ++ x*x*y*y", -2, 2, -2, 2);
           lf->SetFormula(f2);
~~~~
         --The results of the fit are then stored in the function,
           just like when the TH1::Fit or TGraph::Fit is used
         --A linear function of this kind is by no means different
           from any other function, it can be drawn, evaluated, etc.

         --For multidimensional fitting, TFormulas of the form:
           x[0]++...++x[n] can be used
#####   2.2.2 There is no need to create the function if you don't want to,
         the formula can be set by expression:
         --Example:
~~~~
           // 2 is the number of dimensions
           TLinearFitter *lf = new TLinearFitter(2);
           lf->SetFormula("x ++ y ++ x*x*y*y");
~~~~

#####   2.2.3 The fastest functions to compute are polynomials and hyperplanes.
         --Polynomials are set the usual way: "pol1", "pol2",...
         --Hyperplanes are set by expression "hyp3", "hyp4", ...
         ---The "hypN" expressions only work when the linear fitter
            is used directly, not through TH1::Fit or TGraph::Fit.
            To fit a graph or a histogram with a hyperplane, define
            the function as "1++x++y".
         ---A constant term is assumed for a hyperplane, when using
            the "hypN" expression, so "hyp3" is in fact fitting with
            "1++x++y++z" function.
         --Fitting hyperplanes is much faster than fitting other
           expressions so if performance is vital, calculate the
           function values beforehand and give them to the fitter
           as variables
         --Example:
           You want to fit "sin(x)|cos(2*x)" very fast. Calculate
           sin(x) and cos(2*x) beforehand and store them in array *data.
           Then:
           TLinearFitter *lf=new TLinearFitter(2, "hyp2");
           lf->AssignData(npoint, 2, data, y);

#### 2.3 Resetting the formula

#####   2.3.1 If the input data is stored (or added via AssignData() function),
         the fitting formula can be reset without re-adding all the points.
         --Example:
~~~~
           TLinearFitter *lf=new TLinearFitter("1++x++x*x");
           lf->AssignData(n, 1, x, y, e);
           lf->Eval()
           //looking at the parameter significance, you see,
           // that maybe the fit will improve, if you take out
           // the constant term
           lf->SetFormula("x++x*x");
           lf->Eval();
           ...
~~~~

#####   2.3.2 If the input data is not stored, the fitter will have to be
         cleared and the data will have to be added again to try a
         different formula.

### 3.Accessing the fit results

#### 3.1 There are methods in the fitter to access all relevant information:
     --GetParameters, GetCovarianceMatrix, etc
     --the t-values of parameters and their significance can be reached by
       GetParTValue() and GetParSignificance() methods

#### 3.2 If fitting with a pre-defined TF123, the fit results are also
     written into this function.


### 4.Robust fitting - Least Trimmed Squares regression (LTS)
  Outliers are atypical(by definition), infrequant observations; data points
  which do not appear to follow the characteristic distribution of the rest
  of the data. These may reflect genuine properties of the underlying
  phenomenon(variable), or be due to measurement errors or anomalies which
  shouldn't be modelled. (StatSoft electronic textbook)

  Even a single gross outlier can greatly influence the results of least-
  squares fitting procedure, and in this case use of robust(resistant) methods
  is recommended.

  The method implemented here is based on the article and algorithm:
  "Computing LTS Regression for Large Data Sets" by
  P.J.Rousseeuw and Katrien Van Driessen
  The idea of the method is to find the fitting coefficients for a subset
  of h observations (out of n) with the smallest sum of squared residuals.
  The size of the subset h should lie between (npoints + nparameters +1)/2
  and n, and represents the minimal number of good points in the dataset.
  The default value is set to (npoints + nparameters +1)/2, but of course
  if you are sure that the data contains less outliers it's better to change
  h according to your data.

  To perform a robust fit, call EvalRobust() function instead of Eval() after
  adding the points and setting the fitting function.
  Note, that standard errors on parameters are not computed!

*/



////////////////////////////////////////////////////////////////////////////////
///default c-tor, input data is stored
///If you don't want to store the input data,
///run the function StoreData(kFALSE) after constructor

TLinearFitter::TLinearFitter() :
TVirtualFitter(),
   fParams(),
   fParCovar(),
   fTValues(),
   fDesign(),
   fDesignTemp(),
   fDesignTemp2(),
   fDesignTemp3(),
   fAtb(),
   fAtbTemp(),
   fAtbTemp2(),
   fAtbTemp3(),
   fFunctions(),
   fY(),
   fX(),
   fE(),
   fVal()
{
   fChisquare =0;
   fNpoints   =0;
   fNdim      =0;
   fY2        =0;
   fY2Temp    =0;
   fNfixed    =0;
   fIsSet     =kFALSE;
   fFormula   =0;
   fFixedParams=0;
   fSpecial   =0;
   fInputFunction=0;
   fStoreData =kTRUE;
   fRobust    =kFALSE;
   fNfunctions = 0;
   fFormulaSize = 0;
   fH = 0;
}

////////////////////////////////////////////////////////////////////////////////
///The parameter stands for number of dimensions in the fitting formula
///The input data is stored. If you don't want to store the input data,
///run the function StoreData(kFALSE) after constructor

TLinearFitter::TLinearFitter(Int_t ndim) :
   fVal()
{
   fNdim    =ndim;
   fNpoints =0;
   fY2      =0;
   fY2Temp  =0;
   fNfixed  =0;
   fFixedParams=0;
   fFormula =0;
   fIsSet   =kFALSE;
   fChisquare=0;
   fSpecial  =0;
   fInputFunction=0;
   fStoreData=kTRUE;
   fRobust = kFALSE;
   fNfunctions = 0;
   fFormulaSize = 0;
   fH = 0;
}

////////////////////////////////////////////////////////////////////////////////
///First parameter stands for number of dimensions in the fitting formula
///Second parameter is the fitting formula: see class description for formula syntax
///Options:
///The option is to store or not to store the data
///If you don't want to store the data, choose "" for the option, or run
///StoreData(kFalse) member function after the constructor

TLinearFitter::TLinearFitter(Int_t ndim, const char *formula, Option_t *opt)
{
   fNdim=ndim;
   fNpoints=0;
   fChisquare=0;
   fY2=0;
   fNfixed=0;
   fFixedParams=0;
   fSpecial=0;
   fInputFunction=0;
   fFormula = 0;
   TString option=opt;
   option.ToUpper();
   if (option.Contains("D"))
      fStoreData=kTRUE;
   else
      fStoreData=kFALSE;
   fRobust=kFALSE;
   SetFormula(formula);
}

////////////////////////////////////////////////////////////////////////////////
///This constructor uses a linear function. How to create it?
///TFormula now accepts formulas of the following kind:
///TFormula("f", "x++y++z++x*x") or
///TFormula("f", "x[0]++x[1]++x[2]*x[2]");
///Other than the look, it's in no
///way different from the regular formula, it can be evaluated,
///drawn, etc.
///The option is to store or not to store the data
///If you don't want to store the data, choose "" for the option, or run
///StoreData(kFalse) member function after the constructor

TLinearFitter::TLinearFitter(TFormula *function, Option_t *opt)
{
   fNdim=function->GetNdim();
   if (!function->IsLinear()){
      Int_t number=function->GetNumber();
      if (number<299 || number>310){
         Error("TLinearFitter", "Trying to fit with a nonlinear function");
         return;
      }
   }
   fNpoints=0;
   fChisquare=0;
   fY2=0;
   fNfixed=0;
   fFixedParams=0;
   fSpecial=0;
   fFormula = 0;
   TString option=opt;
   option.ToUpper();
   if (option.Contains("D"))
      fStoreData=kTRUE;
   else
      fStoreData=kFALSE;
   fIsSet=kTRUE;
   fRobust=kFALSE;
   fInputFunction=0;

   SetFormula(function);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor

TLinearFitter::TLinearFitter(const TLinearFitter& tlf) :
   TVirtualFitter(tlf),
   fParams(tlf.fParams),
   fParCovar(tlf.fParCovar),
   fTValues(tlf.fTValues),
   fParSign(tlf.fParSign),
   fDesign(tlf.fDesign),
   fDesignTemp(tlf.fDesignTemp),
   fDesignTemp2(tlf.fDesignTemp2),
   fDesignTemp3(tlf.fDesignTemp3),
   fAtb(tlf.fAtb),
   fAtbTemp(tlf.fAtbTemp),
   fAtbTemp2(tlf.fAtbTemp2),
   fAtbTemp3(tlf.fAtbTemp3),
   fFunctions( * (TObjArray *)tlf.fFunctions.Clone()),
   fY(tlf.fY),
   fY2(tlf.fY2),
   fY2Temp(tlf.fY2Temp),
   fX(tlf.fX),
   fE(tlf.fE),
   fInputFunction(tlf.fInputFunction),
   fVal(),
   fNpoints(tlf.fNpoints),
   fNfunctions(tlf.fNfunctions),
   fFormulaSize(tlf.fFormulaSize),
   fNdim(tlf.fNdim),
   fNfixed(tlf.fNfixed),
   fSpecial(tlf.fSpecial),
   fFormula(0),
   fIsSet(tlf.fIsSet),
   fStoreData(tlf.fStoreData),
   fChisquare(tlf.fChisquare),
   fH(tlf.fH),
   fRobust(tlf.fRobust),
   fFitsample(tlf.fFitsample),
   fFixedParams(0)
{
   // make a deep  copy of managed objects
   // fFormula, fFixedParams and fFunctions

   if ( tlf.fFixedParams && fNfixed > 0 ) {
      fFixedParams=new Bool_t[fNfixed];
      for(Int_t i=0; i<fNfixed; ++i)
         fFixedParams[i]=tlf.fFixedParams[i];
   }
   if (tlf.fFormula) {
      fFormula = new char[fFormulaSize+1];
      strlcpy(fFormula,tlf.fFormula,fFormulaSize+1);
   }

}


////////////////////////////////////////////////////////////////////////////////
/// Linear fitter cleanup.

TLinearFitter::~TLinearFitter()
{
   if (fFormula) {
      delete [] fFormula;
      fFormula = 0;
   }
   if (fFixedParams) {
      delete [] fFixedParams;
      fFixedParams = 0;
   }
   fInputFunction = 0;

   //fFunctions.Delete();
   fFunctions.Clear();

}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TLinearFitter& TLinearFitter::operator=(const TLinearFitter& tlf)
{
   if(this!=&tlf) {

      TVirtualFitter::operator=(tlf);
      fParams.ResizeTo(tlf.fParams);      fParams=tlf.fParams;
      fParCovar.ResizeTo(tlf.fParCovar);  fParCovar=tlf.fParCovar;
      fTValues.ResizeTo(tlf.fTValues);    fTValues=tlf.fTValues;
      fParSign.ResizeTo(tlf.fParSign);    fParSign=tlf.fParSign;
      fDesign.ResizeTo(tlf.fDesign);                fDesign=tlf.fDesign;
      fDesignTemp.ResizeTo(tlf.fDesignTemp);        fDesignTemp=tlf.fDesignTemp;
      fDesignTemp2.ResizeTo(tlf.fDesignTemp2);      fDesignTemp2=tlf.fDesignTemp2;
      fDesignTemp3.ResizeTo(tlf.fDesignTemp3);      fDesignTemp3=tlf.fDesignTemp3;

      fAtb.ResizeTo(tlf.fAtb);            fAtb=tlf.fAtb;
      fAtbTemp.ResizeTo(tlf.fAtbTemp);    fAtbTemp=tlf.fAtbTemp;
      fAtbTemp2.ResizeTo(tlf.fAtbTemp2);    fAtbTemp2=tlf.fAtbTemp2;
      fAtbTemp3.ResizeTo(tlf.fAtbTemp3);    fAtbTemp3=tlf.fAtbTemp3;

      // use clear instead of delete
      fFunctions.Clear();
      fFunctions= *(TObjArray*) tlf.fFunctions.Clone();

      fY.ResizeTo(tlf.fY);    fY = tlf.fY;
      fX.ResizeTo(tlf.fX);    fX = tlf.fX;
      fE.ResizeTo(tlf.fE);    fE = tlf.fE;

      fY2 = tlf.fY2;
      fY2Temp = tlf.fY2Temp;
      for(Int_t i = 0; i < 1000; i++) fVal[i] = tlf.fVal[i];

      if(fInputFunction) { delete fInputFunction; fInputFunction = 0; }
      if(tlf.fInputFunction) fInputFunction = new TFormula(*tlf.fInputFunction);

      fNpoints=tlf.fNpoints;
      fNfunctions=tlf.fNfunctions;
      fFormulaSize=tlf.fFormulaSize;
      fNdim=tlf.fNdim;
      fNfixed=tlf.fNfixed;
      fSpecial=tlf.fSpecial;

      if(fFormula) { delete [] fFormula; fFormula = 0; }
      if (tlf.fFormula) {
         fFormula = new char[fFormulaSize+1];
         strlcpy(fFormula,tlf.fFormula,fFormulaSize+1);
      }

      fIsSet=tlf.fIsSet;
      fStoreData=tlf.fStoreData;
      fChisquare=tlf.fChisquare;

      fH=tlf.fH;
      fRobust=tlf.fRobust;
      fFitsample=tlf.fFitsample;

      if(fFixedParams) { delete [] fFixedParams; fFixedParams = 0; }
      if ( tlf.fFixedParams && fNfixed > 0 ) {
         fFixedParams=new Bool_t[fNfixed];
         for(Int_t i=0; i< fNfixed; ++i)
            fFixedParams[i]=tlf.fFixedParams[i];
      }

   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
///Add another linear fitter to this linear fitter. Points and Design matrices
///are added, but the previos fitting results (if any) are deleted.
///Fitters must have same formulas (this is not checked). Fixed parameters are not changed

void TLinearFitter::Add(TLinearFitter *tlf)
{
   fParams.Zero();
   fParCovar.Zero();
   fTValues.Zero();
   fParSign.Zero();

   fDesign += tlf->fDesign;

   fDesignTemp += tlf->fDesignTemp;
   fDesignTemp2 += tlf->fDesignTemp2;
   fDesignTemp3 += tlf->fDesignTemp3;
   fAtb += tlf->fAtb;
   fAtbTemp += tlf->fAtbTemp;
   fAtbTemp2 += tlf->fAtbTemp2;
   fAtbTemp3 += tlf->fAtbTemp3;

   if (fStoreData){
      Int_t size=fY.GetNoElements();
      Int_t newsize = fNpoints+tlf->fNpoints;
      if (size<newsize){
         fY.ResizeTo(newsize);
         fE.ResizeTo(newsize);
         fX.ResizeTo(newsize, fNdim);
      }
      for (Int_t i=fNpoints; i<newsize; i++){
         fY(i)=tlf->fY(i-fNpoints);
         fE(i)=tlf->fE(i-fNpoints);
         for (Int_t j=0; j<fNdim; j++){
            fX(i,j)=tlf->fX(i-fNpoints, j);
         }
      }
   }
   fY2 += tlf->fY2;
   fY2Temp += tlf->fY2Temp;
   fNpoints += tlf->fNpoints;
   //fInputFunction=(TFormula*)tlf.fInputFunction->Clone();

   fChisquare=0;
   fH=0;
   fRobust=0;
}

////////////////////////////////////////////////////////////////////////////////
///Adds 1 point to the fitter.
///First parameter stands for the coordinates of the point, where the function is measured
///Second parameter - the value being fitted
///Third parameter - weight(measurement error) of this point (=1 by default)

void TLinearFitter::AddPoint(Double_t *x, Double_t y, Double_t e)
{
   Int_t size;
   fNpoints++;
   if (fStoreData){
      size=fY.GetNoElements();
      if (size<fNpoints){
         fY.ResizeTo(fNpoints+fNpoints/2);
         fE.ResizeTo(fNpoints+fNpoints/2);
         fX.ResizeTo(fNpoints+fNpoints/2, fNdim);
      }

      Int_t j=fNpoints-1;
      fY(j)=y;
      fE(j)=e;
      for (Int_t i=0; i<fNdim; i++)
         fX(j,i)=x[i];
   }
   //add the point to the design matrix, if the formula has been set
   // (LM: why condition !fRobust ?? - remove it to fix Coverity 11602)
   // when 100< fSpecial < 200 (Polynomial) fFunctions is not empty
   // (see implementation of SetFormula method)
   if (fFunctions.IsEmpty()&&(!fInputFunction)&&(fSpecial<=200)){
      Error("AddPoint", "Point can't be added, because the formula hasn't been set");
      return;
   }
   if (!fRobust) AddToDesign(x, y, e);
}

////////////////////////////////////////////////////////////////////////////////
///This function is to use when you already have all the data in arrays
///and don't want to copy them into the fitter. In this function, the Use() method
///of TVectorD and TMatrixD is used, so no bytes are physically moved around.
///First parameter - number of points to fit
///Second parameter - number of variables in the model
///Third parameter - the variables of the model, stored in the following way:
///(x0(0), x1(0), x2(0), x3(0), x0(1), x1(1), x2(1), x3(1),...

void TLinearFitter::AssignData(Int_t npoints, Int_t xncols, Double_t *x, Double_t *y, Double_t *e)
{
   if (npoints<fNpoints){
      Error("AddData", "Those points are already added");
      return;
   }
   Bool_t same=kFALSE;
   if (fX.GetMatrixArray()==x && fY.GetMatrixArray()==y){
      if (e){
         if (fE.GetMatrixArray()==e)
            same=kTRUE;
      }
   }

   fX.Use(npoints, xncols, x);
   fY.Use(npoints, y);
   if (e)
      fE.Use(npoints, e);
   else {
      fE.ResizeTo(npoints);
      fE=1;
   }
   Int_t xfirst;
   if (!fFunctions.IsEmpty() || fInputFunction ||  fSpecial>200) {
      if (same)
         xfirst=fNpoints;

      else
         xfirst=0;
      for (Int_t i=xfirst; i<npoints; i++)
         AddToDesign(TMatrixDRow(fX, i).GetPtr(), fY(i), fE(i));
   }
   fNpoints=npoints;
}

////////////////////////////////////////////////////////////////////////////////
///Add a point to the AtA matrix and to the Atb vector.

void TLinearFitter::AddToDesign(Double_t *x, Double_t y, Double_t e)
{


   Int_t i, j, ii;
   y/=e;

//   Double_t fVal[1000];

   if ((fSpecial>100)&&(fSpecial<200)){
      //polynomial fitting
      Int_t npar=fSpecial-100;
      fVal[0]=1;
      for (i=1; i<npar; i++)
         fVal[i]=fVal[i-1]*x[0];
      for (i=0; i<npar; i++)
         fVal[i]/=e;
   } else if (fSpecial>200){
      //Hyperplane fitting. Constant term is added
      Int_t npar=fSpecial-201;
      fVal[0]=1./e;
      for (i=0; i<npar; i++)
         fVal[i+1]=x[i]/e;
   } else {
      //general case
      for (ii=0; ii<fNfunctions; ii++){
         if (!fFunctions.IsEmpty()){
            // ffunctions can be TF1 or TFormula depending on how they are created
            TObject * obj = fFunctions.UncheckedAt(ii);
            if (obj->IsA() == TFormula::Class() ) {
               TFormula *f1 = (TFormula*)(obj);
               fVal[ii]=f1->EvalPar(x)/e;
            }
            else if  (obj->IsA() == TF1::Class() ) {
               TF1 *f1 = (TF1*)(obj);
               fVal[ii]=f1->EvalPar(x)/e;
            }
            else {
               Error("AddToDesign","Basis Function %s is of an invalid type %s",obj->GetName(),obj->IsA()->GetName());
               return;
            }
         } else {
            TFormula *f=(TFormula*)fInputFunction->GetLinearPart(ii);
            if (!f){
               Error("AddToDesign","Function %s has no linear parts - maybe missing a ++ in the formula expression",fInputFunction->GetName());
               return;
            }
            fVal[ii]=f->EvalPar(x)/e;
         }
      }
   }
   //additional matrices for numerical stability
   for (i=0; i<fNfunctions; i++){
      for (j=0; j<i; j++)
         fDesignTemp3(j, i)+=fVal[i]*fVal[j];
      fDesignTemp3(i, i)+=fVal[i]*fVal[i];
      fAtbTemp3(i)+=fVal[i]*y;

   }
   fY2Temp+=y*y;
   fIsSet=kTRUE;

   if (fNpoints % 100 == 0 && fNpoints>100){
      fDesignTemp2+=fDesignTemp3;
      fDesignTemp3.Zero();
      fAtbTemp2+=fAtbTemp3;
      fAtbTemp3.Zero();
      if (fNpoints % 10000 == 0 && fNpoints>10000){
         fDesignTemp+=fDesignTemp2;
         fDesignTemp2.Zero();
         fAtbTemp+=fAtbTemp2;
         fAtbTemp2.Zero();
         fY2+=fY2Temp;
         fY2Temp=0;
         if (fNpoints % 1000000 == 0 && fNpoints>1000000){
            fDesign+=fDesignTemp;
            fDesignTemp.Zero();
            fAtb+=fAtbTemp;
            fAtbTemp.Zero();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TLinearFitter::AddTempMatrices()
{
   if (fDesignTemp3.GetNrows()>0){
      fDesignTemp2+=fDesignTemp3;
      fDesignTemp+=fDesignTemp2;
      fDesign+=fDesignTemp;
      fDesignTemp3.Zero();
      fDesignTemp2.Zero();
      fDesignTemp.Zero();
      fAtbTemp2+=fAtbTemp3;
      fAtbTemp+=fAtbTemp2;
      fAtb+=fAtbTemp;
      fAtbTemp3.Zero();
      fAtbTemp2.Zero();
      fAtbTemp.Zero();

      fY2+=fY2Temp;
      fY2Temp=0;
      }
}

////////////////////////////////////////////////////////////////////////////////
///Clears everything. Used in TH1::Fit and TGraph::Fit().

void TLinearFitter::Clear(Option_t * /*option*/)
{
   fParams.Clear();
   fParCovar.Clear();
   fTValues.Clear();
   fParSign.Clear();
   fDesign.Clear();
   fDesignTemp.Clear();
   fDesignTemp2.Clear();
   fDesignTemp3.Clear();
   fAtb.Clear();
   fAtbTemp.Clear();
   fAtbTemp2.Clear();
   fAtbTemp3.Clear();
   fFunctions.Clear();
   fInputFunction=0;
   fY.Clear();
   fX.Clear();
   fE.Clear();

   fNpoints=0;
   fNfunctions=0;
   fFormulaSize=0;
   fNdim=0;
   if (fFormula) delete [] fFormula;
   fFormula=0;
   fIsSet=0;
   if (fFixedParams) delete [] fFixedParams;
   fFixedParams=0;

   fChisquare=0;
   fY2=0;
   fSpecial=0;
   fRobust=kFALSE;
   fFitsample.Clear();
}

////////////////////////////////////////////////////////////////////////////////
///To be used when different sets of points are fitted with the same formula.

void TLinearFitter::ClearPoints()
{
   fDesign.Zero();
   fAtb.Zero();
   fDesignTemp.Zero();
   fDesignTemp2.Zero();
   fDesignTemp3.Zero();
   fAtbTemp.Zero();
   fAtbTemp2.Zero();
   fAtbTemp3.Zero();

   fParams.Zero();
   fParCovar.Zero();
   fTValues.Zero();
   fParSign.Zero();

   for (Int_t i=0; i<fNfunctions; i++)
      fFixedParams[i]=0;
   fChisquare=0;
   fNpoints=0;

}

////////////////////////////////////////////////////////////////////////////////
///Calculates the chisquare.

void TLinearFitter::Chisquare()
{
   Int_t i, j;
   Double_t sumtotal2;
   Double_t temp, temp2;

   if (!fStoreData){
      sumtotal2 = 0;
      for (i=0; i<fNfunctions; i++){
         for (j=0; j<i; j++){
            sumtotal2 += 2*fParams(i)*fParams(j)*fDesign(j, i);
         }
         sumtotal2 += fParams(i)*fParams(i)*fDesign(i, i);
         sumtotal2 -= 2*fParams(i)*fAtb(i);
      }
      sumtotal2 += fY2;
   } else {
      sumtotal2 = 0;
      if (fInputFunction){
         for (i=0; i<fNpoints; i++){
            temp = fInputFunction->EvalPar(TMatrixDRow(fX, i).GetPtr());
            temp2 = (fY(i)-temp)*(fY(i)-temp);
            temp2 /= fE(i)*fE(i);
            sumtotal2 += temp2;
         }
      } else {
         sumtotal2 = 0;
         Double_t val[100];
         for (Int_t point=0; point<fNpoints; point++){
            temp = 0;
            if ((fSpecial>100)&&(fSpecial<200)){
               Int_t npar = fSpecial-100;
               val[0] = 1;
               for (i=1; i<npar; i++)
                  val[i] = val[i-1]*fX(point, 0);
               for (i=0; i<npar; i++)
                  temp += fParams(i)*val[i];
            } else {
               if (fSpecial>200) {
                  //hyperplane case
                  Int_t npar = fSpecial-201;
                  temp+=fParams(0);
                  for (i=0; i<npar; i++)
                     temp += fParams(i+1)*fX(point, i);
               } else {
                  for (j=0; j<fNfunctions; j++) {
                     TFormula *f1 = (TFormula*)(fFunctions.UncheckedAt(j));
                     val[j] = f1->EvalPar(TMatrixDRow(fX, point).GetPtr());
                     temp += fParams(j)*val[j];
                  }
               }
            }
         temp2 = (fY(point)-temp)*(fY(point)-temp);
         temp2 /= fE(point)*fE(point);
         sumtotal2 += temp2;
         }
      }
   }
   fChisquare = sumtotal2;

}

////////////////////////////////////////////////////////////////////////////////
/// Computes parameters' t-values and significance

void TLinearFitter::ComputeTValues()
{
   for (Int_t i=0; i<fNfunctions; i++){
      fTValues(i) = fParams(i)/(TMath::Sqrt(fParCovar(i, i)));
      fParSign(i) = 2*(1-TMath::StudentI(TMath::Abs(fTValues(i)),fNpoints-fNfunctions+fNfixed));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the fit and evaluate the parameters
/// Returns 0 if the fit is ok, 1 if there are errors

Int_t TLinearFitter::Eval()
{
   Double_t e;
   if (fFunctions.IsEmpty()&&(!fInputFunction)&&(fSpecial<=200)){
      Error("TLinearFitter::Eval", "The formula hasn't been set");
      return 1;
   }
   //
   fParams.ResizeTo(fNfunctions);
   fTValues.ResizeTo(fNfunctions);
   fParSign.ResizeTo(fNfunctions);
   fParCovar.ResizeTo(fNfunctions,fNfunctions);

   fChisquare=0;

   if (!fIsSet){
      Bool_t update = UpdateMatrix();
      if (!update){
         //no points to fit
         fParams.Zero();
         fParCovar.Zero();
         fTValues.Zero();
         fParSign.Zero();
         fChisquare=0;
         if (fInputFunction){
            fInputFunction->SetParameters(fParams.GetMatrixArray());
            for (Int_t i=0; i<fNfunctions; i++)
               ((TF1*)fInputFunction)->SetParError(i, 0);
            ((TF1*)fInputFunction)->SetChisquare(0);
            ((TF1*)fInputFunction)->SetNDF(0);
            ((TF1*)fInputFunction)->SetNumberFitPoints(0);
         }
         return 1;
      }
   }

   AddTempMatrices();

//fixing fixed parameters, if there are any
   Int_t i, ii, j=0;
   if (fNfixed>0){
      for (ii=0; ii<fNfunctions; ii++)
         fDesignTemp(ii, fNfixed) = fAtb(ii);
      for (i=0; i<fNfunctions; i++){
         if (fFixedParams[i]){
            for (ii=0; ii<i; ii++)
               fDesignTemp(ii, j) = fDesign(ii, i);
            for (ii=i; ii<fNfunctions; ii++)
               fDesignTemp(ii, j) = fDesign(i, ii);
            j++;
            for (ii=0; ii<fNfunctions; ii++){
               fAtb(ii)-=fParams(i)*(fDesignTemp(ii, j-1));
            }
         }
      }
      for (i=0; i<fNfunctions; i++){
         if (fFixedParams[i]){
            for (ii=0; ii<fNfunctions; ii++){
               fDesign(ii, i) = 0;
               fDesign(i, ii) = 0;
            }
            fDesign (i, i) = 1;
            fAtb(i) = fParams(i);
         }
      }
   }

   TDecompChol chol(fDesign);
   Bool_t ok;
   TVectorD coef(fNfunctions);
   coef=chol.Solve(fAtb, ok);
   if (!ok){
      Error("Eval","Matrix inversion failed");
      fParams.Zero();
      fParCovar.Zero();
      fTValues.Zero();
      fParSign.Zero();
      if (fInputFunction){
         fInputFunction->SetParameters(fParams.GetMatrixArray());
         if (fInputFunction->InheritsFrom(TF1::Class())){
            ((TF1*)fInputFunction)->SetChisquare(0);
            ((TF1*)fInputFunction)->SetNDF(fNpoints-fNfunctions+fNfixed);
            ((TF1*)fInputFunction)->SetNumberFitPoints(fNpoints);
         }
      }
      return 1;
   }
   fParams=coef;
   fParCovar=chol.Invert();

   if (fInputFunction){
      fInputFunction->SetParameters(fParams.GetMatrixArray());
      if (fInputFunction->InheritsFrom(TF1::Class())){
         for (i=0; i<fNfunctions; i++){
            e = TMath::Sqrt(fParCovar(i, i));
            ((TF1*)fInputFunction)->SetParError(i, e);
         }
         if (!fObjectFit)
            ((TF1*)fInputFunction)->SetChisquare(GetChisquare());
         ((TF1*)fInputFunction)->SetNDF(fNpoints-fNfunctions+fNfixed);
         ((TF1*)fInputFunction)->SetNumberFitPoints(fNpoints);
         }
   }

   //if parameters were fixed, change the design matrix back as it was before fixing
   j = 0;
   if (fNfixed>0){
      for (i=0; i<fNfunctions; i++){
         if (fFixedParams[i]){
            for (ii=0; ii<i; ii++){
               fDesign(ii, i) = fDesignTemp(ii, j);
               fAtb(ii) = fDesignTemp(ii, fNfixed);
            }
            for (ii=i; ii<fNfunctions; ii++){
               fDesign(i, ii) = fDesignTemp(ii, j);
               fAtb(ii) = fDesignTemp(ii, fNfixed);
            }
            j++;
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///Fixes paramter #ipar at its current value.

void TLinearFitter::FixParameter(Int_t ipar)
{
   if (fParams.NonZeros()<1){
      Error("FixParameter", "no value available to fix the parameter");
      return;
   }
   if (ipar>fNfunctions || ipar<0){
      Error("FixParameter", "illegal parameter value");
      return;
   }
   if (fNfixed==fNfunctions) {
      Error("FixParameter", "no free parameters left");
      return;
   }
   if (!fFixedParams)
      fFixedParams = new Bool_t[fNfunctions];
   fFixedParams[ipar] = 1;
   fNfixed++;
}

////////////////////////////////////////////////////////////////////////////////
///Fixes parameter #ipar at value parvalue.

void TLinearFitter::FixParameter(Int_t ipar, Double_t parvalue)
{
   if (ipar>fNfunctions || ipar<0){
      Error("FixParameter", "illegal parameter value");
      return;
   }
   if (fNfixed==fNfunctions) {
      Error("FixParameter", "no free parameters left");
      return;
   }
   if(!fFixedParams)
      fFixedParams = new Bool_t[fNfunctions];
   fFixedParams[ipar] = 1;
   if (fParams.GetNoElements()<fNfunctions)
      fParams.ResizeTo(fNfunctions);
   fParams(ipar) = parvalue;
   fNfixed++;
}

////////////////////////////////////////////////////////////////////////////////
///Releases parameter #ipar.

void TLinearFitter::ReleaseParameter(Int_t ipar)
{
   if (ipar>fNfunctions || ipar<0){
      Error("ReleaseParameter", "illegal parameter value");
      return;
   }
   if (!fFixedParams[ipar]){
      Warning("ReleaseParameter","This parameter is not fixed\n");
      return;
   } else {
      fFixedParams[ipar] = 0;
      fNfixed--;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Get the Atb vector - a vector, used for internal computations

void TLinearFitter::GetAtbVector(TVectorD &v)
{
   if (v.GetNoElements()!=fAtb.GetNoElements())
      v.ResizeTo(fAtb.GetNoElements());
   v = fAtb;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the Chisquare.

Double_t TLinearFitter::GetChisquare()
{
   if (fChisquare > 1e-16)
      return fChisquare;
   else {
      Chisquare();
      return fChisquare;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Computes point-by-point confidence intervals for the fitted function
///Parameters:
///n - number of points
///ndim - dimensions of points
///x - points, at which to compute the intervals, for ndim > 1
///    should be in order: (x0,y0, x1, y1, ... xn, yn)
///ci - computed intervals are returned in this array
///cl - confidence level, default=0.95
///
///NOTE, that this method can only be used when the fitting function inherits from a TF1,
///so it's not possible when the fitting function was set as a string or as a pure TFormula

void TLinearFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)
{
   if (fInputFunction){
      Double_t *grad = new Double_t[fNfunctions];
      Double_t *sum_vector = new Double_t[fNfunctions];
      Double_t c=0;
      Int_t df = fNpoints-fNfunctions+fNfixed;
      Double_t t = TMath::StudentQuantile(0.5 + cl/2, df);
      Double_t chidf = TMath::Sqrt(fChisquare/df);

      for (Int_t ipoint=0; ipoint<n; ipoint++){
         c=0;
         ((TF1*)(fInputFunction))->GradientPar(x+ndim*ipoint, grad);
         //multiply the covariance matrix by gradient
         for (Int_t irow=0; irow<fNfunctions; irow++){
            sum_vector[irow]=0;
            for (Int_t icol=0; icol<fNfunctions; icol++){
               sum_vector[irow]+=fParCovar(irow,icol)*grad[icol];
            }
         }
         for (Int_t i=0; i<fNfunctions; i++)
            c+=grad[i]*sum_vector[i];
         c=TMath::Sqrt(c);
         ci[ipoint]=c*t*chidf;
      }

      delete [] grad;
      delete [] sum_vector;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Computes confidence intervals at level cl. Default is 0.95
///The TObject parameter can be a TGraphErrors, a TGraph2DErrors or a TH123.
///For Graphs, confidence intervals are computed for each point,
///the value of the graph at that point is set to the function value at that
///point, and the graph y-errors (or z-errors) are set to the value of
///the confidence interval at that point
///For Histograms, confidence intervals are computed for each bin center
///The bin content of this bin is then set to the function value at the bin
///center, and the bin error is set to the confidence interval value.
///Allowed combinations:
///Fitted object               Passed object
///TGraph                      TGraphErrors, TH1
///TGraphErrors, AsymmErrors   TGraphErrors, TH1
///TH1                         TGraphErrors, TH1
///TGraph2D                    TGraph2DErrors, TH2
///TGraph2DErrors              TGraph2DErrors, TH2
///TH2                         TGraph2DErrors, TH2
///TH3                         TH3

void TLinearFitter::GetConfidenceIntervals(TObject *obj, Double_t cl)
{
   if (!fInputFunction) {
      Error("GetConfidenceIntervals", "The case of fitting not with a TFormula is not yet implemented");
      return;
   }

   //TGraph//////////////////

   if (obj->InheritsFrom(TGraph::Class())) {
      TGraph *gr = (TGraph*)obj;
      if (!gr->GetEY()){
         Error("GetConfidenceIntervals", "A TGraphErrors should be passed instead of a graph");
         return;
      }
      if (fObjectFit->InheritsFrom(TGraph2D::Class())){
         Error("GetConfidenceIntervals", "A TGraph2DErrors should be passed instead of a graph");
         return;
      }
      if (fObjectFit->InheritsFrom(TH1::Class())){
         if (((TH1*)(fObjectFit))->GetDimension()>1){
            Error("GetConfidenceIntervals", "A TGraph2DErrors or a TH23 should be passed instead of a graph");
            return;
         }
      }

      GetConfidenceIntervals(gr->GetN(),1,gr->GetX(), gr->GetEY(), cl);
      for (Int_t i=0; i<gr->GetN(); i++)
         gr->SetPoint(i, gr->GetX()[i], fInputFunction->Eval(gr->GetX()[i]));
   }

   //TGraph2D///////////////
   else if (obj->InheritsFrom(TGraph2D::Class())) {
      TGraph2D *gr2 = (TGraph2D*)obj;
      if (!gr2->GetEZ()){
         Error("GetConfidenceIntervals", "A TGraph2DErrors should be passed instead of a TGraph2D");
         return;
      }
      if (fObjectFit->InheritsFrom(TGraph::Class())){
         Error("GetConfidenceIntervals", "A TGraphErrors should be passed instead of a TGraph2D");
         return;
      }
      if (fObjectFit->InheritsFrom(TH1::Class())){
         if (((TH1*)(fObjectFit))->GetDimension()==1){
            Error("GetConfidenceIntervals", "A TGraphErrors or a TH1 should be passed instead of a graph");
            return;
         }
      }
      Double_t xy[2];
      Int_t np = gr2->GetN();
      Double_t *grad = new Double_t[fNfunctions];
      Double_t *sum_vector = new Double_t[fNfunctions];
      Double_t *x = gr2->GetX();
      Double_t *y = gr2->GetY();
      Double_t t = TMath::StudentQuantile(0.5 + cl/2, ((TF1*)(fInputFunction))->GetNDF());
      Double_t chidf = TMath::Sqrt(fChisquare/((TF1*)(fInputFunction))->GetNDF());
      Double_t c = 0;
      for (Int_t ipoint=0; ipoint<np; ipoint++){
         c=0;
         xy[0]=x[ipoint];
         xy[1]=y[ipoint];
         ((TF1*)(fInputFunction))->GradientPar(xy, grad);
         for (Int_t irow=0; irow<fNfunctions; irow++){
            sum_vector[irow]=0;
            for (Int_t icol=0; icol<fNfunctions; icol++)
               sum_vector[irow]+=fParCovar(irow, icol)*grad[icol];
         }
         for (Int_t i=0; i<fNfunctions; i++)
            c+=grad[i]*sum_vector[i];
         c=TMath::Sqrt(c);
         gr2->SetPoint(ipoint, xy[0], xy[1], fInputFunction->EvalPar(xy));
         gr2->GetEZ()[ipoint]=c*t*chidf;
      }
      delete [] grad;
      delete [] sum_vector;
   }

   //TH1////////////////////////
   else if (obj->InheritsFrom(TH1::Class())) {
      if (fObjectFit->InheritsFrom(TGraph::Class())){
         if (((TH1*)obj)->GetDimension()>1){
            Error("GetConfidenceIntervals", "Fitted graph and passed histogram have different number of dimensions");
            return;
         }
      }
      if (fObjectFit->InheritsFrom(TGraph2D::Class())){
         if (((TH1*)obj)->GetDimension()!=2){
            Error("GetConfidenceIntervals", "Fitted graph and passed histogram have different number of dimensions");
            return;
         }
      }
      if (fObjectFit->InheritsFrom(TH1::Class())){
         if (((TH1*)(fObjectFit))->GetDimension()!=((TH1*)(obj))->GetDimension()){
            Error("GetConfidenceIntervals", "Fitted and passed histograms have different number of dimensions");
            return;
         }
      }


      TH1 *hfit = (TH1*)obj;
      Double_t *grad = new Double_t[fNfunctions];
      Double_t *sum_vector = new Double_t[fNfunctions];
      Double_t x[3];
      Int_t hxfirst = hfit->GetXaxis()->GetFirst();
      Int_t hxlast  = hfit->GetXaxis()->GetLast();
      Int_t hyfirst = hfit->GetYaxis()->GetFirst();
      Int_t hylast  = hfit->GetYaxis()->GetLast();
      Int_t hzfirst = hfit->GetZaxis()->GetFirst();
      Int_t hzlast  = hfit->GetZaxis()->GetLast();

      TAxis *xaxis  = hfit->GetXaxis();
      TAxis *yaxis  = hfit->GetYaxis();
      TAxis *zaxis  = hfit->GetZaxis();
      Double_t t = TMath::StudentQuantile(0.5 + cl/2, ((TF1*)(fInputFunction))->GetNDF());
      Double_t chidf = TMath::Sqrt(fChisquare/((TF1*)(fInputFunction))->GetNDF());
      Double_t c=0;
      for (Int_t binz=hzfirst; binz<=hzlast; binz++){
         x[2]=zaxis->GetBinCenter(binz);
         for (Int_t biny=hyfirst; biny<=hylast; biny++) {
            x[1]=yaxis->GetBinCenter(biny);
            for (Int_t binx=hxfirst; binx<=hxlast; binx++) {
               x[0]=xaxis->GetBinCenter(binx);
               ((TF1*)(fInputFunction))->GradientPar(x, grad);
               c=0;
               for (Int_t irow=0; irow<fNfunctions; irow++){
                  sum_vector[irow]=0;
                  for (Int_t icol=0; icol<fNfunctions; icol++)
                     sum_vector[irow]+=fParCovar(irow, icol)*grad[icol];
               }
               for (Int_t i=0; i<fNfunctions; i++)
                  c+=grad[i]*sum_vector[i];
               c=TMath::Sqrt(c);
               hfit->SetBinContent(binx, biny, binz, fInputFunction->EvalPar(x));
               hfit->SetBinError(binx, biny, binz, c*t*chidf);
            }
         }
      }
      delete [] grad;
      delete [] sum_vector;
   }
   else {
      Error("GetConfidenceIntervals", "This object type is not supported");
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Returns covariance matrix

Double_t* TLinearFitter::GetCovarianceMatrix() const
{
   Double_t *p = const_cast<Double_t*>(fParCovar.GetMatrixArray());
   return p;
}

////////////////////////////////////////////////////////////////////////////////
///Returns covariance matrix

void TLinearFitter::GetCovarianceMatrix(TMatrixD &matr)
{
   if (matr.GetNrows()!=fNfunctions || matr.GetNcols()!=fNfunctions){
      matr.ResizeTo(fNfunctions, fNfunctions);
   }
   matr = fParCovar;
}

////////////////////////////////////////////////////////////////////////////////
///Returns the internal design matrix

void TLinearFitter::GetDesignMatrix(TMatrixD &matr)
{
   if (matr.GetNrows()!=fNfunctions || matr.GetNcols()!=fNfunctions){
      matr.ResizeTo(fNfunctions, fNfunctions);
   }
   matr = fDesign;
}

////////////////////////////////////////////////////////////////////////////////
///Returns parameter errors

void TLinearFitter::GetErrors(TVectorD &vpar)
{
   if (vpar.GetNoElements()!=fNfunctions) {
      vpar.ResizeTo(fNfunctions);
   }
   for (Int_t i=0; i<fNfunctions; i++)
      vpar(i) = TMath::Sqrt(fParCovar(i, i));

}

////////////////////////////////////////////////////////////////////////////////
///Returns parameter values

void TLinearFitter::GetParameters(TVectorD &vpar)
{
   if (vpar.GetNoElements()!=fNfunctions) {
      vpar.ResizeTo(fNfunctions);
   }
   vpar=fParams;
}

////////////////////////////////////////////////////////////////////////////////
///Returns the value and the name of the parameter #ipar
///NB: In the calling function the argument name must be set large enough

Int_t TLinearFitter::GetParameter(Int_t ipar,char* name,Double_t& value,Double_t& /*verr*/,Double_t& /*vlow*/, Double_t& /*vhigh*/) const
{
   if (ipar<0 || ipar>fNfunctions) {
      Error("GetParError", "illegal value of parameter");
      return 0;
   }
   value = fParams(ipar);
   if (fInputFunction)
      strcpy(name, fInputFunction->GetParName(ipar));
   else
      name = 0;
   return 1;
}


////////////////////////////////////////////////////////////////////////////////
///Returns the error of parameter #ipar

Double_t TLinearFitter::GetParError(Int_t ipar) const
{
   if (ipar<0 || ipar>fNfunctions) {
      Error("GetParError", "illegal value of parameter");
      return 0;
   }

   return TMath::Sqrt(fParCovar(ipar, ipar));
}


////////////////////////////////////////////////////////////////////////////////
///Returns name of parameter #ipar

const char *TLinearFitter::GetParName(Int_t ipar) const
{
   if (ipar<0 || ipar>fNfunctions) {
      Error("GetParError", "illegal value of parameter");
      return 0;
   }
   if (fInputFunction)
      return fInputFunction->GetParName(ipar);
   return "";
}

////////////////////////////////////////////////////////////////////////////////
///Returns the t-value for parameter #ipar

Double_t TLinearFitter::GetParTValue(Int_t ipar)
{
   if (ipar<0 || ipar>fNfunctions) {
      Error("GetParTValue", "illegal value of parameter");
      return 0;
   }
   if (!fTValues.NonZeros())
      ComputeTValues();
   return fTValues(ipar);
}

////////////////////////////////////////////////////////////////////////////////
///Returns the significance of parameter #ipar

Double_t TLinearFitter::GetParSignificance(Int_t ipar)
{
   if (ipar<0 || ipar>fNfunctions) {
      Error("GetParSignificance", "illegal value of parameter");
      return 0;
   }
   if (!fParSign.NonZeros())
      ComputeTValues();
   return fParSign(ipar);
}

////////////////////////////////////////////////////////////////////////////////
///For robust lts fitting, returns the sample, on which the best fit was based

void TLinearFitter::GetFitSample(TBits &bits)
{
   if (!fRobust){
      Error("GetFitSample", "there is no fit sample in ordinary least-squares fit");
      return;
   }
   for (Int_t i=0; i<fNpoints; i++)
      bits.SetBitNumber(i, fFitsample.TestBitNumber(i));

}

////////////////////////////////////////////////////////////////////////////////
///Merge objects in list

Int_t TLinearFitter::Merge(TCollection *list)
{
   if (!list) return -1;
   TIter next(list);
   TLinearFitter *lfit = 0;
   while ((lfit = (TLinearFitter*)next())) {
      if (!lfit->InheritsFrom(TLinearFitter::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s",lfit->ClassName(),this->ClassName());
         return -1;
      }
      Add(lfit);
   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
///set the basis functions in case the fitting function is not
/// set directly
/// The TLinearFitter will manage and delete the functions contained in the list

void TLinearFitter::SetBasisFunctions(TObjArray * functions)
{
   fFunctions = *(functions);
   fFunctions.SetOwner(kTRUE);
   int size = fFunctions.GetEntries();

   fNfunctions=size;
   //change the size of design matrix
   fDesign.ResizeTo(size, size);
   fAtb.ResizeTo(size);
   fDesignTemp.ResizeTo(size, size);
   fDesignTemp2.ResizeTo(size, size);
   fDesignTemp3.ResizeTo(size, size);
   fAtbTemp.ResizeTo(size);
   fAtbTemp2.ResizeTo(size);
   fAtbTemp3.ResizeTo(size);
   if (fFixedParams)
      delete [] fFixedParams;
   fFixedParams=new Bool_t[size];
   fDesign.Zero();
   fAtb.Zero();
   fDesignTemp.Zero();
   fDesignTemp2.Zero();
   fDesignTemp3.Zero();
   fAtbTemp.Zero();
   fAtbTemp2.Zero();
   fAtbTemp3.Zero();
   fY2Temp=0;
   fY2=0;
   for (int i=0; i<size; i++)
      fFixedParams[i]=0;
   fIsSet=kFALSE;
   fChisquare=0;

}


////////////////////////////////////////////////////////////////////////////////
///set the number of dimensions

void TLinearFitter::SetDim(Int_t ndim)
{
   fNdim=ndim;
   fY.ResizeTo(ndim+1);
   fX.ResizeTo(ndim+1, ndim);
   fE.ResizeTo(ndim+1);

   fNpoints=0;
   fIsSet=kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///Additive parts should be separated by "++".
///Examples (ai are parameters to fit):
///1.fitting function: a0*x0 + a1*x1 + a2*x2
///  input formula "x[0]++x[1]++x[2]"
///2.TMath functions can be used:
///  fitting function: a0*TMath::Gaus(x, 0, 1) + a1*y
///  input formula:    "TMath::Gaus(x, 0, 1)++y"
///fills the array of functions

void TLinearFitter::SetFormula(const char *formula)
{
   Int_t size = 0, special = 0;
   Int_t i;
   //Int_t len = strlen(formula);
   if (fInputFunction)
      fInputFunction = nullptr;
   if (fFormula)
      delete [] fFormula;
   fFormulaSize = strlen(formula);
   fFormula = new char[fFormulaSize+1];
   strlcpy(fFormula, formula,fFormulaSize+1);
   fSpecial = 0;
   //in case of a hyperplane:
   char *fstring;
   fstring = (char *)strstr(fFormula, "hyp");
   if (fstring!=0){
      // isHyper = kTRUE;
      fstring+=3;
      sscanf(fstring, "%d", &size);
      //+1 for the constant term
      size++;
      fSpecial=200+size;
   }

   if (fSpecial==0) {
      //in case it's not a hyperplane
      TString sstring(fFormula);
      sstring = sstring.ReplaceAll("++", 2, "|", 1);
      TString replaceformula;

      //count the number of functions

      TObjArray *oa = sstring.Tokenize("|");

      //change the size of functions array and clear it
      if (!fFunctions.IsEmpty())
         fFunctions.Clear();

      // do not own the functions in this case
      fFunctions.SetOwner(kFALSE);

      fNfunctions = oa->GetEntriesFast();
      fFunctions.Expand(fNfunctions);

      //check if the old notation of xi is used somewhere instead of x[i]
      char pattern[12];
      char replacement[14];
      for (i=0; i<fNdim; i++){
         snprintf(pattern,sizeof(pattern), "x%d", i);
         snprintf(replacement,sizeof(replacement), "x[%d]", i);
         sstring = sstring.ReplaceAll(pattern, Int_t(i/10)+2, replacement, Int_t(i/10)+4);
      }

      //fill the array of functions
      oa = sstring.Tokenize("|");
      for (i=0; i<fNfunctions; i++) {
         replaceformula = ((TObjString *)oa->UncheckedAt(i))->GetString();
         // look first if exists in the map
         TFormula * f = nullptr;
         auto elem = fgFormulaMap.find(replaceformula );
         if (elem != fgFormulaMap.end() )
            f = elem->second;
         else {
            // create a new formula and add in the static map
            f = new TFormula("f", replaceformula.Data());
            {
               R__LOCKGUARD(gROOTMutex);
               fgFormulaMap[replaceformula]=f;
            }
         }
         if (!f) {
            Error("TLinearFitter", "f_linear not allocated");
            return;
         }
         special=f->GetNumber();
         fFunctions.Add(f);
      }

      if ((fNfunctions==1)&&(special>299)&&(special<310)){
         //if fitting a polynomial
         size=special-299;
         fSpecial=100+size;
      } else
         size=fNfunctions;
      oa->Delete();
      delete oa;
   }
   fNfunctions=size;
   //change the size of design matrix
   fDesign.ResizeTo(size, size);
   fAtb.ResizeTo(size);
   fDesignTemp.ResizeTo(size, size);
   fDesignTemp2.ResizeTo(size, size);
   fDesignTemp3.ResizeTo(size, size);
   fAtbTemp.ResizeTo(size);
   fAtbTemp2.ResizeTo(size);
   fAtbTemp3.ResizeTo(size);
   if (fFixedParams)
      delete [] fFixedParams;
   fFixedParams=new Bool_t[size];
   fDesign.Zero();
   fAtb.Zero();
   fDesignTemp.Zero();
   fDesignTemp2.Zero();
   fDesignTemp3.Zero();
   fAtbTemp.Zero();
   fAtbTemp2.Zero();
   fAtbTemp3.Zero();
   fY2Temp=0;
   fY2=0;
   for (i=0; i<size; i++)
      fFixedParams[i]=0;
   fIsSet=kFALSE;
   fChisquare=0;

}

////////////////////////////////////////////////////////////////////////////////
///Set the fitting function.

void TLinearFitter::SetFormula(TFormula *function)
{
   Int_t special, size;
   fInputFunction=function;
   fNfunctions=fInputFunction->GetNpar();
   fSpecial = 0;
   special=fInputFunction->GetNumber();
   if (!fFunctions.IsEmpty())
      fFunctions.Delete();

   if ((special>299)&&(special<310)){
      //if fitting a polynomial
      size=special-299;
      fSpecial=100+size;
   } else
      size=fNfunctions;

   fNfunctions=size;
   //change the size of design matrix
   fDesign.ResizeTo(size, size);
   fAtb.ResizeTo(size);
   fDesignTemp.ResizeTo(size, size);
   fAtbTemp.ResizeTo(size);

   fDesignTemp2.ResizeTo(size, size);
   fDesignTemp3.ResizeTo(size, size);

   fAtbTemp2.ResizeTo(size);
   fAtbTemp3.ResizeTo(size);
   //
   if (fFixedParams)
      delete [] fFixedParams;
   fFixedParams=new Bool_t[size];
   fDesign.Zero();
   fAtb.Zero();
   fDesignTemp.Zero();
   fAtbTemp.Zero();

   fDesignTemp2.Zero();
   fDesignTemp3.Zero();

   fAtbTemp2.Zero();
   fAtbTemp3.Zero();
   fY2Temp=0;
   fY2=0;
   for (Int_t i=0; i<size; i++)
      fFixedParams[i]=0;
   //check if any parameters are fixed (not for pure TFormula)

   if (function->InheritsFrom(TF1::Class())){
      Double_t al,bl;
      for (Int_t i=0;i<fNfunctions;i++) {
         ((TF1*)function)->GetParLimits(i,al,bl);
         if (al*bl !=0 && al >= bl) {
            FixParameter(i, function->GetParameter(i));
         }
      }
   }

   fIsSet=kFALSE;
   fChisquare=0;

}

////////////////////////////////////////////////////////////////////////////////
///Update the design matrix after the formula has been changed.

Bool_t TLinearFitter::UpdateMatrix()
{
   if (fStoreData) {
      for (Int_t i=0; i<fNpoints; i++) {
         AddToDesign(TMatrixDRow(fX, i).GetPtr(), fY(i), fE(i));
      }
      return 1;
   } else
      return 0;

}

////////////////////////////////////////////////////////////////////////////////
///To use in TGraph::Fit and TH1::Fit().

Int_t TLinearFitter::ExecuteCommand(const char *command, Double_t *args, Int_t /*nargs*/)
{
   if (!strcmp(command, "FitGraph")){
      if (args)      return GraphLinearFitter(args[0]);
      else           return GraphLinearFitter(0);
   }
   if (!strcmp(command, "FitGraph2D")){
      if (args)      return Graph2DLinearFitter(args[0]);
      else           return Graph2DLinearFitter(0);
   }
   if (!strcmp(command, "FitMultiGraph")){
      if (args)     return  MultiGraphLinearFitter(args[0]);
      else          return  MultiGraphLinearFitter(0);
   }
   if (!strcmp(command, "FitHist"))       return HistLinearFitter();
//   if (!strcmp(command, "FitMultiGraph")) MultiGraphLinearFitter();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Level = 3 (to be consistent with minuit)  prints parameters and parameter
/// errors.

void TLinearFitter::PrintResults(Int_t level, Double_t /*amin*/) const
{
   if (level==3){
      if (!fRobust){
         printf("Fitting results:\nParameters:\nNO.\t\tVALUE\t\tERROR\n");
         for (Int_t i=0; i<fNfunctions; i++){
            printf("%d\t%e\t%e\n", i, fParams(i), TMath::Sqrt(fParCovar(i, i)));
         }
      } else {
         printf("Fitting results:\nParameters:\nNO.\t\tVALUE\n");
         for (Int_t i=0; i<fNfunctions; i++){
            printf("%d\t%e\n", i, fParams(i));
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///Used in TGraph::Fit().

Int_t TLinearFitter::GraphLinearFitter(Double_t h)
{
   StoreData(kFALSE);
   TGraph *grr=(TGraph*)GetObjectFit();
   TF1 *f1=(TF1*)GetUserFunc();
   Foption_t fitOption=GetFitOption();

   //Int_t np=0;
   Double_t *x=grr->GetX();
   Double_t *y=grr->GetY();
   Double_t e;

   Int_t fitResult = 0;
   //set the fitting formula
   SetDim(1);
   SetFormula(f1->GetFormula());

   if (fitOption.Robust){
      fRobust=kTRUE;
      StoreData(kTRUE);
   }
   //put the points into the fitter
   Int_t n=grr->GetN();
   for (Int_t i=0; i<n; i++){
      if (!f1->IsInside(&x[i])) continue;
      e=grr->GetErrorY(i);
      if (e<0 || fitOption.W1)
         e=1;
      AddPoint(&x[i], y[i], e);
   }

   if (fitOption.Robust){
      return EvalRobust(h);
   }

   fitResult = Eval();

   //calculate the precise chisquare
   if (!fitResult){
      if (!fitOption.Nochisq){
         Double_t temp, temp2, sumtotal=0;
         for (Int_t i=0; i<n; i++){
            if (!f1->IsInside(&x[i])) continue;
            temp=f1->Eval(x[i]);
            temp2=(y[i]-temp)*(y[i]-temp);
            e=grr->GetErrorY(i);
            if (e<0 || fitOption.W1)
               e=1;
            temp2/=(e*e);

            sumtotal+=temp2;
         }
         fChisquare=sumtotal;
         f1->SetChisquare(fChisquare);
      }
   }
   return fitResult;
}

////////////////////////////////////////////////////////////////////////////////
///Minimisation function for a TGraph2D

Int_t TLinearFitter::Graph2DLinearFitter(Double_t h)
{
   StoreData(kFALSE);

   TGraph2D *gr=(TGraph2D*)GetObjectFit();
   TF2 *f2=(TF2*)GetUserFunc();

   Foption_t fitOption=GetFitOption();
   Int_t n        = gr->GetN();
   Double_t *gx   = gr->GetX();
   Double_t *gy   = gr->GetY();
   Double_t *gz   = gr->GetZ();
   Double_t x[2];
   Double_t z, e;
   Int_t fitResult=0;
   SetDim(2);
   SetFormula(f2->GetFormula());

   if (fitOption.Robust){
      fRobust=kTRUE;
      StoreData(kTRUE);
   }

   for (Int_t bin=0;bin<n;bin++) {
      x[0] = gx[bin];
      x[1] = gy[bin];
      if (!f2->IsInside(x)) {
         continue;
      }
      z   = gz[bin];
      e=gr->GetErrorZ(bin);
      if (e<0 || fitOption.W1)
         e=1;
      AddPoint(x, z, e);
   }

   if (fitOption.Robust){
      return EvalRobust(h);
   }

   fitResult = Eval();

   if (!fitResult){
      if (!fitOption.Nochisq){
         //calculate the precise chisquare
         Double_t temp, temp2, sumtotal=0;
         for (Int_t bin=0; bin<n; bin++){
            x[0] = gx[bin];
            x[1] = gy[bin];
            if (!f2->IsInside(x)) {
               continue;
            }
            z   = gz[bin];

            temp=f2->Eval(x[0], x[1]);
            temp2=(z-temp)*(z-temp);
            e=gr->GetErrorZ(bin);
            if (e<0 || fitOption.W1)
               e=1;
            temp2/=(e*e);

            sumtotal+=temp2;
         }
         fChisquare=sumtotal;
         f2->SetChisquare(fChisquare);
      }
   }
   return fitResult;
}

////////////////////////////////////////////////////////////////////////////////
///Minimisation function for a TMultiGraph

Int_t TLinearFitter::MultiGraphLinearFitter(Double_t h)
{
   Int_t n, i;
   Double_t *gx, *gy;
   Double_t e;
   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TMultiGraph *mg     = (TMultiGraph*)grFitter->GetObjectFit();
   TF1 *f1   = (TF1*)grFitter->GetUserFunc();
   Foption_t fitOption = grFitter->GetFitOption();
   Int_t fitResult=0;
   SetDim(1);

   if (fitOption.Robust){
      fRobust=kTRUE;
      StoreData(kTRUE);
   }
   SetFormula(f1->GetFormula());

   TGraph *gr;
   TIter next(mg->GetListOfGraphs());
   while ((gr = (TGraph*) next())) {
      n        = gr->GetN();
      gx   = gr->GetX();
      gy   = gr->GetY();
      for (i=0; i<n; i++){
         if (!f1->IsInside(&gx[i])) continue;
         e=gr->GetErrorY(i);
         if (e<0 || fitOption.W1)
            e=1;
         AddPoint(&gx[i], gy[i], e);
      }
   }

   if (fitOption.Robust){
      return EvalRobust(h);
   }

   fitResult = Eval();

   //calculate the chisquare
   if (!fitResult){
      if (!fitOption.Nochisq){
         Double_t temp, temp2, sumtotal=0;
         next.Reset();
         while((gr = (TGraph*)next())) {
            n        = gr->GetN();
            gx   = gr->GetX();
            gy   = gr->GetY();
            for (i=0; i<n; i++){
               if (!f1->IsInside(&gx[i])) continue;
               temp=f1->Eval(gx[i]);
               temp2=(gy[i]-temp)*(gy[i]-temp);
               e=gr->GetErrorY(i);
               if (e<0 || fitOption.W1)
                  e=1;
               temp2/=(e*e);

               sumtotal+=temp2;
            }

         }
         fChisquare=sumtotal;
         f1->SetChisquare(fChisquare);
      }
   }
   return fitResult;
}

////////////////////////////////////////////////////////////////////////////////
/// Minimization function for H1s using a Chisquare method.

Int_t TLinearFitter::HistLinearFitter()
{
   StoreData(kFALSE);
   Double_t cu,eu;
   // Double_t dersum[100], grad[100];
   Double_t x[3];
   Int_t bin,binx,biny,binz;
   //   Axis_t binlow, binup, binsize;

   TH1 *hfit = (TH1*)GetObjectFit();
   TF1 *f1   = (TF1*)GetUserFunc();
   Int_t fitResult;
   Foption_t fitOption = GetFitOption();
   //SetDim(hfit->GetDimension());
   SetDim(f1->GetNdim());
   SetFormula(f1->GetFormula());

   Int_t hxfirst = GetXfirst();
   Int_t hxlast  = GetXlast();
   Int_t hyfirst = GetYfirst();
   Int_t hylast  = GetYlast();
   Int_t hzfirst = GetZfirst();
   Int_t hzlast  = GetZlast();
   TAxis *xaxis  = hfit->GetXaxis();
   TAxis *yaxis  = hfit->GetYaxis();
   TAxis *zaxis  = hfit->GetZaxis();

   for (binz=hzfirst;binz<=hzlast;binz++) {
      x[2]  = zaxis->GetBinCenter(binz);
      for (biny=hyfirst;biny<=hylast;biny++) {
         x[1]  = yaxis->GetBinCenter(biny);
         for (binx=hxfirst;binx<=hxlast;binx++) {
            x[0]  = xaxis->GetBinCenter(binx);
            if (!f1->IsInside(x)) continue;
            bin = hfit->GetBin(binx,biny,binz);
            cu  = hfit->GetBinContent(bin);
            if (f1->GetNdim() < hfit->GetDimension()) {
               if (hfit->GetDimension() > 2) cu = x[2];
               else                          cu = x[1];
            }
            if (fitOption.W1) {
               if (fitOption.W1==1 && cu == 0) continue;
               eu = 1;
            } else {
               eu  = hfit->GetBinError(bin);
               if (eu <= 0) continue;
            }
            AddPoint(x, cu, eu);
         }
      }
   }

   fitResult = Eval();

   if (!fitResult){
      if (!fitOption.Nochisq){
         Double_t temp, temp2, sumtotal=0;
         for (binz=hzfirst;binz<=hzlast;binz++) {
            x[2]  = zaxis->GetBinCenter(binz);
            for (biny=hyfirst;biny<=hylast;biny++) {
               x[1]  = yaxis->GetBinCenter(biny);
               for (binx=hxfirst;binx<=hxlast;binx++) {
                  x[0]  = xaxis->GetBinCenter(binx);
                  if (!f1->IsInside(x)) continue;
                  bin = hfit->GetBin(binx,biny,binz);
                  cu  = hfit->GetBinContent(bin);

                  if (fitOption.W1) {
                     if (fitOption.W1==1 && cu == 0) continue;
                     eu = 1;
                  } else {
                     eu  = hfit->GetBinError(bin);
                     if (eu <= 0) continue;
                  }
                  temp=f1->EvalPar(x);
                  temp2=(cu-temp)*(cu-temp);
                  temp2/=(eu*eu);
                  sumtotal+=temp2;
               }
            }
         }

         fChisquare=sumtotal;
         f1->SetChisquare(fChisquare);
      }
   }
   return fitResult;
}

////////////////////////////////////////////////////////////////////////////////

void TLinearFitter::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      Int_t old_matr_size = fNfunctions;
      R__b.ReadClassBuffer(TLinearFitter::Class(),this);
      if (old_matr_size < fNfunctions) {
         fDesignTemp.ResizeTo(fNfunctions, fNfunctions);
         fAtbTemp.ResizeTo(fNfunctions);

         fDesignTemp2.ResizeTo(fNfunctions, fNfunctions);
         fDesignTemp3.ResizeTo(fNfunctions, fNfunctions);

         fAtbTemp2.ResizeTo(fNfunctions);
         fAtbTemp3.ResizeTo(fNfunctions);
      }
   } else {
      if (fAtb.NonZeros()==0) AddTempMatrices();
      R__b.WriteClassBuffer(TLinearFitter::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Finds the parameters of the fitted function in case data contains
///outliers.
///Parameter h stands for the minimal fraction of good points in the
///dataset (h < 1, i.e. for 70% of good points take h=0.7).
///The default value of h*Npoints is  (Npoints + Nparameters+1)/2
///If the user provides a value of h smaller than above, default is taken
///See class description for the algorithm details

Int_t TLinearFitter::EvalRobust(Double_t h)
{
   fRobust = kTRUE;
   Double_t kEps = 1e-13;
   Int_t nmini = 300;
   Int_t i, j, maxind=0, k, k1 = 500;
   Int_t nbest = 10;
   Double_t chi2 = -1;

   if (fFunctions.IsEmpty()&&(!fInputFunction)&&(fSpecial<=200)){
      Error("TLinearFitter::EvalRobust", "The formula hasn't been set");
      return 1;
   }


   Double_t *bestchi2 = new Double_t[nbest];
   for (i=0; i<nbest; i++)
      bestchi2[i]=1e30;

   Int_t hdef=Int_t((fNpoints+fNfunctions+1)/2);

   if (h>0.000001 && h<1 && fNpoints*h > hdef)
      fH = Int_t(fNpoints*h);
   else {
      // print message only when h is not zero
      if (h>0) Warning("Fitting:", "illegal value of H, default is taken, h = %3.2f",double(hdef)/fNpoints);
      fH=hdef;
   }

   fDesign.ResizeTo(fNfunctions, fNfunctions);
   fAtb.ResizeTo(fNfunctions);
   fParams.ResizeTo(fNfunctions);

   Int_t *index = new Int_t[fNpoints];
   Double_t *residuals = new Double_t[fNpoints];

   if (fNpoints < 2*nmini) {
      //when number of cases is small

      //to store the best coefficients (columnwise)
      TMatrixD cstock(fNfunctions, nbest);
      for (k = 0; k < k1; k++) {
         CreateSubset(fNpoints, fH, index);
         chi2 = CStep(1, fH, residuals,index, index, -1, -1);
         chi2 = CStep(2, fH, residuals,index, index, -1, -1);
         maxind = TMath::LocMax(nbest, bestchi2);
         if (chi2 < bestchi2[maxind]) {
            bestchi2[maxind] = chi2;
            for (i=0; i<fNfunctions; i++)
               cstock(i, maxind) = fParams(i);
         }
      }

      //for the nbest best results, perform CSteps until convergence
      Int_t *bestindex = new Int_t[fH];
      Double_t currentbest;
      for (i=0; i<nbest; i++) {
         for (j=0; j<fNfunctions; j++)
            fParams(j) = cstock(j, i);
         chi2 = 1;
         while (chi2 > kEps) {
            chi2 = CStep(2, fH, residuals,index, index, -1, -1);
            if (TMath::Abs(chi2 - bestchi2[i]) < kEps)
               break;
            else
               bestchi2[i] = chi2;
         }
         currentbest = TMath::MinElement(nbest, bestchi2);
         if (chi2 <= currentbest + kEps) {
            for (j=0; j<fH; j++){
               bestindex[j]=index[j];
            }
            maxind = i;
         }
         for (j=0; j<fNfunctions; j++)
            cstock(j, i) = fParams(j);
      }
      //report the result with the lowest chisquare
      for (j=0; j<fNfunctions; j++)
         fParams(j) = cstock(j, maxind);
      fFitsample.SetBitNumber(fNpoints, kFALSE);
      for (j=0; j<fH; j++){
         //printf("bestindex[%d]=%d\n", j, bestindex[j]);
         fFitsample.SetBitNumber(bestindex[j]);
      }
      if (fInputFunction && fInputFunction->InheritsFrom(TF1::Class())) {
         ((TF1*)fInputFunction)->SetChisquare(bestchi2[maxind]);
         ((TF1*)fInputFunction)->SetNumberFitPoints(fH);
         ((TF1*)fInputFunction)->SetNDF(fH-fNfunctions);
      }
      delete [] index;
      delete [] bestindex;
      delete [] residuals;
      delete [] bestchi2;
      return 0;
   }
   //if n is large, the dataset should be partitioned
   Int_t indsubdat[5];
   for (i=0; i<5; i++)
      indsubdat[i] = 0;

   Int_t nsub = Partition(nmini, indsubdat);
   Int_t hsub;

   Int_t sum = TMath::Min(nmini*5, fNpoints);

   Int_t *subdat = new Int_t[sum]; //to store the indices of selected cases
   RDraw(subdat, indsubdat);

   TMatrixD cstockbig(fNfunctions, nbest*5);
   Int_t *beststock = new Int_t[nbest];
   Int_t i_start = 0;
   Int_t i_end = indsubdat[0];
   Int_t k2 = Int_t(k1/nsub);
   for (Int_t kgroup = 0; kgroup < nsub; kgroup++) {

      hsub = Int_t(fH * indsubdat[kgroup]/fNpoints);
      for (i=0; i<nbest; i++)
         bestchi2[i] = 1e16;
      for (k=0; k<k2; k++) {
         CreateSubset(indsubdat[kgroup], hsub, index);
         chi2 = CStep(1, hsub, residuals, index, subdat, i_start, i_end);
         chi2 = CStep(2, hsub, residuals, index, subdat, i_start, i_end);
         maxind = TMath::LocMax(nbest, bestchi2);
         if (chi2 < bestchi2[maxind]){
            for (i=0; i<fNfunctions; i++)
               cstockbig(i, nbest*kgroup + maxind) = fParams(i);
            bestchi2[maxind] = chi2;
         }
      }
      if (kgroup != nsub - 1){
         i_start += indsubdat[kgroup];
         i_end += indsubdat[kgroup+1];
      }
   }

   for (i=0; i<nbest; i++)
      bestchi2[i] = 1e30;
   //on the pooled subset
   Int_t hsub2 = Int_t(fH*sum/fNpoints);
   for (k=0; k<nbest*5; k++) {
      for (i=0; i<fNfunctions; i++)
         fParams(i)=cstockbig(i, k);
      chi2 = CStep(1, hsub2, residuals, index, subdat, 0, sum);
      chi2 = CStep(2, hsub2, residuals, index, subdat, 0, sum);
      maxind = TMath::LocMax(nbest, bestchi2);
      if (chi2 < bestchi2[maxind]){
         beststock[maxind] = k;
         bestchi2[maxind] = chi2;
      }
   }

   //now the array beststock keeps indices of 10 best candidates in cstockbig matrix
   for (k=0; k<nbest; k++) {
      for (i=0; i<fNfunctions; i++)
         fParams(i) = cstockbig(i, beststock[k]);
      chi2 = CStep(1, fH, residuals, index, index, -1, -1);
      chi2 = CStep(2, fH, residuals, index, index, -1, -1);
      bestchi2[k] = chi2;
   }

   maxind = TMath::LocMin(nbest, bestchi2);
   for (i=0; i<fNfunctions; i++)
      fParams(i)=cstockbig(i, beststock[maxind]);

   chi2 = 1;
   while (chi2 > kEps) {
      chi2 = CStep(2, fH, residuals, index, index, -1, -1);
      if (TMath::Abs(chi2 - bestchi2[maxind]) < kEps)
         break;
      else
         bestchi2[maxind] = chi2;
   }

   fFitsample.SetBitNumber(fNpoints, kFALSE);
   for (j=0; j<fH; j++)
      fFitsample.SetBitNumber(index[j]);
   if (fInputFunction){
      ((TF1*)fInputFunction)->SetChisquare(bestchi2[maxind]);
      ((TF1*)fInputFunction)->SetNumberFitPoints(fH);
      ((TF1*)fInputFunction)->SetNDF(fH-fNfunctions);
   }

   delete [] subdat;
   delete [] beststock;
   delete [] bestchi2;
   delete [] residuals;
   delete [] index;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///Creates a p-subset to start
///ntotal - total number of points from which the subset is chosen

void TLinearFitter::CreateSubset(Int_t ntotal, Int_t h, Int_t *index)
{
   Int_t i, j;
   Bool_t repeat=kFALSE;
   Int_t nindex=0;
   Int_t num;
   for(i=0; i<ntotal; i++)
      index[i] = ntotal+1;

   TRandom2 r;
   //create a p-subset
   for (i=0; i<fNfunctions; i++) {
      num=Int_t(r.Uniform(0, 1)*(ntotal-1));
      if (i>0){
         for(j=0; j<=i-1; j++) {
            if(index[j]==num)
               repeat = kTRUE;
         }
      }
      if(repeat==kTRUE) {
         i--;
         repeat = kFALSE;
      } else {
         index[i] = num;
         nindex++;
      }
   }

   //compute the coefficients of a hyperplane through the p-subset
   fDesign.Zero();
   fAtb.Zero();
   for (i=0; i<fNfunctions; i++){
      AddToDesign(TMatrixDRow(fX, index[i]).GetPtr(), fY(index[i]), fE(index[i]));
   }
   Bool_t ok;

   ok = Linf();

   //if the chosen points don't define a hyperplane, add more
   while (!ok && (nindex < h)) {
      repeat=kFALSE;
      do{
         num=Int_t(r.Uniform(0,1)*(ntotal-1));
         repeat=kFALSE;
         for(i=0; i<nindex; i++) {
            if(index[i]==num) {
               repeat=kTRUE;
               break;
            }
         }
      } while(repeat==kTRUE);

      index[nindex] = num;
      nindex++;
      //check if the system is of full rank now
      AddToDesign(TMatrixDRow(fX, index[nindex-1]).GetPtr(), fY(index[nindex-1]), fE(index[nindex-1]));
      ok = Linf();
   }
}

////////////////////////////////////////////////////////////////////////////////
///The CStep procedure, as described in the article

Double_t TLinearFitter::CStep(Int_t step, Int_t h, Double_t *residuals, Int_t *index, Int_t *subdat, Int_t start, Int_t end)
{
   R__ASSERT( !fFunctions.IsEmpty() || fInputFunction ||  fSpecial>200);

   Int_t i, j, itemp, n;
   Double_t func;
   Double_t val[100];
   Int_t npar;
   if (start > -1) {
      n = end - start;
      for (i=0; i<n; i++) {
         func = 0;
         itemp = subdat[start+i];
         if (fInputFunction){
            fInputFunction->SetParameters(fParams.GetMatrixArray());
            func=fInputFunction->EvalPar(TMatrixDRow(fX, itemp).GetPtr());
         } else {
            func=0;
            if ((fSpecial>100)&&(fSpecial<200)){
                  npar = fSpecial-100;
                  val[0] = 1;
                  for (j=1; j<npar; j++)
                     val[j] = val[j-1]*fX(itemp, 0);
                  for (j=0; j<npar; j++)
                     func += fParams(j)*val[j];
            } else if (fSpecial>200) {
               //hyperplane case
               npar = fSpecial-201;
               func+=fParams(0);
               for (j=0; j<npar; j++)
                  func += fParams(j+1)*fX(itemp, j);
            } else {
               // general case
               for (j=0; j<fNfunctions; j++) {
                  TF1 *f1 = (TF1*)(fFunctions.UncheckedAt(j));
                  val[j] = f1->EvalPar(TMatrixDRow(fX, itemp).GetPtr());
                  func += fParams(j)*val[j];
               }
            }
         }
         residuals[i] = (fY(itemp) - func)*(fY(itemp) - func)/(fE(i)*fE(i));
      }
   } else {
      n=fNpoints;
      for (i=0; i<fNpoints; i++) {
         func = 0;
         if (fInputFunction){
            fInputFunction->SetParameters(fParams.GetMatrixArray());
            func=fInputFunction->EvalPar(TMatrixDRow(fX, i).GetPtr());
         } else {
            func=0;
            if ((fSpecial>100)&&(fSpecial<200)){
               npar = fSpecial-100;
               val[0] = 1;
               for (j=1; j<npar; j++)
                  val[j] = val[j-1]*fX(i, 0);
               for (j=0; j<npar; j++)
                  func += fParams(j)*val[j];
            } else if (fSpecial>200) {
               //hyperplane case
               npar = fSpecial-201;
               func+=fParams(0);
               for (j=0; j<npar; j++)
                  func += fParams(j+1)*fX(i, j);
            } else {
               for (j=0; j<fNfunctions; j++) {
                  TF1 *f1 = (TF1*)(fFunctions.UncheckedAt(j));
                  val[j] = f1->EvalPar(TMatrixDRow(fX, i).GetPtr());
                  func += fParams(j)*val[j];
               }
            }
         }
         residuals[i] = (fY(i) - func)*(fY(i) - func)/(fE(i)*fE(i));
      }
   }
   //take h with smallest residuals
   TMath::KOrdStat(n, residuals, h-1, index);
   //add them to the design matrix
   fDesign.Zero();
   fAtb.Zero();
   for (i=0; i<h; i++)
      AddToDesign(TMatrixDRow(fX, index[i]).GetPtr(), fY(index[i]), fE(index[i]));

   Linf();

   //don't calculate the chisquare at the 1st cstep
   if (step==1) return 0;
   Double_t sum=0;


   if (start > -1) {
      for (i=0; i<h; i++) {
         itemp = subdat[start+index[i]];
         if (fInputFunction){
            fInputFunction->SetParameters(fParams.GetMatrixArray());
            func=fInputFunction->EvalPar(TMatrixDRow(fX, itemp).GetPtr());
         } else {
            func=0;
            if ((fSpecial>100)&&(fSpecial<200)){
                  npar = fSpecial-100;
                  val[0] = 1;
                  for (j=1; j<npar; j++)
                     val[j] = val[j-1]*fX(itemp, 0);
                  for (j=0; j<npar; j++)
                     func += fParams(j)*val[j];
            } else if (fSpecial>200) {
               //hyperplane case
               npar = fSpecial-201;
               func+=fParams(0);
               for (j=0; j<npar; j++)
                  func += fParams(j+1)*fX(itemp, j);
            } else {
               for (j=0; j<fNfunctions; j++) {
                  TF1 *f1 = (TF1*)(fFunctions.UncheckedAt(j));
                  val[j] = f1->EvalPar(TMatrixDRow(fX, itemp).GetPtr());
                  func += fParams(j)*val[j];
               }
            }
         }
         sum+=(fY(itemp)-func)*(fY(itemp)-func)/(fE(itemp)*fE(itemp));
      }
   } else {
      for (i=0; i<h; i++) {
         if (fInputFunction){
            fInputFunction->SetParameters(fParams.GetMatrixArray());
            func=fInputFunction->EvalPar(TMatrixDRow(fX, index[i]).GetPtr());
         } else {
            func=0;
            if ((fSpecial>100)&&(fSpecial<200)){
               npar = fSpecial-100;
               val[0] = 1;
               for (j=1; j<npar; j++)
                  val[j] = val[j-1]*fX(index[i], 0);
               for (j=0; j<npar; j++)
                  func += fParams(j)*val[j];
            } else {
               if (fSpecial>200) {
                  //hyperplane case
                  npar = fSpecial-201;
                  func+=fParams(0);
                  for (j=0; j<npar; j++)
                     func += fParams(j+1)*fX(index[i], j);
               } else {
                  for (j=0; j<fNfunctions; j++) {
                     TF1 *f1 = (TF1*)(fFunctions.UncheckedAt(j));
                     val[j] = f1->EvalPar(TMatrixDRow(fX, index[i]).GetPtr());
                     func += fParams(j)*val[j];
                  }
               }
            }
         }

         sum+=(fY(index[i])-func)*(fY(index[i])-func)/(fE(index[i])*fE(index[i]));
      }
   }

   return sum;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TLinearFitter::Linf()
{
   //currently without the intercept term
   fDesignTemp2+=fDesignTemp3;
   fDesignTemp+=fDesignTemp2;
   fDesign+=fDesignTemp;
   fDesignTemp3.Zero();
   fDesignTemp2.Zero();
   fDesignTemp.Zero();
   fAtbTemp2+=fAtbTemp3;
   fAtbTemp+=fAtbTemp2;
   fAtb+=fAtbTemp;
   fAtbTemp3.Zero();
   fAtbTemp2.Zero();
   fAtbTemp.Zero();

   fY2+=fY2Temp;
   fY2Temp=0;


   TDecompChol chol(fDesign);
   TVectorD temp(fNfunctions);
   Bool_t ok;
   temp = chol.Solve(fAtb, ok);
   if (!ok){
      Error("Linf","Matrix inversion failed");
      //fDesign.Print();
      fParams.Zero();
      return kFALSE;
   }
   fParams = temp;
   return ok;
}

////////////////////////////////////////////////////////////////////////////////
///divides the elements into approximately equal subgroups
///number of elements in each subgroup is stored in indsubdat
///number of subgroups is returned

Int_t TLinearFitter::Partition(Int_t nmini, Int_t *indsubdat)
{
   Int_t nsub;

   if ((fNpoints>=2*nmini) && (fNpoints<=(3*nmini-1))) {
      if (fNpoints%2==1){
         indsubdat[0]=Int_t(fNpoints*0.5);
         indsubdat[1]=Int_t(fNpoints*0.5)+1;
      } else
         indsubdat[0]=indsubdat[1]=Int_t(fNpoints/2);
      nsub=2;
   }
   else{
      if((fNpoints>=3*nmini) && (fNpoints<(4*nmini -1))) {
         if(fNpoints%3==0){
            indsubdat[0]=indsubdat[1]=indsubdat[2]=Int_t(fNpoints/3);
         } else {
            indsubdat[0]=Int_t(fNpoints/3);
            indsubdat[1]=Int_t(fNpoints/3)+1;
            if (fNpoints%3==1) indsubdat[2]=Int_t(fNpoints/3);
            else indsubdat[2]=Int_t(fNpoints/3)+1;
         }
         nsub=3;
      }
      else{
         if((fNpoints>=4*nmini)&&(fNpoints<=(5*nmini-1))){
            if (fNpoints%4==0) indsubdat[0]=indsubdat[1]=indsubdat[2]=indsubdat[3]=Int_t(fNpoints/4);
            else {
               indsubdat[0]=Int_t(fNpoints/4);
               indsubdat[1]=Int_t(fNpoints/4)+1;
               if(fNpoints%4==1) indsubdat[2]=indsubdat[3]=Int_t(fNpoints/4);
               if(fNpoints%4==2) {
                  indsubdat[2]=Int_t(fNpoints/4)+1;
                  indsubdat[3]=Int_t(fNpoints/4);
               }
               if(fNpoints%4==3) indsubdat[2]=indsubdat[3]=Int_t(fNpoints/4)+1;
            }
            nsub=4;
         } else {
            for(Int_t i=0; i<5; i++)
               indsubdat[i]=nmini;
            nsub=5;
         }
      }
   }
   return nsub;
}

////////////////////////////////////////////////////////////////////////////////
///Draws ngroup nonoverlapping subdatasets out of a dataset of size n
///such that the selected case numbers are uniformly distributed from 1 to n

void TLinearFitter::RDraw(Int_t *subdat, Int_t *indsubdat)
{
   Int_t jndex = 0;
   Int_t nrand;
   Int_t i, k, m, j;
   Int_t ngroup=0;
   for (i=0; i<5; i++) {
      if (indsubdat[i]!=0)
         ngroup++;
   }
   TRandom r;
   for (k=1; k<=ngroup; k++) {
      for (m=1; m<=indsubdat[k-1]; m++) {
         nrand = Int_t(r.Uniform(0, 1) * (fNpoints-jndex)) + 1;
         jndex++;
         if (jndex==1) {
            subdat[0] = nrand;
         } else {
            subdat[jndex-1] = nrand + jndex - 2;
            for (i=1; i<=jndex-1; i++) {
               if(subdat[i-1] > nrand+i-2) {
                  for(j=jndex; j>=i+1; j--) {
                     subdat[j-1] = subdat[j-2];
                  }
                  subdat[i-1] = nrand+i-2;
                  break;  //breaking the loop for(i=1...
               }
            }
         }
      }
   }
}


