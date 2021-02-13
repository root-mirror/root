// @(#)root/hist:$Id$
// Author: Rene Brun   18/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "strlcpy.h"
#include "snprintf.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TObjString.h"
#include "TInterpreter.h"
#include "TPluginManager.h"
#include "TBrowser.h"
#include "TColor.h"
#include "TMethodCall.h"
#include "TF1Helper.h"
#include "TF1NormSum.h"
#include "TF1Convolution.h"
#include "TVirtualMutex.h"
#include "Math/WrappedFunction.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/BrentMethods.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IntegratorOptions.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/RichardsonDerivator.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "Math/Factory.h"
#include "Math/ChebyshevPol.h"
#include "Fit/FitResult.h"
// for I/O backward compatibility
#include "v5/TF1Data.h"

#include "AnalyticalIntegrals.h"

std::atomic<Bool_t> TF1::fgAbsValue(kFALSE);
Bool_t TF1::fgRejectPoint = kFALSE;
std::atomic<Bool_t> TF1::fgAddToGlobList(kTRUE);
static Double_t gErrorTF1 = 0;

using TF1Updater_t = void (*)(Int_t nobjects, TObject **from, TObject **to);
bool R__SetClonesArrayTF1Updater(TF1Updater_t func);


namespace {
struct TF1v5Convert : public TF1 {
public:
   void Convert(ROOT::v5::TF1Data &from)
   {
      // convert old TF1 to new one
      fNpar = from.GetNpar();
      fNdim = from.GetNdim();
      if (from.fType == 0) {
         // formula functions
         // if ndim is not 1  set xmin max to zero to avoid error in ctor
         double xmin = from.fXmin;
         double xmax = from.fXmax;
         if (fNdim > 1) {
            xmin = 0;
            xmax = 0;
         }
         TF1 fnew(from.GetName(), from.GetExpFormula(), xmin, xmax);
         if (fNdim > 1) {
            fnew.SetRange(from.fXmin, from.fXmax);
         }
         fnew.Copy(*this);
         // need to set parameter values
         if (from.GetParameters())
            fFormula->SetParameters(from.GetParameters());
      } else {
         // case of a function pointers
         fParams.reset(new TF1Parameters(fNpar));
         fName = from.GetName();
         fTitle = from.GetTitle();
         // need to set parameter values
         if (from.GetParameters())
            fParams->SetParameters(from.GetParameters());
      }
      // copy the other data members
      fNpx = from.fNpx;
      fType = (EFType)from.fType;
      fNpfits = from.fNpfits;
      fNDF = from.fNDF;
      fChisquare = from.fChisquare;
      fMaximum = from.fMaximum;
      fMinimum = from.fMinimum;
      fXmin = from.fXmin;
      fXmax = from.fXmax;

      if (from.fParErrors)
         fParErrors = std::vector<Double_t>(from.fParErrors, from.fParErrors + fNpar);
      if (from.fParMin)
         fParMin = std::vector<Double_t>(from.fParMin, from.fParMin + fNpar);
      if (from.fParMax)
         fParMax = std::vector<Double_t>(from.fParMax, from.fParMax + fNpar);
      if (from.fNsave > 0) {
         assert(from.fSave);
         fSave = std::vector<Double_t>(from.fSave, from.fSave + from.fNsave);
      }
      // set the bits
      for (int ibit = 0; ibit < 24; ++ibit)
         if (from.TestBit(BIT(ibit)))
            SetBit(BIT(ibit));

      // copy the graph attributes
      auto &fromLine = static_cast<TAttLine &>(from);
      fromLine.Copy(*this);
      auto &fromFill = static_cast<TAttFill &>(from);
      fromFill.Copy(*this);
      auto &fromMarker = static_cast<TAttMarker &>(from);
      fromMarker.Copy(*this);
   }
};
} // unnamed namespace

static void R__v5TF1Updater(Int_t nobjects, TObject **from, TObject **to)
{
   auto **fromv5 = (ROOT::v5::TF1Data **)from;
   auto **target = (TF1v5Convert **)to;

   for (int i = 0; i < nobjects; ++i) {
      if (fromv5[i] && target[i])
         target[i]->Convert(*fromv5[i]);
   }
}

int R__RegisterTF1UpdaterTrigger = R__SetClonesArrayTF1Updater(R__v5TF1Updater);

ClassImp(TF1);

// class wrapping evaluation of TF1(x) - y0
class GFunc {
   const TF1 *fFunction;
   const double fY0;
public:
   GFunc(const TF1 *function , double y): fFunction(function), fY0(y) {}
   double operator()(double x) const
   {
      return fFunction->Eval(x) - fY0;
   }
};

// class wrapping evaluation of -TF1(x)
class GInverseFunc {
   const TF1 *fFunction;
public:
   GInverseFunc(const TF1 *function): fFunction(function) {}

   double operator()(double x) const
   {
      return - fFunction->Eval(x);
   }
};
// class wrapping evaluation of -TF1(x) for multi-dimension
class GInverseFuncNdim {
   TF1 *fFunction;
public:
   GInverseFuncNdim(TF1 *function): fFunction(function) {}

   double operator()(const double *x) const
   {
      return - fFunction->EvalPar(x, (Double_t *)0);
   }
};

// class wrapping function evaluation directly in 1D interface (used for integration)
// and implementing the methods for the momentum calculations

class  TF1_EvalWrapper : public ROOT::Math::IGenFunction {
public:
   TF1_EvalWrapper(TF1 *f, const Double_t *par, bool useAbsVal, Double_t n = 1, Double_t x0 = 0) :
      fFunc(f),
      fPar(((par) ? par : f->GetParameters())),
      fAbsVal(useAbsVal),
      fN(n),
      fX0(x0)
   {
      fFunc->InitArgs(fX, fPar);
      if (par) fFunc->SetParameters(par);
   }

   ROOT::Math::IGenFunction *Clone()  const
   {
      // use default copy constructor
      TF1_EvalWrapper *f =  new TF1_EvalWrapper(*this);
      f->fFunc->InitArgs(f->fX, f->fPar);
      return f;
   }
   // evaluate |f(x)|
   Double_t DoEval(Double_t x) const
   {
      // use evaluation with stored parameters (i.e. pass zero)
      fX[0] = x;
      Double_t fval = fFunc->EvalPar(fX, 0);
      if (fAbsVal && fval < 0)  return -fval;
      return fval;
   }
   // evaluate x * |f(x)|
   Double_t EvalFirstMom(Double_t x)
   {
      fX[0] = x;
      return fX[0] * TMath::Abs(fFunc->EvalPar(fX, 0));
   }
   // evaluate (x - x0) ^n * f(x)
   Double_t EvalNMom(Double_t x) const
   {
      fX[0] = x;
      return TMath::Power(fX[0] - fX0, fN) * TMath::Abs(fFunc->EvalPar(fX, 0));
   }

   TF1 *fFunc;
   mutable Double_t fX[1];
   const double *fPar;
   Bool_t fAbsVal;
   Double_t fN;
   Double_t fX0;
};

////////////////////////////////////////////////////////////////////////////////
/** \class TF1
    \ingroup Hist
    \brief 1-Dim function class


## TF1: 1-Dim function class

A TF1 object is a 1-Dim function defined between a lower and upper limit.
The function may be a simple function based on a TFormula expression or a precompiled user function.
The function may have associated parameters.
TF1 graphics function is via the TH1 and TGraph drawing functions.

The following types of functions can be created:

1.  [Expression using variable x and no parameters]([#F1)
2.  [Expression using variable x with parameters](#F2)
3.  [Lambda Expression with variable x and parameters](#F3)
4.  [A general C function with parameters](#F4)
5.  [A general C++ function object (functor) with parameters](#F5)
6.  [A member function with parameters of a general C++ class](#F6)



### <a name="F1"></a> 1 - Expression using variable x and no parameters

#### Case 1: inline expression using standard C++ functions/operators

Begin_Macro(source)
{
   TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10);
   fa1->Draw();
}
End_Macro

#### Case 2: inline expression using a ROOT function (e.g. from TMath) without parameters


Begin_Macro(source)
{
   TF1 *fa2 = new TF1("fa2","TMath::DiLog(x)",0,10);
   fa2->Draw();
}
End_Macro

#### Case 3: inline expression using a user defined CLING function by name

~~~~{.cpp}
Double_t myFunc(double x) { return x+sin(x); }
....
TF1 *fa3 = new TF1("fa3","myFunc(x)",-3,5);
fa3->Draw();
~~~~

### <a name="F2"></a> 2 - Expression using variable x with parameters

#### Case 1: inline expression using standard C++ functions/operators

* Example a:


~~~~{.cpp}
TF1 *fa = new TF1("fa","[0]*x*sin([1]*x)",-3,3);
~~~~

This creates a function of variable x with 2 parameters. The parameters must be initialized via:

~~~~{.cpp}
   fa->SetParameter(0,value_first_parameter);
   fa->SetParameter(1,value_second_parameter);
~~~~


Parameters may be given a name:

~~~~{.cpp}
    fa->SetParName(0,"Constant");
~~~~

* Example b:

~~~~{.cpp}
    TF1 *fb = new TF1("fb","gaus(0)*expo(3)",0,10);
~~~~


``gaus(0)`` is a substitute for ``[0]*exp(-0.5*((x-[1])/[2])**2)`` and ``(0)`` means start numbering parameters at ``0``. ``expo(3)`` is a substitute for ``exp([3]+[4]*x)``.

#### Case 2: inline expression using TMath functions with parameters

~~~~{.cpp}
   TF1 *fb2 = new TF1("fa3","TMath::Landau(x,[0],[1],0)",-5,10);
   fb2->SetParameters(0.2,1.3);
   fb2->Draw();
~~~~



Begin_Macro
{
    TCanvas *c = new TCanvas("c","c",0,0,500,300);
    TF1 *fb2 = new TF1("fa3","TMath::Landau(x,[0],[1],0)",-5,10);
    fb2->SetParameters(0.2,1.3); fb2->Draw();
    return c;
}
End_Macro

###<a name="F3"></a> 3 - A lambda expression with variables and parameters

\since **6.00/00:**
TF1 supports using lambda expressions in the formula. This allows, by using a full C++ syntax the full power of lambda
functions and still maintain the capability of storing the function in a file which cannot be done with
function pointer or lambda written not as expression, but as code (see items below).

Example on how using lambda to define a sum of two functions.
Note that is necessary to provide the number of parameters

~~~~{.cpp}
TF1 f1("f1","sin(x)",0,10);
TF1 f2("f2","cos(x)",0,10);
TF1 fsum("f1","[&](double *x, double *p){ return p[0]*f1(x) + p[1]*f2(x); }",0,10,2);
~~~~

###<a name="F4"></a> 4 - A general C function with parameters

Consider the macro myfunc.C below:

~~~~{.cpp}
   // Macro myfunc.C
   Double_t myfunction(Double_t *x, Double_t *par)
   {
      Float_t xx =x[0];
      Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
      return f;
   }
   void myfunc()
   {
      TF1 *f1 = new TF1("myfunc",myfunction,0,10,2);
      f1->SetParameters(2,1);
      f1->SetParNames("constant","coefficient");
      f1->Draw();
   }
   void myfit()
   {
      TH1F *h1=new TH1F("h1","test",100,0,10);
      h1->FillRandom("myfunc",20000);
      TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
      f1->SetParameters(800,1);
      h1->Fit("myfunc");
   }
~~~~



In an interactive session you can do:

~~~~
   Root > .L myfunc.C
   Root > myfunc();
   Root > myfit();
~~~~



TF1 objects can reference other TF1 objects of type A or B defined above. This excludes CLing or compiled functions. However, there is a restriction. A function cannot reference a basic function if the basic function is a polynomial polN.

Example:

~~~~{.cpp}
{
      TF1 *fcos = new TF1 ("fcos", "[0]*cos(x)", 0., 10.);
      fcos->SetParNames( "cos");
      fcos->SetParameter( 0, 1.1);

      TF1 *fsin = new TF1 ("fsin", "[0]*sin(x)", 0., 10.);
      fsin->SetParNames( "sin");
      fsin->SetParameter( 0, 2.1);

      TF1 *fsincos = new TF1 ("fsc", "fcos+fsin");

      TF1 *fs2 = new TF1 ("fs2", "fsc+fsc");
}
~~~~


### <a name="F5"></a> 5 - A general C++ function object (functor) with parameters

A TF1 can be created from any C++ class implementing the operator()(double *x, double *p). The advantage of the function object is that he can have a state and reference therefore what-ever other object. In this way the user can customize his function.

Example:


~~~~{.cpp}
class  MyFunctionObject {
 public:
   // use constructor to customize your function object

   double operator() (double *x, double *p) {
      // function implementation using class data members
   }
};
{
    ....
   MyFunctionObject fobj;
   TF1 * f = new TF1("f",fobj,0,1,npar);    // create TF1 class.
   .....
}
~~~~

#### Using a lambda function as a general C++ functor object

From C++11 we can use both std::function or even better lambda functions to create the TF1.
As above the lambda must have the right signature but can capture whatever we want. For example we can make
a TF1 from the TGraph::Eval function as shown below where we use as function parameter the graph normalization.

~~~~{.cpp}
TGraph * g = new TGraph(npointx, xvec, yvec);
TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*g->Eval(x[0]); }, xmin, xmax, 1);
~~~~


### <a name="F6"></a> 6 - A member function with parameters of a general C++ class

A TF1 can be created in this case from any member function of a class which has the signature of (double * , double *) and returning a double.

Example:

~~~~{.cpp}
class  MyFunction {
 public:
   ...
   double Evaluate() (double *x, double *p) {
      // function implementation
   }
};
{
    ....
   MyFunction * fptr = new MyFunction(....);  // create the user function class
   TF1 * f = new TF1("f",fptr,&MyFunction::Evaluate,0,1,npar,"MyFunction","Evaluate");   // create TF1 class.

   .....
}
~~~~

See also the tutorial __math/exampleFunctor.C__ for a running example.
*/
////////////////////////////////////////////////////////////////////////////

TF1 *TF1::fgCurrent = 0;


////////////////////////////////////////////////////////////////////////////////
/// TF1 default constructor.

TF1::TF1():
   TNamed(), TAttLine(), TAttFill(), TAttMarker(),
   fXmin(0), fXmax(0), fNpar(0), fNdim(0), fType(EFType::kFormula)
{
   SetFillStyle(0);
}


////////////////////////////////////////////////////////////////////////////////
/// F1 constructor using a formula definition
///
/// See TFormula constructor for explanation of the formula syntax.
///
/// See tutorials: fillrandom, first, fit1, formula1, multifit
/// for real examples.
///
/// Creates a function of type A or B between xmin and xmax
///
/// if formula has the form "fffffff;xxxx;yyyy", it is assumed that
/// the formula string is "fffffff" and "xxxx" and "yyyy" are the
/// titles for the X and Y axis respectively.

TF1::TF1(const char *name, const char *formula, Double_t xmin, Double_t xmax, EAddToList addToGlobList, bool vectorize) :
   TNamed(name, formula), TAttLine(), TAttFill(), TAttMarker(), fType(EFType::kFormula)
{
   if (xmin < xmax) {
      fXmin      = xmin;
      fXmax      = xmax;
   } else {
      fXmin = xmax; //when called from TF2,TF3
      fXmax = xmin;
   }
   // Create rep formula (no need to add to gROOT list since we will add the TF1 object)

   // First check if we are making a convolution
   if (TString(formula, 5) == "CONV(" && formula[strlen(formula) - 1] == ')') {
      // Look for single ',' delimiter
      int delimPosition = -1;
      int parenCount = 0;
      for (unsigned int i = 5; i < strlen(formula) - 1; i++) {
         if (formula[i] == '(')
            parenCount++;
         else if (formula[i] == ')')
            parenCount--;
         else if (formula[i] == ',' && parenCount == 0) {
            if (delimPosition == -1)
               delimPosition = i;
            else
               Error("TF1", "CONV takes 2 arguments. Too many arguments found in : %s", formula);
         }
      }
      if (delimPosition == -1)
         Error("TF1", "CONV takes 2 arguments. Only one argument found in : %s", formula);

      // Having found the delimiter, define the first and second formulas
      TString formula1 = TString(TString(formula)(5, delimPosition - 5));
      TString formula2 = TString(TString(formula)(delimPosition + 1, strlen(formula) - 1 - (delimPosition + 1)));
      // remove spaces from these formulas
      formula1.ReplaceAll(' ', "");
      formula2.ReplaceAll(' ', "");

      TF1 *function1 = (TF1 *)(gROOT->GetListOfFunctions()->FindObject(formula1));
      if (function1 == nullptr)
         function1 = new TF1((const char *)formula1, (const char *)formula1, xmin, xmax);
      TF1 *function2 = (TF1 *)(gROOT->GetListOfFunctions()->FindObject(formula2));
      if (function2 == nullptr)
         function2 = new TF1((const char *)formula2, (const char *)formula2, xmin, xmax);

      // std::cout << "functions have been defined" << std::endl;

      TF1Convolution *conv = new TF1Convolution(function1, function2);

      // (note: currently ignoring `useFFT` option)
      fNpar = conv->GetNpar();
      fNdim = 1;                         // (note: may want to extend this in the future?)

      fType = EFType::kCompositionFcn;
      fComposition = std::unique_ptr<TF1AbsComposition>(conv);

      fParams = std::unique_ptr<TF1Parameters>(new TF1Parameters(fNpar)); // default to zeros (TF1Convolution has no GetParameters())
      // set parameter names
      for (int i = 0; i < fNpar; i++)
         this->SetParName(i, conv->GetParName(i));
      //  set parameters to default values
      int f1Npar = function1->GetNpar();
      int f2Npar = function2->GetNpar();
      // first, copy parameters from function1
      for (int i = 0; i < f1Npar; i++)
         this->SetParameter(i, function1->GetParameter(i));
      // then, check if the "Constant" parameters were combined
      // (this code assumes function2 has at most one parameter named "Constant")
      if (conv->GetNpar() == f1Npar + f2Npar - 1) {
         int cst1 = function1->GetParNumber("Constant");
         int cst2 = function2->GetParNumber("Constant");
         this->SetParameter(cst1, function1->GetParameter(cst1) * function2->GetParameter(cst2));
         // and copy parameters from function2
         for (int i = 0; i < f2Npar; i++)
            if (i < cst2)
               this->SetParameter(f1Npar + i, function2->GetParameter(i));
            else if (i > cst2)
               this->SetParameter(f1Npar + i - 1, function2->GetParameter(i));
      } else {
         // or if no constant, simply copy parameters from function2
         for (int i = 0; i < f2Npar; i++)
            this->SetParameter(i + f1Npar, function2->GetParameter(i));
      }

      // Then check if we need NSUM syntax:
   } else if (TString(formula, 5) == "NSUM(" && formula[strlen(formula) - 1] == ')') {
      // using comma as delimiter
      char delimiter = ',';
      // first, remove "NSUM(" and ")" and spaces
      TString formDense = TString(formula)(5,strlen(formula)-5-1);
      formDense.ReplaceAll(' ', "");

      // make sure standard functions are defined (e.g. gaus, expo)
      InitStandardFunctions();

      // Go char-by-char to split terms and define the relevant functions
      int parenCount = 0;
      int termStart = 0;
      TObjArray *newFuncs = new TObjArray();
      newFuncs->SetOwner(kTRUE);
      TObjArray *coeffNames = new TObjArray();
      coeffNames->SetOwner(kTRUE);
      TString fullFormula("");
      for (int i = 0; i < formDense.Length(); ++i) {
         if (formDense[i] == '(')
            parenCount++;
         else if (formDense[i] == ')')
            parenCount--;
         else if (formDense[i] == delimiter && parenCount == 0) {
            // term goes from termStart to i
            DefineNSUMTerm(newFuncs, coeffNames, fullFormula, formDense, termStart, i, xmin, xmax);
            termStart = i + 1;
         }
      }
      DefineNSUMTerm(newFuncs, coeffNames, fullFormula, formDense, termStart, formDense.Length(), xmin, xmax);

      TF1NormSum *normSum = new TF1NormSum(fullFormula, xmin, xmax);

      if (xmin == 0 && xmax == 1.) Info("TF1","Created TF1NormSum object using the default [0,1] range");

      fNpar = normSum->GetNpar();
      fNdim = 1; // (note: may want to extend functionality in the future)

      fType = EFType::kCompositionFcn;
      fComposition = std::unique_ptr<TF1AbsComposition>(normSum);

      fParams = std::unique_ptr<TF1Parameters>(new TF1Parameters(fNpar));
      fParams->SetParameters(&(normSum->GetParameters())[0]); // inherit default parameters from normSum

      // Parameter names
      for (int i = 0; i < fNpar; i++) {
         if (coeffNames->At(i) != nullptr) {
            TString coeffName = ((TObjString *)coeffNames->At(i))->GetString();
            this->SetParName(i, (const char *)coeffName);
         } else {
            this->SetParName(i, normSum->GetParName(i));
         }
      }

   } else { // regular TFormula
      fFormula = std::unique_ptr<TFormula>(new TFormula(name, formula, false, vectorize));
      fNpar = fFormula->GetNpar();
      // TFormula can have dimension zero, but since this is a TF1 minimal dim is 1
      fNdim = fFormula->GetNdim() == 0 ? 1 : fFormula->GetNdim();
   }
   if (fNpar) {
      fParErrors.resize(fNpar);
      fParMin.resize(fNpar);
      fParMax.resize(fNpar);
   }
   // do we want really to have this un-documented feature where we accept cases where dim > 1
   // by setting xmin >= xmax ??
   if (fNdim > 1 && xmin < xmax) {
      Error("TF1", "function: %s/%s has dimension %d instead of 1", name, formula, fNdim);
      MakeZombie();
   }

   DoInitialize(addToGlobList);
}
TF1::EAddToList GetGlobalListOption(Option_t * opt)  {
   if (opt == nullptr) return TF1::EAddToList::kDefault;
   TString option(opt);
   option.ToUpper();
   if (option.Contains("NL")) return TF1::EAddToList::kNo;
   if (option.Contains("GL")) return TF1::EAddToList::kAdd;
   return TF1::EAddToList::kDefault;
}
bool GetVectorizedOption(Option_t * opt)  {
   if (opt == nullptr) return false;
   TString option(opt);
   option.ToUpper();
   if (option.Contains("VEC")) return true;
   return false;
}
TF1::TF1(const char *name, const char *formula, Double_t xmin, Double_t xmax, Option_t * opt) :
////////////////////////////////////////////////////////////////////////////////
/// Same constructor as above (for TFormula based function) but passing an option strings
///  available options
///  VEC -  vectorize the formula expressions (not possible for lambda based expressions)
///  NL   - function is not stores in the global list of functions
///  GL   -  function will be always stored in the global list of functions ,
///         independently of the global setting of TF1::DefaultAddToGlobalList
///////////////////////////////////////////////////////////////////////////////////
   TF1(name, formula, xmin, xmax, GetGlobalListOption(opt), GetVectorizedOption(opt) )
{}
////////////////////////////////////////////////////////////////////////////////
/// F1 constructor using name of an interpreted function.
///
///  Creates a function of type C between xmin and xmax.
///  name is the name of an interpreted C++ function.
///  The function is defined with npar parameters
///  fcn must be a function of type:
///
///     Double_t fcn(Double_t *x, Double_t *params)
///
///  This constructor is called for functions of type C by the C++ interpreter.
///
/// WARNING! A function created with this constructor cannot be Cloned.

TF1::TF1(const char *name, Double_t xmin, Double_t xmax, Int_t npar, Int_t ndim, EAddToList addToGlobList) :
   TF1(EFType::kInterpreted, name, xmin, xmax, npar, ndim, addToGlobList, new TF1Parameters(npar))
{
   if (fName.Data()[0] == '*') {  // case TF1 name starts with a *
      Info("TF1", "TF1 has a name starting with a \'*\' - it is for saved TF1 objects in a .C file");
      return; //case happens via SavePrimitive
   } else if (fName.IsNull()) {
      Error("TF1", "requires a proper function name!");
      return;
   }

   fMethodCall = std::unique_ptr<TMethodCall>(new TMethodCall());
   fMethodCall->InitWithPrototype(fName, "Double_t*,Double_t*");

   if (! fMethodCall->IsValid()) {
      Error("TF1", "No function found with the signature %s(Double_t*,Double_t*)", name);
      return;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor using a pointer to a real function.
///
/// \param npar is the number of free parameters used by the function
///
/// This constructor creates a function of type C when invoked
/// with the normal C++ compiler.
///
/// see test program test/stress.cxx (function stress1) for an example.
/// note the interface with an intermediate pointer.
///
/// WARNING! A function created with this constructor cannot be Cloned.

TF1::TF1(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin, Double_t xmax, Int_t npar, Int_t ndim, EAddToList addToGlobList) :
   TF1(EFType::kPtrScalarFreeFcn, name, xmin, xmax, npar, ndim, addToGlobList, new TF1Parameters(npar), new TF1FunctorPointerImpl<double>(ROOT::Math::ParamFunctor(fcn)))
{}

////////////////////////////////////////////////////////////////////////////////
/// Constructor using a pointer to real function.
///
/// \param npar is the number of free parameters used by the function
///
/// This constructor creates a function of type C when invoked
/// with the normal C++ compiler.
///
/// see test program test/stress.cxx (function stress1) for an example.
/// note the interface with an intermediate pointer.
///
/// WARNING! A function created with this constructor cannot be Cloned.

TF1::TF1(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin, Double_t xmax, Int_t npar, Int_t ndim, EAddToList addToGlobList) :
   TF1(EFType::kPtrScalarFreeFcn, name, xmin, xmax, npar, ndim, addToGlobList, new TF1Parameters(npar), new TF1FunctorPointerImpl<double>(ROOT::Math::ParamFunctor(fcn)))
{}

////////////////////////////////////////////////////////////////////////////////
/// Constructor using the Functor class.
///
/// \param xmin and
/// \param xmax define the plotting range of the function
/// \param npar is the number of free parameters used by the function
///
/// This constructor can be used only in compiled code
///
/// WARNING! A function created with this constructor cannot be Cloned.

TF1::TF1(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin, Double_t xmax, Int_t npar, Int_t ndim, EAddToList addToGlobList) :
   TF1(EFType::kPtrScalarFreeFcn, name, xmin, xmax, npar, ndim, addToGlobList, new TF1Parameters(npar), new TF1FunctorPointerImpl<double>(ROOT::Math::ParamFunctor(f)))
{}

////////////////////////////////////////////////////////////////////////////////
/// Common initialization of the TF1. Add to the global list and
/// set the default style

void TF1::DoInitialize(EAddToList addToGlobalList)
{
   // add to global list of functions if default adding is on OR if bit is set
   bool doAdd = ((addToGlobalList == EAddToList::kDefault && fgAddToGlobList)
                 || addToGlobalList == EAddToList::kAdd);
   if (doAdd && gROOT) {
      SetBit(kNotGlobal, kFALSE);
      R__LOCKGUARD(gROOTMutex);
      // Store formula in linked list of formula in ROOT
      TF1 *f1old = (TF1 *)gROOT->GetListOfFunctions()->FindObject(fName);
      if (f1old) {
         gROOT->GetListOfFunctions()->Remove(f1old);
         // We removed f1old from the list, it is not longer global.
         // (See TF1::AddToGlobalList which requires this flag to be correct).
         f1old->SetBit(kNotGlobal, kTRUE);
      }
      gROOT->GetListOfFunctions()->Add(this);
   } else
      SetBit(kNotGlobal, kTRUE);

   if (gStyle) {
      SetLineColor(gStyle->GetFuncColor());
      SetLineWidth(gStyle->GetFuncWidth());
      SetLineStyle(gStyle->GetFuncStyle());
   }
   SetFillStyle(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to add/avoid to add automatically functions to the global list  (gROOT->GetListOfFunctions() )
/// After having called this static method, all the functions created afterwards will follow the
/// desired behaviour.
///
/// By default the functions are added automatically
/// It returns the previous status (true if the functions are added automatically)

Bool_t TF1::DefaultAddToGlobalList(Bool_t on)
{
   return fgAddToGlobList.exchange(on);
}

////////////////////////////////////////////////////////////////////////////////
/// Add to global list of functions (gROOT->GetListOfFunctions() )
/// return previous status (true if the function was already in the list false if not)

Bool_t TF1::AddToGlobalList(Bool_t on)
{
   if (!gROOT) return false;

   bool prevStatus = !TestBit(kNotGlobal);
   if (on)  {
      if (prevStatus) {
         R__LOCKGUARD(gROOTMutex);
         assert(gROOT->GetListOfFunctions()->FindObject(this) != nullptr);
         return on; // do nothing
      }
      // do I need to delete previous one with the same name ???
      //TF1 * old = dynamic_cast<TF1*>( gROOT->GetListOfFunctions()->FindObject(GetName()) );
      //if (old) { gROOT->GetListOfFunctions()->Remove(old); old->SetBit(kNotGlobal, kTRUE); }
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfFunctions()->Add(this);
      SetBit(kNotGlobal, kFALSE);
   } else if (prevStatus) {
      // if previous status was on and now is off we need to remove the function
      SetBit(kNotGlobal, kTRUE);
      R__LOCKGUARD(gROOTMutex);
      TF1 *old = dynamic_cast<TF1 *>(gROOT->GetListOfFunctions()->FindObject(GetName()));
      if (!old) {
         Warning("AddToGlobalList", "Function is supposed to be in the global list but it is not present");
         return kFALSE;
      }
      gROOT->GetListOfFunctions()->Remove(this);
   }
   return prevStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper functions for NSUM parsing

// Defines the formula that a given term uses, if not already defined,
// and appends "sanitized" formula to `fullFormula` string
void TF1::DefineNSUMTerm(TObjArray *newFuncs, TObjArray *coeffNames, TString &fullFormula, TString &formula,
                         int termStart, int termEnd, Double_t xmin, Double_t xmax)
{
   TString originalTerm = formula(termStart, termEnd-termStart);
   int coeffLength = TermCoeffLength(originalTerm);
   if (coeffLength != -1)
      termStart += coeffLength + 1;

   // `originalFunc` is the real formula and `cleanedFunc` is the
   // sanitized version that will not confuse the TF1NormSum
   // constructor
   TString originalFunc = formula(termStart, termEnd-termStart);
   TString cleanedFunc = TString(formula(termStart, termEnd-termStart))
      .ReplaceAll('+', "<plus>")
      .ReplaceAll('*',"<times>");

   // define function (if necessary)
   if (!gROOT->GetListOfFunctions()->FindObject(cleanedFunc))
      newFuncs->Add(new TF1(cleanedFunc, originalFunc, xmin, xmax));

   // append sanitized term to `fullFormula`
   if (fullFormula.Length() != 0)
      fullFormula.Append('+');

   // include numerical coefficient
   if (coeffLength != -1 && originalTerm[0] != '[')
      fullFormula.Append(originalTerm(0, coeffLength+1));

   // add coefficient name
   if (coeffLength != -1 && originalTerm[0] == '[')
      coeffNames->Add(new TObjString(TString(originalTerm(1,coeffLength-2))));
   else
      coeffNames->Add(nullptr);

   fullFormula.Append(cleanedFunc);
}


// Returns length of coeff at beginning of a given term, not counting the '*'
// Returns -1 if no coeff found
// Coeff can be either a number or parameter name
int TF1::TermCoeffLength(TString &term) {
  int firstAsterisk = term.First('*');
  if (firstAsterisk == -1) // no asterisk found
    return -1;

  if (TString(term(0,firstAsterisk)).IsFloat())
     return firstAsterisk;

  if (term[0] == '[' && term[firstAsterisk-1] == ']'
      && TString(term(1,firstAsterisk-2)).IsAlnum())
     return firstAsterisk;

  return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TF1 &TF1::operator=(const TF1 &rhs)
{
   if (this != &rhs) {
      rhs.Copy(*this);
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// TF1 default destructor.

TF1::~TF1()
{
   if (fHistogram) delete fHistogram;

   // this was before in TFormula destructor
   {
      R__LOCKGUARD(gROOTMutex);
      if (gROOT) gROOT->GetListOfFunctions()->Remove(this);
   }

   if (fParent) fParent->RecursiveRemove(this);

}


////////////////////////////////////////////////////////////////////////////////

TF1::TF1(const TF1 &f1) :
   TNamed(f1), TAttLine(f1), TAttFill(f1), TAttMarker(f1),
   fXmin(0), fXmax(0), fNpar(0), fNdim(0), fType(EFType::kFormula)
{
   ((TF1 &)f1).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Static function: set the fgAbsValue flag.
/// By default TF1::Integral uses the original function value to compute the integral
/// However, TF1::Moment, CentralMoment require to compute the integral
/// using the absolute value of the function.

void TF1::AbsValue(Bool_t flag)
{
   fgAbsValue = flag;
}


////////////////////////////////////////////////////////////////////////////////
/// Browse.

void TF1::Browse(TBrowser *b)
{
   Draw(b ? b->GetDrawOption() : "");
   gPad->Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this F1 to a new F1.
/// Note that the cached integral with its related arrays are not copied
/// (they are also set as transient data members)

void TF1::Copy(TObject &obj) const
{
   delete((TF1 &)obj).fHistogram;

   TNamed::Copy((TF1 &)obj);
   TAttLine::Copy((TF1 &)obj);
   TAttFill::Copy((TF1 &)obj);
   TAttMarker::Copy((TF1 &)obj);
   ((TF1 &)obj).fXmin = fXmin;
   ((TF1 &)obj).fXmax = fXmax;
   ((TF1 &)obj).fNpx  = fNpx;
   ((TF1 &)obj).fNpar = fNpar;
   ((TF1 &)obj).fNdim = fNdim;
   ((TF1 &)obj).fType = fType;
   ((TF1 &)obj).fChisquare = fChisquare;
   ((TF1 &)obj).fNpfits  = fNpfits;
   ((TF1 &)obj).fNDF     = fNDF;
   ((TF1 &)obj).fMinimum = fMinimum;
   ((TF1 &)obj).fMaximum = fMaximum;

   ((TF1 &)obj).fParErrors = fParErrors;
   ((TF1 &)obj).fParMin    = fParMin;
   ((TF1 &)obj).fParMax    = fParMax;
   ((TF1 &)obj).fParent    = fParent;
   ((TF1 &)obj).fSave      = fSave;
   ((TF1 &)obj).fHistogram = 0;
   ((TF1 &)obj).fMethodCall = 0;
   ((TF1 &)obj).fNormalized = fNormalized;
   ((TF1 &)obj).fNormIntegral = fNormIntegral;
   ((TF1 &)obj).fFormula   = 0;

   if (fFormula) assert(fFormula->GetNpar() == fNpar);

   // use copy-constructor of TMethodCall
   TMethodCall *m = (fMethodCall) ? new TMethodCall(*fMethodCall) : nullptr;
   ((TF1 &)obj).fMethodCall.reset(m);

   TFormula *formulaToCopy = (fFormula) ? new TFormula(*fFormula) : nullptr;
   ((TF1 &)obj).fFormula.reset(formulaToCopy);

   TF1Parameters *paramsToCopy = (fParams) ? new TF1Parameters(*fParams) : nullptr;
   ((TF1 &)obj).fParams.reset(paramsToCopy);

   TF1FunctorPointer *functorToCopy = (fFunctor) ? fFunctor->Clone() : nullptr;
   ((TF1 &)obj).fFunctor.reset(functorToCopy);

   TF1AbsComposition *comp = nullptr;
   if (fComposition) {
      comp = (TF1AbsComposition *)fComposition->IsA()->New();
      fComposition->Copy(*comp);
   }
   ((TF1 &)obj).fComposition.reset(comp);
}


////////////////////////////////////////////////////////////////////////////////
/// Make a complete copy of the underlying object.  If 'newname' is set,
/// the copy's name will be set to that name.

TObject* TF1::Clone(const char* newname) const
{

   TF1* obj = (TF1*) TNamed::Clone(newname);

   if (fHistogram) {
      obj->fHistogram = (TH1*)fHistogram->Clone();
      obj->fHistogram->SetDirectory(0);
   }

   return obj;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the first derivative of the function at point x,
/// computed by Richardson's extrapolation method (use 2 derivative estimates
/// to compute a third, more accurate estimation)
/// first, derivatives with steps h and h/2 are computed by central difference formulas
/// \f[
///   D(h) = \frac{f(x+h) - f(x-h)}{2h}
/// \f]
/// the final estimate
/// \f[
///      D = \frac{4D(h/2) - D(h)}{3}
/// \f]
/// "Numerical Methods for Scientists and Engineers", H.M.Antia, 2nd edition"
///
/// if the argument params is null, the current function parameters are used,
/// otherwise the parameters in params are used.
///
/// the argument eps may be specified to control the step size (precision).
/// the step size is taken as eps*(xmax-xmin).
/// the default value (0.001) should be good enough for the vast majority
/// of functions. Give a smaller value if your function has many changes
/// of the second derivative in the function range.
///
/// Getting the error via TF1::DerivativeError:
///   (total error = roundoff error + interpolation error)
/// the estimate of the roundoff error is taken as follows:
/// \f[
///    err = k\sqrt{f(x)^{2} + x^{2}deriv^{2}}\sqrt{\sum ai^{2}},
/// \f]
/// where k is the double precision, ai are coefficients used in
/// central difference formulas
/// interpolation error is decreased by making the step size h smaller.
///
/// \author Anna Kreshuk

Double_t TF1::Derivative(Double_t x, Double_t *params, Double_t eps) const
{
   if (GetNdim() > 1) {
      Warning("Derivative", "Function dimension is larger than one");
   }

   ROOT::Math::RichardsonDerivator rd;
   double xmin, xmax;
   GetRange(xmin, xmax);
   // this is not optimal (should be used the average x instead of the range)
   double h = eps * std::abs(xmax - xmin);
   if (h <= 0) h = 0.001;
   double der = 0;
   if (params) {
      ROOT::Math::WrappedTF1 wtf(*(const_cast<TF1 *>(this)));
      wtf.SetParameters(params);
      der = rd.Derivative1(wtf, x, h);
   } else {
      // no need to set parameters used a non-parametric wrapper to avoid allocating
      // an array with parameter values
      ROOT::Math::WrappedFunction<const TF1 & > wf(*this);
      der = rd.Derivative1(wf, x, h);
   }

   gErrorTF1 = rd.Error();
   return der;

}


////////////////////////////////////////////////////////////////////////////////
/// Returns the second derivative of the function at point x,
/// computed by Richardson's extrapolation method (use 2 derivative estimates
/// to compute a third, more accurate estimation)
/// first, derivatives with steps h and h/2 are computed by central difference formulas
/// \f[
///    D(h) = \frac{f(x+h) - 2f(x) + f(x-h)}{h^{2}}
/// \f]
/// the final estimate
/// \f[
///    D = \frac{4D(h/2) - D(h)}{3}
/// \f]
///  "Numerical Methods for Scientists and Engineers", H.M.Antia, 2nd edition"
///
/// if the argument params is null, the current function parameters are used,
/// otherwise the parameters in params are used.
///
/// the argument eps may be specified to control the step size (precision).
/// the step size is taken as eps*(xmax-xmin).
/// the default value (0.001) should be good enough for the vast majority
/// of functions. Give a smaller value if your function has many changes
/// of the second derivative in the function range.
///
/// Getting the error via TF1::DerivativeError:
///   (total error = roundoff error + interpolation error)
/// the estimate of the roundoff error is taken as follows:
/// \f[
///    err = k\sqrt{f(x)^{2} + x^{2}deriv^{2}}\sqrt{\sum ai^{2}},
/// \f]
/// where k is the double precision, ai are coefficients used in
/// central difference formulas
/// interpolation error is decreased by making the step size h smaller.
///
/// \author Anna Kreshuk

Double_t TF1::Derivative2(Double_t x, Double_t *params, Double_t eps) const
{
   if (GetNdim() > 1) {
      Warning("Derivative2", "Function dimension is larger than one");
   }

   ROOT::Math::RichardsonDerivator rd;
   double xmin, xmax;
   GetRange(xmin, xmax);
   // this is not optimal (should be used the average x instead of the range)
   double h = eps * std::abs(xmax - xmin);
   if (h <= 0) h = 0.001;
   double der = 0;
   if (params) {
      ROOT::Math::WrappedTF1 wtf(*(const_cast<TF1 *>(this)));
      wtf.SetParameters(params);
      der = rd.Derivative2(wtf, x, h);
   } else {
      // no need to set parameters used a non-parametric wrapper to avoid allocating
      // an array with parameter values
      ROOT::Math::WrappedFunction<const TF1 & > wf(*this);
      der = rd.Derivative2(wf, x, h);
   }

   gErrorTF1 = rd.Error();

   return der;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the third derivative of the function at point x,
/// computed by Richardson's extrapolation method (use 2 derivative estimates
/// to compute a third, more accurate estimation)
/// first, derivatives with steps h and h/2 are computed by central difference formulas
/// \f[
///    D(h) = \frac{f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)}{2h^{3}}
/// \f]
/// the final estimate
/// \f[
///    D = \frac{4D(h/2) - D(h)}{3}
/// \f]
///  "Numerical Methods for Scientists and Engineers", H.M.Antia, 2nd edition"
///
/// if the argument params is null, the current function parameters are used,
/// otherwise the parameters in params are used.
///
/// the argument eps may be specified to control the step size (precision).
/// the step size is taken as eps*(xmax-xmin).
/// the default value (0.001) should be good enough for the vast majority
/// of functions. Give a smaller value if your function has many changes
/// of the second derivative in the function range.
///
/// Getting the error via TF1::DerivativeError:
///   (total error = roundoff error + interpolation error)
/// the estimate of the roundoff error is taken as follows:
/// \f[
///    err = k\sqrt{f(x)^{2} + x^{2}deriv^{2}}\sqrt{\sum ai^{2}},
/// \f]
/// where k is the double precision, ai are coefficients used in
/// central difference formulas
/// interpolation error is decreased by making the step size h smaller.
///
/// \author Anna Kreshuk

Double_t TF1::Derivative3(Double_t x, Double_t *params, Double_t eps) const
{
   if (GetNdim() > 1) {
      Warning("Derivative3", "Function dimension is larger than one");
   }

   ROOT::Math::RichardsonDerivator rd;
   double xmin, xmax;
   GetRange(xmin, xmax);
   // this is not optimal (should be used the average x instead of the range)
   double h = eps * std::abs(xmax - xmin);
   if (h <= 0) h = 0.001;
   double der = 0;
   if (params) {
      ROOT::Math::WrappedTF1 wtf(*(const_cast<TF1 *>(this)));
      wtf.SetParameters(params);
      der = rd.Derivative3(wtf, x, h);
   } else {
      // no need to set parameters used a non-parametric wrapper to avoid allocating
      // an array with parameter values
      ROOT::Math::WrappedFunction<const TF1 & > wf(*this);
      der = rd.Derivative3(wf, x, h);
   }

   gErrorTF1 = rd.Error();
   return der;

}


////////////////////////////////////////////////////////////////////////////////
/// Static function returning the error of the last call to the of Derivative's
/// functions

Double_t TF1::DerivativeError()
{
   return gErrorTF1;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a function.
///
/// Compute the closest distance of approach from point px,py to this
/// function. The distance is computed in pixels units.
///
/// Note that px is called with a negative value when the TF1 is in
/// TGraph or TH1 list of functions. In this case there is no point
/// looking at the histogram axis.

Int_t TF1::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (!fHistogram) return 9999;
   Int_t distance = 9999;
   if (px >= 0) {
      distance = fHistogram->DistancetoPrimitive(px, py);
      if (distance <= 1) return distance;
   } else {
      px = -px;
   }

   Double_t xx[1];
   Double_t x    = gPad->AbsPixeltoX(px);
   xx[0]         = gPad->PadtoX(x);
   if (xx[0] < fXmin || xx[0] > fXmax) return distance;
   Double_t fval = Eval(xx[0]);
   Double_t y    = gPad->YtoPad(fval);
   Int_t pybin   = gPad->YtoAbsPixel(y);
   return TMath::Abs(py - pybin);
}


////////////////////////////////////////////////////////////////////////////////
/// Draw this function with its current attributes.
///
/// Possible option values are:
///
/// option | description
/// -------|----------------------------------------
/// "SAME" | superimpose on top of existing picture
/// "L"    | connect all computed points with a straight line
/// "C"    | connect all computed points with a smooth curve
/// "FC"   | draw a fill area below a smooth curve
///
/// Note that the default value is "L". Therefore to draw on top
/// of an existing picture, specify option "LSAME"
///
/// NB. You must use DrawCopy if you want to draw several times the same
///     function in the current canvas.

void TF1::Draw(Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   if (gPad && !opt.Contains("same")) gPad->Clear();

   AppendPad(option);

   gPad->IncrementPaletteColor(1, opt);
}


////////////////////////////////////////////////////////////////////////////////
/// Draw a copy of this function with its current attributes.
///
///  This function MUST be used instead of Draw when you want to draw
///  the same function with different parameters settings in the same canvas.
///
/// Possible option values are:
///
/// option | description
/// -------|----------------------------------------
/// "SAME" | superimpose on top of existing picture
/// "L"    | connect all computed points with a straight line
/// "C"    | connect all computed points with a smooth curve
/// "FC"   | draw a fill area below a smooth curve
///
/// Note that the default value is "L". Therefore to draw on top
/// of an existing picture, specify option "LSAME"

TF1 *TF1::DrawCopy(Option_t *option) const
{
   TF1 *newf1 = (TF1 *)this->IsA()->New();
   Copy(*newf1);
   newf1->AppendPad(option);
   newf1->SetBit(kCanDelete);
   return newf1;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw derivative of this function
///
/// An intermediate TGraph object is built and drawn with option.
/// The function returns a pointer to the TGraph object. Do:
///
///     TGraph *g = (TGraph*)myfunc.DrawDerivative(option);
///
/// The resulting graph will be drawn into the current pad.
/// If this function is used via the context menu, it recommended
/// to create a new canvas/pad before invoking this function.

TObject *TF1::DrawDerivative(Option_t *option)
{
   TVirtualPad *pad = gROOT->GetSelectedPad();
   TVirtualPad *padsav = gPad;
   if (pad) pad->cd();

   TGraph *gr = new TGraph(this, "d");
   gr->Draw(option);
   if (padsav) padsav->cd();
   return gr;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw integral of this function
///
/// An intermediate TGraph object is built and drawn with option.
/// The function returns a pointer to the TGraph object. Do:
///
///     TGraph *g = (TGraph*)myfunc.DrawIntegral(option);
///
/// The resulting graph will be drawn into the current pad.
/// If this function is used via the context menu, it recommended
/// to create a new canvas/pad before invoking this function.

TObject *TF1::DrawIntegral(Option_t *option)
{
   TVirtualPad *pad = gROOT->GetSelectedPad();
   TVirtualPad *padsav = gPad;
   if (pad) pad->cd();

   TGraph *gr = new TGraph(this, "i");
   gr->Draw(option);
   if (padsav) padsav->cd();
   return gr;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw function between xmin and xmax.

void TF1::DrawF1(Double_t xmin, Double_t xmax, Option_t *option)
{
//    //if(Compile(formula)) return ;
   SetRange(xmin, xmax);

   Draw(option);
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate this function.
///
/// Computes the value of this function (general case for a 3-d function)
/// at point x,y,z.
/// For a 1-d function give y=0 and z=0
/// The current value of variables x,y,z is passed through x, y and z.
/// The parameters used will be the ones in the array params if params is given
/// otherwise parameters will be taken from the stored data members fParams

Double_t TF1::Eval(Double_t x, Double_t y, Double_t z, Double_t t) const
{
   if (fType == EFType::kFormula) return fFormula->Eval(x, y, z, t);

   Double_t xx[4] = {x, y, z, t};
   Double_t *pp = (Double_t *)fParams->GetParameters();
   // if (fType == EFType::kInterpreted)((TF1 *)this)->InitArgs(xx, pp);
   return ((TF1 *)this)->EvalPar(xx, pp);
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate function with given coordinates and parameters.
///
/// Compute the value of this function at point defined by array x
/// and current values of parameters in array params.
/// If argument params is omitted or equal 0, the internal values
/// of parameters (array fParams) will be used instead.
/// For a 1-D function only x[0] must be given.
/// In case of a multi-dimensional function, the arrays x must be
/// filled with the corresponding number of dimensions.
///
/// WARNING. In case of an interpreted function (fType=2), it is the
/// user's responsibility to initialize the parameters via InitArgs
/// before calling this function.
/// InitArgs should be called at least once to specify the addresses
/// of the arguments x and params.
/// InitArgs should be called every time these addresses change.

Double_t TF1::EvalPar(const Double_t *x, const Double_t *params)
{
   //fgCurrent = this;

   if (fType == EFType::kFormula) {
      assert(fFormula);

      if (fNormalized && fNormIntegral != 0)
         return fFormula->EvalPar(x, params) / fNormIntegral;
      else
         return fFormula->EvalPar(x, params);
   }
   Double_t result = 0;
   if (fType == EFType::kPtrScalarFreeFcn || fType == EFType::kTemplScalar)  {
      if (fFunctor) {
         assert(fParams);
         if (params) result = ((TF1FunctorPointerImpl<Double_t> *)fFunctor.get())->fImpl((Double_t *)x, (Double_t *)params);
         else        result = ((TF1FunctorPointerImpl<Double_t> *)fFunctor.get())->fImpl((Double_t *)x, (Double_t *)fParams->GetParameters());

      } else          result = GetSave(x);

      if (fNormalized && fNormIntegral != 0)
         result = result / fNormIntegral;

      return result;
   }
   if (fType == EFType::kInterpreted) {
      if (fMethodCall) fMethodCall->Execute(result);
      else             result = GetSave(x);

      if (fNormalized && fNormIntegral != 0)
         result = result / fNormIntegral;

      return result;
   }

#ifdef R__HAS_VECCORE
   if (fType == EFType::kTemplVec) {
      if (fFunctor) {
         if (params) result =  EvalParVec(x, params);
         else result =  EvalParVec(x, (Double_t *) fParams->GetParameters());
      }
      else {
         result = GetSave(x);
      }

      if (fNormalized && fNormIntegral != 0)
         result = result / fNormIntegral;

      return result;
   }
#endif

   if (fType == EFType::kCompositionFcn) {
      if (!fComposition)
         Error("EvalPar", "Composition function not found");

      result = (*fComposition)(x, params);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when a F1 is clicked with the locator

void TF1::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   if (fHistogram) fHistogram->ExecuteEvent(event, px, py);

   if (!gPad->GetView()) {
      if (event == kMouseMotion)  gPad->SetCursor(kHand);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Fix the value of a parameter
/// The specified value will be used in a fit operation

void TF1::FixParameter(Int_t ipar, Double_t value)
{
   if (ipar < 0 || ipar > GetNpar() - 1) return;
   SetParameter(ipar, value);
   if (value != 0) SetParLimits(ipar, value, value);
   else            SetParLimits(ipar, 1, 1);
}


////////////////////////////////////////////////////////////////////////////////
/// Static function returning the current function being processed

TF1 *TF1::GetCurrent()
{
   ::Warning("TF1::GetCurrent", "This function is obsolete and is working only for the current painted functions");
   return fgCurrent;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the histogram used to visualise the function
/// Note that this histogram is managed by the function and
/// in same case it is automatically deleted when some TF1 functions are called
///  such as TF1::SetParameters, TF1::SetNpx, TF1::SetRange
/// It is then reccomended either to clone the return object or calling again teh GetHistogram
/// function whenever is needed

TH1 *TF1::GetHistogram() const
{
   if (fHistogram) return fHistogram;

   // histogram has not been yet created - create it
   // should not we make this function not const ??
   const_cast<TF1 *>(this)->fHistogram  = const_cast<TF1 *>(this)->CreateHistogram();
   if (!fHistogram) Error("GetHistogram", "Error creating histogram for function %s of type %s", GetName(), IsA()->GetName());
   return fHistogram;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum value of the function
///
/// Method:
///  First, the grid search is used to bracket the maximum
///  with the step size = (xmax-xmin)/fNpx.
///  This way, the step size can be controlled via the SetNpx() function.
///  If the function is unimodal or if its extrema are far apart, setting
///  the fNpx to a small value speeds the algorithm up many times.
///  Then, Brent's method is applied on the bracketed interval
///  epsilon (default = 1.E-10) controls the relative accuracy (if |x| > 1 )
///  and absolute (if |x| < 1)  and maxiter (default = 100) controls the maximum number
///  of iteration of the Brent algorithm
///  If the flag logx is set the grid search is done in log step size
///  This is done automatically if the log scale is set in the current Pad
///
/// NOTE: see also TF1::GetMaximumX and TF1::GetX

Double_t TF1::GetMaximum(Double_t xmin, Double_t xmax, Double_t epsilon, Int_t maxiter, Bool_t logx) const
{
   if (xmin >= xmax) {
      xmin = fXmin;
      xmax = fXmax;
   }

   if (!logx && gPad != 0) logx = gPad->GetLogx();

   ROOT::Math::BrentMinimizer1D bm;
   GInverseFunc g(this);
   ROOT::Math::WrappedFunction<GInverseFunc> wf1(g);
   bm.SetFunction(wf1, xmin, xmax);
   bm.SetNpx(fNpx);
   bm.SetLogScan(logx);
   bm.Minimize(maxiter, epsilon, epsilon);
   Double_t x;
   x = - bm.FValMinimum();

   return x;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the X value corresponding to the maximum value of the function
///
/// Method:
///  First, the grid search is used to bracket the maximum
///  with the step size = (xmax-xmin)/fNpx.
///  This way, the step size can be controlled via the SetNpx() function.
///  If the function is unimodal or if its extrema are far apart, setting
///  the fNpx to a small value speeds the algorithm up many times.
///  Then, Brent's method is applied on the bracketed interval
///  epsilon (default = 1.E-10) controls the relative accuracy (if |x| > 1 )
///  and absolute (if |x| < 1)  and maxiter (default = 100) controls the maximum number
///  of iteration of the Brent algorithm
///  If the flag logx is set the grid search is done in log step size
///  This is done automatically if the log scale is set in the current Pad
///
/// NOTE: see also TF1::GetX

Double_t TF1::GetMaximumX(Double_t xmin, Double_t xmax, Double_t epsilon, Int_t maxiter, Bool_t logx) const
{
   if (xmin >= xmax) {
      xmin = fXmin;
      xmax = fXmax;
   }

   if (!logx && gPad != 0) logx = gPad->GetLogx();

   ROOT::Math::BrentMinimizer1D bm;
   GInverseFunc g(this);
   ROOT::Math::WrappedFunction<GInverseFunc> wf1(g);
   bm.SetFunction(wf1, xmin, xmax);
   bm.SetNpx(fNpx);
   bm.SetLogScan(logx);
   bm.Minimize(maxiter, epsilon, epsilon);
   Double_t x;
   x = bm.XMinimum();

   return x;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the minimum value of the function on the (xmin, xmax) interval
///
/// Method:
///  First, the grid search is used to bracket the maximum
///  with the step size = (xmax-xmin)/fNpx. This way, the step size
///  can be controlled via the SetNpx() function. If the function is
///  unimodal or if its extrema are far apart, setting the fNpx to
///  a small value speeds the algorithm up many times.
///  Then, Brent's method is applied on the bracketed interval
///  epsilon (default = 1.E-10) controls the relative accuracy (if |x| > 1 )
///  and absolute (if |x| < 1)  and maxiter (default = 100) controls the maximum number
///  of iteration of the Brent algorithm
///  If the flag logx is set the grid search is done in log step size
///  This is done automatically if the log scale is set in the current Pad
///
/// NOTE: see also TF1::GetMaximumX and TF1::GetX

Double_t TF1::GetMinimum(Double_t xmin, Double_t xmax, Double_t epsilon, Int_t maxiter, Bool_t logx) const
{
   if (xmin >= xmax) {
      xmin = fXmin;
      xmax = fXmax;
   }

   if (!logx && gPad != 0) logx = gPad->GetLogx();

   ROOT::Math::BrentMinimizer1D bm;
   ROOT::Math::WrappedFunction<const TF1 &> wf1(*this);
   bm.SetFunction(wf1, xmin, xmax);
   bm.SetNpx(fNpx);
   bm.SetLogScan(logx);
   bm.Minimize(maxiter, epsilon, epsilon);
   Double_t x;
   x = bm.FValMinimum();

   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the minimum of a function of whatever dimension.
/// While GetMinimum works only for 1D function , GetMinimumNDim works for all dimensions
/// since it uses the minimizer interface
/// vector x at beginning will contained the initial point, on exit will contain the result

Double_t TF1::GetMinMaxNDim(Double_t *x , bool findmax, Double_t epsilon, Int_t maxiter) const
{
   R__ASSERT(x != 0);

   int ndim = GetNdim();
   if (ndim == 0) {
      Error("GetMinimumNDim", "Function of dimension 0 - return Eval(x)");
      return (const_cast<TF1 &>(*this))(x);
   }

   // create minimizer class
   const char *minimName = ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
   const char *minimAlgo = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
   ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(minimName, minimAlgo);

   if (min == 0) {
      Error("GetMinimumNDim", "Error creating minimizer %s", minimName);
      return 0;
   }

   // minimizer will be set using default values
   if (epsilon > 0) min->SetTolerance(epsilon);
   if (maxiter > 0) min->SetMaxFunctionCalls(maxiter);

   // create wrapper class from TF1 (cannot use Functor, t.b.i.)
   ROOT::Math::WrappedMultiFunction<TF1 &> objFunc(const_cast<TF1 &>(*this), ndim);
   // create -f(x) when searching for the maximum
   GInverseFuncNdim invFunc(const_cast<TF1 *>(this));
   ROOT::Math::WrappedMultiFunction<GInverseFuncNdim &> objFuncInv(invFunc, ndim);
   if (!findmax)
      min->SetFunction(objFunc);
   else
      min->SetFunction(objFuncInv);

   std::vector<double> rmin(ndim);
   std::vector<double> rmax(ndim);
   GetRange(&rmin[0], &rmax[0]);
   for (int i = 0; i < ndim; ++i) {
      const char *xname =  0;
      double stepSize = 0.1;
      // use range for step size or give some value depending on x if range is not defined
      if (rmax[i] > rmin[i])
         stepSize = (rmax[i] - rmin[i]) / 100;
      else if (std::abs(x[i]) > 1.)
         stepSize = 0.1 * x[i];

      // set variable names
      if (ndim <= 3) {
         if (i == 0) {
            xname = "x";
         } else if (i == 1) {
            xname = "y";
         } else {
            xname = "z";
         }
      } else {
         xname = TString::Format("x_%d", i);
         // arbitrary step sie (should be computed from range)
      }

      if (rmin[i] < rmax[i]) {
         //Info("GetMinMax","setting limits on %s - [ %f , %f ]",xname,rmin[i],rmax[i]);
         min->SetLimitedVariable(i, xname, x[i], stepSize, rmin[i], rmax[i]);
      } else {
         min->SetVariable(i, xname, x[i], stepSize);
      }
   }

   bool ret = min->Minimize();
   if (!ret) {
      Error("GetMinimumNDim", "Error minimizing function %s", GetName());
   }
   if (min->X()) std::copy(min->X(), min->X() + ndim, x);
   double fmin = min->MinValue();
   delete min;
   // need to revert sign in case looking for maximum
   return (findmax) ? -fmin : fmin;

}


////////////////////////////////////////////////////////////////////////////////
/// Returns the X value corresponding to the minimum value of the function
/// on the (xmin, xmax) interval
///
/// Method:
///  First, the grid search is used to bracket the maximum
///  with the step size = (xmax-xmin)/fNpx. This way, the step size
///  can be controlled via the SetNpx() function. If the function is
///  unimodal or if its extrema are far apart, setting the fNpx to
///  a small value speeds the algorithm up many times.
///  Then, Brent's method is applied on the bracketed interval
///  epsilon (default = 1.E-10) controls the relative accuracy (if |x| > 1 )
///  and absolute (if |x| < 1)  and maxiter (default = 100) controls the maximum number
///  of iteration of the Brent algorithm
///  If the flag logx is set the grid search is done in log step size
///  This is done automatically if the log scale is set in the current Pad
///
/// NOTE: see also TF1::GetX

Double_t TF1::GetMinimumX(Double_t xmin, Double_t xmax, Double_t epsilon, Int_t maxiter, Bool_t logx) const
{
   if (xmin >= xmax) {
      xmin = fXmin;
      xmax = fXmax;
   }

   ROOT::Math::BrentMinimizer1D bm;
   ROOT::Math::WrappedFunction<const TF1 &> wf1(*this);
   bm.SetFunction(wf1, xmin, xmax);
   bm.SetNpx(fNpx);
   bm.SetLogScan(logx);
   bm.Minimize(maxiter, epsilon, epsilon);
   Double_t x;
   x = bm.XMinimum();

   return x;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the X value corresponding to the function value fy for (xmin<x<xmax).
/// in other words it can find the roots of the function when fy=0 and successive calls
/// by changing the next call to [xmin+eps,xmax] where xmin is the previous root.
///
/// Method:
///  First, the grid search is used to bracket the maximum
///  with the step size = (xmax-xmin)/fNpx. This way, the step size
///  can be controlled via the SetNpx() function. If the function is
///  unimodal or if its extrema are far apart, setting the fNpx to
///  a small value speeds the algorithm up many times.
///  Then, Brent's method is applied on the bracketed interval
///  epsilon (default = 1.E-10) controls the relative accuracy (if |x| > 1 )
///  and absolute (if |x| < 1)  and maxiter (default = 100) controls the maximum number
///  of iteration of the Brent algorithm
///  If the flag logx is set the grid search is done in log step size
///  This is done automatically if the log scale is set in the current Pad
///
/// NOTE: see also TF1::GetMaximumX, TF1::GetMinimumX

Double_t TF1::GetX(Double_t fy, Double_t xmin, Double_t xmax, Double_t epsilon, Int_t maxiter, Bool_t logx) const
{
   if (xmin >= xmax) {
      xmin = fXmin;
      xmax = fXmax;
   }

   if (!logx && gPad != 0) logx = gPad->GetLogx();

   GFunc g(this, fy);
   ROOT::Math::WrappedFunction<GFunc> wf1(g);
   ROOT::Math::BrentRootFinder brf;
   brf.SetFunction(wf1, xmin, xmax);
   brf.SetNpx(fNpx);
   brf.SetLogScan(logx);
   bool ret = brf.Solve(maxiter, epsilon, epsilon);
   if (!ret) Error("GetX","[%f,%f] is not a valid interval",xmin,xmax);
   return (ret) ? brf.Root() : TMath::QuietNaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of degrees of freedom in the fit
/// the fNDF parameter has been previously computed during a fit.
/// The number of degrees of freedom corresponds to the number of points
/// used in the fit minus the number of free parameters.

Int_t TF1::GetNDF() const
{
   Int_t npar = GetNpar();
   if (fNDF == 0 && (fNpfits > npar)) return fNpfits - npar;
   return fNDF;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the number of free parameters

Int_t TF1::GetNumberFreeParameters() const
{
   Int_t ntot = GetNpar();
   Int_t nfree = ntot;
   Double_t al, bl;
   for (Int_t i = 0; i < ntot; i++) {
      ((TF1 *)this)->GetParLimits(i, al, bl);
      if (al * bl != 0 && al >= bl) nfree--;
   }
   return nfree;
}


////////////////////////////////////////////////////////////////////////////////
/// Redefines TObject::GetObjectInfo.
/// Displays the function info (x, function value)
/// corresponding to cursor position px,py

char *TF1::GetObjectInfo(Int_t px, Int_t /* py */) const
{
   static char info[64];
   Double_t x = gPad->PadtoX(gPad->AbsPixeltoX(px));
   snprintf(info, 64, "(x=%g, f=%g)", x, ((TF1 *)this)->Eval(x));
   return info;
}


////////////////////////////////////////////////////////////////////////////////
/// Return value of parameter number ipar

Double_t TF1::GetParError(Int_t ipar) const
{
   if (ipar < 0 || ipar > GetNpar() - 1) return 0;
   return fParErrors[ipar];
}


////////////////////////////////////////////////////////////////////////////////
/// Return limits for parameter ipar.

void TF1::GetParLimits(Int_t ipar, Double_t &parmin, Double_t &parmax) const
{
   parmin = 0;
   parmax = 0;
   int n = fParMin.size();
   assert(n == int(fParMax.size()) && n <= fNpar);
   if (ipar < 0 || ipar > n - 1) return;
   parmin = fParMin[ipar];
   parmax = fParMax[ipar];
}


////////////////////////////////////////////////////////////////////////////////
/// Return the fit probability

Double_t TF1::GetProb() const
{
   if (fNDF <= 0) return 0;
   return TMath::Prob(fChisquare, fNDF);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute Quantiles for density distribution of this function
///
/// Quantile x_q of a probability distribution Function F is defined as
/// \f[
///     F(x_{q}) = \int_{xmin}^{x_{q}} f dx = q with 0 <= q <= 1.
/// \f]
/// For instance the median \f$ x_{\frac{1}{2}} \f$ of a distribution is defined as that value
/// of the random variable for which the distribution function equals 0.5:
/// \f[
///        F(x_{\frac{1}{2}}) = \prod(x < x_{\frac{1}{2}}) = \frac{1}{2}
/// \f]
///
/// \param[in] nprobSum maximum size of array q and size of array probSum
/// \param[out] q array filled with nq quantiles
/// \param[in] probSum array of positions where quantiles will be computed.
///     It is assumed to contain at least nprobSum values.
/// \return value nq (<=nprobSum) with the number of quantiles computed
///
///  Getting quantiles from two histograms and storing results in a TGraph,
///  a so-called QQ-plot
///
///      TGraph *gr = new TGraph(nprob);
///      f1->GetQuantiles(nprob,gr->GetX());
///      f2->GetQuantiles(nprob,gr->GetY());
///      gr->Draw("alp");
///
/// \author Eddy Offermann


Int_t TF1::GetQuantiles(Int_t nprobSum, Double_t *q, const Double_t *probSum)
{
   // LM: change to use fNpx
   // should we change code to use a root finder ?
   // It should be more precise and more efficient
   const Int_t npx     = TMath::Max(fNpx, 2 * nprobSum);
   const Double_t xMin = GetXmin();
   const Double_t xMax = GetXmax();
   const Double_t dx   = (xMax - xMin) / npx;

   TArrayD integral(npx + 1);
   TArrayD alpha(npx);
   TArrayD beta(npx);
   TArrayD gamma(npx);

   integral[0] = 0;
   Int_t intNegative = 0;
   Int_t i;
   for (i = 0; i < npx; i++) {
      Double_t integ = Integral(Double_t(xMin + i * dx), Double_t(xMin + i * dx + dx), 0.0);
      if (integ < 0) {
         intNegative++;
         integ = -integ;
      }
      integral[i + 1] = integral[i] + integ;
   }

   if (intNegative > 0)
      Warning("GetQuantiles", "function:%s has %d negative values: abs assumed",
              GetName(), intNegative);
   if (integral[npx] == 0) {
      Error("GetQuantiles", "Integral of function is zero");
      return 0;
   }

   const Double_t total = integral[npx];
   for (i = 1; i <= npx; i++) integral[i] /= total;
   //the integral r for each bin is approximated by a parabola
   //  x = alpha + beta*r +gamma*r**2
   // compute the coefficients alpha, beta, gamma for each bin
   for (i = 0; i < npx; i++) {
      const Double_t x0 = xMin + dx * i;
      const Double_t r2 = integral[i + 1] - integral[i];
      const Double_t r1 = Integral(x0, x0 + 0.5 * dx, 0.0) / total;
      gamma[i] = (2 * r2 - 4 * r1) / (dx * dx);
      beta[i]  = r2 / dx - gamma[i] * dx;
      alpha[i] = x0;
      gamma[i] *= 2;
   }

   // Be careful because of finite precision in the integral; Use the fact that the integral
   // is monotone increasing
   for (i = 0; i < nprobSum; i++) {
      const Double_t r = probSum[i];
      Int_t bin  = TMath::Max(TMath::BinarySearch(npx + 1, integral.GetArray(), r), (Long64_t)0);
      // in case the prob is 1
      if (bin == npx) {
         q[i] = xMax;
         continue;
      }
      // LM use a tolerance 1.E-12 (integral precision)
      while (bin < npx - 1 && TMath::AreEqualRel(integral[bin + 1], r, 1E-12)) {
         if (TMath::AreEqualRel(integral[bin + 2], r, 1E-12)) bin++;
         else break;
      }

      const Double_t rr = r - integral[bin];
      if (rr != 0.0) {
         Double_t xx = 0.0;
         const Double_t fac = -2.*gamma[bin] * rr / beta[bin] / beta[bin];
         if (fac != 0 && fac <= 1)
            xx = (-beta[bin] + TMath::Sqrt(beta[bin] * beta[bin] + 2 * gamma[bin] * rr)) / gamma[bin];
         else if (beta[bin] != 0.)
            xx = rr / beta[bin];
         q[i] = alpha[bin] + xx;
      } else {
         q[i] = alpha[bin];
         if (integral[bin + 1] == r) q[i] += dx;
      }
   }

   return nprobSum;
}
////////////////////////////////////////////////////////////////////////////////
///
///  Compute the cumulative function at fNpx points between fXmin and fXmax.
///  Option can be used to force a log scale (option = "log"), linear (option = "lin") or automatic if empty.
Bool_t TF1::ComputeCdfTable(Option_t * option) {

   fIntegral.resize(fNpx + 1);
   fAlpha.resize(fNpx + 1);
   fBeta.resize(fNpx);
   fGamma.resize(fNpx);
   fIntegral[0] = 0;
   fAlpha[fNpx] = 0;
   Double_t integ;
   Int_t intNegative = 0;
   Int_t i;
   Bool_t logbin = kFALSE;
   Double_t dx;
   Double_t xmin = fXmin;
   Double_t xmax = fXmax;
   TString opt(option);
   opt.ToUpper();
   // perform a log binning if specified by user (option="Log") or if some conditions are met
   // and the user explicitly does not specify a Linear binning option
   if (opt.Contains("LOG") || ((xmin > 0 && xmax / xmin > fNpx) && !opt.Contains("LIN"))) {
      logbin = kTRUE;
      fAlpha[fNpx] = 1;
      xmin = TMath::Log10(fXmin);
      xmax = TMath::Log10(fXmax);
      if (gDebug)
         Info("GetRandom", "Use log scale for tabulating the integral in [%f,%f] with %d points", fXmin, fXmax, fNpx);
   }
   dx = (xmax - xmin) / fNpx;

   std::vector<Double_t> xx(fNpx + 1);
   for (i = 0; i < fNpx; i++) {
      xx[i] = xmin + i * dx;
   }
   xx[fNpx] = xmax;
   for (i = 0; i < fNpx; i++) {
      if (logbin) {
         integ = Integral(TMath::Power(10, xx[i]), TMath::Power(10, xx[i + 1]), 0.0);
      } else {
         integ = Integral(xx[i], xx[i + 1], 0.0);
      }
      if (integ < 0) {
         intNegative++;
         integ = -integ;
      }
      fIntegral[i + 1] = fIntegral[i] + integ;
   }
   if (intNegative > 0) {
      Warning("GetRandom", "function:%s has %d negative values: abs assumed", GetName(), intNegative);
   }
   if (fIntegral[fNpx] == 0) {
      Error("GetRandom", "Integral of function is zero");
      return kFALSE;
   }
   Double_t total = fIntegral[fNpx];
   for (i = 1; i <= fNpx; i++) { // normalize integral to 1
      fIntegral[i] /= total;
   }
   // the integral r for each bin is approximated by a parabola
   //  x = alpha + beta*r +gamma*r**2
   // compute the coefficients alpha, beta, gamma for each bin
   Double_t x0, r1, r2, r3;
   for (i = 0; i < fNpx; i++) {
      x0 = xx[i];
      r2 = fIntegral[i + 1] - fIntegral[i];
      if (logbin)
         r1 = Integral(TMath::Power(10, x0), TMath::Power(10, x0 + 0.5 * dx), 0.0) / total;
      else
         r1 = Integral(x0, x0 + 0.5 * dx, 0.0) / total;
      r3 = 2 * r2 - 4 * r1;
      if (TMath::Abs(r3) > 1e-8)
         fGamma[i] = r3 / (dx * dx);
      else
         fGamma[i] = 0;
      fBeta[i] = r2 / dx - fGamma[i] * dx;
      fAlpha[i] = x0;
      fGamma[i] *= 2;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a random number following this function shape.
///
/// @param rng  Random number generator. By default (or when passing a nullptr) the global gRandom is used
/// @param option Option string which controls the binning used to compute the integral. Default mode is automatic depending of
///               xmax, xmin and Npx (function points).
///               Possible values are:
///               -  "LOG" to force usage of log scale for tabulating the integral
///               -  "LIN" to force usage of linear scale when tabulating the integral
///
/// The distribution contained in the function fname (TF1) is integrated
/// over the channel contents.
/// It is normalized to 1.
/// For each bin the integral is approximated by a parabola.
/// The parabola coefficients are stored as non persistent data members
/// Getting one random number implies:
///  - Generating a random number between 0 and 1 (say r1)
///  - Look in which bin in the normalized integral r1 corresponds to
///  - Evaluate the parabolic curve in the selected bin to find the corresponding X value.
///
/// The user can provide as optional parameter a Random number generator.
/// By default gRandom is used
///
/// If the ratio fXmax/fXmin > fNpx the integral is tabulated in log scale in x
/// A log scale for the intergral is also always used if a user specifies the "LOG" option
/// Instead if a user requestes a "LIN" option the integral binning is never done in log scale
/// whatever the fXmax/fXmin ratio is
///
/// Note that the parabolic approximation is very good as soon as the number of bins is greater than 50.


Double_t TF1::GetRandom(TRandom * rng, Option_t * option)
{
   //  Check if integral array must be build
   if (fIntegral.size() == 0) {
      Bool_t ret = ComputeCdfTable(option);
      if (!ret) return TMath::QuietNaN();
   }


   // return random number
   Double_t r  = (rng) ? rng->Rndm() : gRandom->Rndm();
   Int_t bin  = TMath::BinarySearch(fNpx, fIntegral.data(), r);
   Double_t rr = r - fIntegral[bin];

   Double_t yy;
   if (fGamma[bin] != 0)
      yy = (-fBeta[bin] + TMath::Sqrt(fBeta[bin] * fBeta[bin] + 2 * fGamma[bin] * rr)) / fGamma[bin];
   else
      yy = rr / fBeta[bin];
   Double_t x = fAlpha[bin] + yy;
   if (fAlpha[fNpx] > 0) return TMath::Power(10, x);
   return x;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a random number following this function shape in [xmin,xmax]
///
///   The distribution contained in the function fname (TF1) is integrated
///   over the channel contents.
///   It is normalized to 1.
///   For each bin the integral is approximated by a parabola.
///   The parabola coefficients are stored as non persistent data members
///   Getting one random number implies:
///     - Generating a random number between 0 and 1 (say r1)
///     - Look in which bin in the normalized integral r1 corresponds to
///     - Evaluate the parabolic curve in the selected bin to find
///       the corresponding X value.
///
///   The parabolic approximation is very good as soon as the number
///   of bins is greater than 50.
///
///  @param  xmin    minimum value for  generated random numbers
///  @param  xmax    maximum value for  generated random numbers
///  @param  rng     (optional) random number generator pointer
///  @param  option  (optional) : `LOG` or `LIN` to force the usage of a log or linear scale for computing the cumulative integral table
///
///  IMPORTANT NOTE
///
///  The integral of the function is computed at fNpx points. If the function
///  has sharp peaks, you should increase the number of points (SetNpx)
///  such that the peak is correctly tabulated at several points.

Double_t TF1::GetRandom(Double_t xmin, Double_t xmax, TRandom * rng, Option_t * option)
{
   //  Check if integral array must be build
   if (fIntegral.size() == 0) {
      Bool_t ret = ComputeCdfTable(option);
      if (!ret) return TMath::QuietNaN();
   }

   // return random number
   Double_t dx   = (fXmax - fXmin) / fNpx;
   Int_t nbinmin = (Int_t)((xmin - fXmin) / dx);
   Int_t nbinmax = (Int_t)((xmax - fXmin) / dx) + 2;
   if (nbinmax > fNpx) nbinmax = fNpx;

   Double_t pmin = fIntegral[nbinmin];
   Double_t pmax = fIntegral[nbinmax];

   Double_t r, x, xx, rr;
   do {
      r  = (rng) ? rng->Uniform(pmin, pmax) : gRandom->Uniform(pmin, pmax);

      Int_t bin  = TMath::BinarySearch(fNpx, fIntegral.data(), r);
      rr = r - fIntegral[bin];

      if (fGamma[bin] != 0)
         xx = (-fBeta[bin] + TMath::Sqrt(fBeta[bin] * fBeta[bin] + 2 * fGamma[bin] * rr)) / fGamma[bin];
      else
         xx = rr / fBeta[bin];
      x = fAlpha[bin] + xx;
   } while (x < xmin || x > xmax);
   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Return range of a generic N-D function.

void TF1::GetRange(Double_t *rmin, Double_t *rmax) const
{
   int ndim = GetNdim();

   double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
   GetRange(xmin, ymin, zmin, xmax, ymax, zmax);
   for (int i = 0; i < ndim; ++i) {
      if (i == 0) {
         rmin[0] = xmin;
         rmax[0] = xmax;
      } else if (i == 1) {
         rmin[1] = ymin;
         rmax[1] = ymax;
      } else if (i == 2) {
         rmin[2] = zmin;
         rmax[2] = zmax;
      } else {
         rmin[i] = 0;
         rmax[i] = 0;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return range of a 1-D function.

void TF1::GetRange(Double_t &xmin, Double_t &xmax) const
{
   xmin = fXmin;
   xmax = fXmax;
}


////////////////////////////////////////////////////////////////////////////////
/// Return range of a 2-D function.

void TF1::GetRange(Double_t &xmin, Double_t &ymin,  Double_t &xmax, Double_t &ymax) const
{
   xmin = fXmin;
   xmax = fXmax;
   ymin = 0;
   ymax = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Return range of function.

void TF1::GetRange(Double_t &xmin, Double_t &ymin, Double_t &zmin, Double_t &xmax, Double_t &ymax, Double_t &zmax) const
{
   xmin = fXmin;
   xmax = fXmax;
   ymin = 0;
   ymax = 0;
   zmin = 0;
   zmax = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Get value corresponding to X in array of fSave values

Double_t TF1::GetSave(const Double_t *xx)
{
   if (fSave.size() == 0) return 0;
   //if (fSave == 0) return 0;
   int fNsave = fSave.size();
   Double_t x    = Double_t(xx[0]);
   Double_t y, dx, xmin, xmax, xlow, xup, ylow, yup;
   if (fParent && fParent->InheritsFrom(TH1::Class())) {
      //if parent is a histogram the function had been saved at the center of the bins
      //we make a linear interpolation between the saved values
      xmin = fSave[fNsave - 3];
      xmax = fSave[fNsave - 2];
      if (fSave[fNsave - 1] == xmax) {
         TH1 *h = (TH1 *)fParent;
         TAxis *xaxis = h->GetXaxis();
         Int_t bin1  = xaxis->FindBin(xmin);
         Int_t binup = xaxis->FindBin(xmax);
         Int_t bin   = xaxis->FindBin(x);
         if (bin < binup) {
            xlow = xaxis->GetBinCenter(bin);
            xup  = xaxis->GetBinCenter(bin + 1);
            ylow = fSave[bin - bin1];
            yup  = fSave[bin - bin1 + 1];
         } else {
            xlow = xaxis->GetBinCenter(bin - 1);
            xup  = xaxis->GetBinCenter(bin);
            ylow = fSave[bin - bin1 - 1];
            yup  = fSave[bin - bin1];
         }
         dx = xup - xlow;
         y  = ((xup * ylow - xlow * yup) + x * (yup - ylow)) / dx;
         return y;
      }
   }
   Int_t np = fNsave - 3;
   xmin = Double_t(fSave[np + 1]);
   xmax = Double_t(fSave[np + 2]);
   dx   = (xmax - xmin) / np;
   if (x < xmin || x > xmax) return 0;
   // return a Nan in case of x=nan, otherwise will crash later
   if (TMath::IsNaN(x)) return x;
   if (dx <= 0) return 0;

   Int_t bin     = Int_t((x - xmin) / dx);
   xlow = xmin + bin * dx;
   xup  = xlow + dx;
   ylow = fSave[bin];
   yup  = fSave[bin + 1];
   y    = ((xup * ylow - xlow * yup) + x * (yup - ylow)) / dx;
   return y;
}


////////////////////////////////////////////////////////////////////////////////
/// Get x axis of the function.

TAxis *TF1::GetXaxis() const
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetXaxis();
}


////////////////////////////////////////////////////////////////////////////////
/// Get y axis of the function.

TAxis *TF1::GetYaxis() const
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetYaxis();
}


////////////////////////////////////////////////////////////////////////////////
/// Get z axis of the function. (In case this object is a TF2 or TF3)

TAxis *TF1::GetZaxis() const
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetZaxis();
}



////////////////////////////////////////////////////////////////////////////////
/// Compute the gradient (derivative) wrt a parameter ipar
///
/// \param ipar index of parameter for which the derivative is computed
/// \param x point, where the derivative is computed
/// \param eps - if the errors of parameters have been computed, the step used in
/// numerical differentiation is eps*parameter_error.
///
/// if the errors have not been computed, step=eps is used
/// default value of eps = 0.01
/// Method is the same as in Derivative() function
///
/// If a parameter is fixed, the gradient on this parameter = 0

Double_t TF1::GradientPar(Int_t ipar, const Double_t *x, Double_t eps)
{
   return GradientParTempl<Double_t>(ipar, x, eps);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gradient wrt parameters
/// If the TF1 object is based on a formula expression (TFormula)
/// and TFormula::GenerateGradientPar() has been successfully called
/// automatic differentiation using CLAD is used instead of the default
/// numerical differentiation
///
/// \param x  point, were the gradient is computed
/// \param grad  used to return the computed gradient, assumed to be of at least fNpar size
/// \param eps if the errors of parameters have been computed, the step used in
/// numerical differentiation is eps*parameter_error.
///
/// if the errors have not been computed, step=eps is used
/// default value of eps = 0.01
/// Method is the same as in Derivative() function
///
/// If a parameter is fixed, the gradient on this parameter = 0

void TF1::GradientPar(const Double_t *x, Double_t *grad, Double_t eps)
{
   if (fFormula && fFormula->HasGeneratedGradient())
      fFormula->GradientPar(x,grad);
   else
      GradientParTempl<Double_t>(x, grad, eps);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize parameters addresses.

void TF1::InitArgs(const Double_t *x, const Double_t *params)
{
   if (fMethodCall) {
      Long_t args[2];
      args[0] = (Long_t)x;
      if (params) args[1] = (Long_t)params;
      else        args[1] = (Long_t)GetParameters();
      fMethodCall->SetParamPtrs(args);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Create the basic function objects

void TF1::InitStandardFunctions()
{
   TF1 *f1;
   R__LOCKGUARD(gROOTMutex);
   if (!gROOT->GetListOfFunctions()->FindObject("gaus")) {
      f1 = new TF1("gaus", "gaus", -1, 1);
      f1->SetParameters(1, 0, 1);
      f1 = new TF1("gausn", "gausn", -1, 1);
      f1->SetParameters(1, 0, 1);
      f1 = new TF1("landau", "landau", -1, 1);
      f1->SetParameters(1, 0, 1);
      f1 = new TF1("landaun", "landaun", -1, 1);
      f1->SetParameters(1, 0, 1);
      f1 = new TF1("expo", "expo", -1, 1);
      f1->SetParameters(1, 1);
      for (Int_t i = 0; i < 10; i++) {
         f1 = new TF1(Form("pol%d", i), Form("pol%d", i), -1, 1);
         f1->SetParameters(1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
         // create also chebyshev polynomial
         // (note polynomial object will not be deleted)
         // note that these functions cannot be stored
         ROOT::Math::ChebyshevPol *pol = new ROOT::Math::ChebyshevPol(i);
         Double_t min = -1;
         Double_t max = 1;
         f1 = new TF1(TString::Format("chebyshev%d", i), pol, min, max, i + 1, 1);
         f1->SetParameters(1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
      }

   }
}
////////////////////////////////////////////////////////////////////////////////
/// IntegralOneDim or analytical integral

Double_t TF1::Integral(Double_t a, Double_t b,  Double_t epsrel)
{
   Double_t error = 0;
   if (GetNumber() > 0) {
      Double_t result = 0.;
      if (gDebug) {
         Info("computing analytical integral for function %s with number %d", GetName(), GetNumber());
      }
      result = AnalyticalIntegral(this, a, b);
      // if it is a formula that havent been implemented in analytical integral a NaN is return
      if (!TMath::IsNaN(result)) return result;
      if (gDebug)
         Warning("analytical integral not available for %s - with number %d  compute numerical integral", GetName(), GetNumber());
   }
   return IntegralOneDim(a, b, epsrel, epsrel, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Return Integral of function between a and b using the given parameter values and
/// relative and absolute tolerance.
///
/// The default integrator defined in ROOT::Math::IntegratorOneDimOptions::DefaultIntegrator() is used
/// If ROOT contains the MathMore library the default integrator is set to be
/// the adaptive ROOT::Math::GSLIntegrator (based on QUADPACK) or otherwise the
/// ROOT::Math::GaussIntegrator is used
/// See the reference documentation of these classes for more information about the
/// integration algorithms
/// To change integration algorithm just do :
/// ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator(IntegratorName);
/// Valid integrator names are:
///   - Gauss  :               for ROOT::Math::GaussIntegrator
///   - GaussLegendre    :     for ROOT::Math::GaussLegendreIntegrator
///   - Adaptive         :     for ROOT::Math::GSLIntegrator adaptive method (QAG)
///   - AdaptiveSingular :     for ROOT::Math::GSLIntegrator adaptive singular method (QAGS)
///   - NonAdaptive      :     for ROOT::Math::GSLIntegrator non adaptive (QNG)
///
/// In order to use the GSL integrators one needs to have the MathMore library installed
///
/// Note 1:
///
///   Values of the function f(x) at the interval end-points A and B are not
///   required. The subprogram may therefore be used when these values are
///   undefined.
///
/// Note 2:
///
///   Instead of TF1::Integral, you may want to use the combination of
///   TF1::CalcGaussLegendreSamplingPoints and TF1::IntegralFast.
///   See an example with the following script:
///
/// ~~~ {.cpp}
///   void gint() {
///      TF1 *g = new TF1("g","gaus",-5,5);
///      g->SetParameters(1,0,1);
///      //default gaus integration method uses 6 points
///      //not suitable to integrate on a large domain
///      double r1 = g->Integral(0,5);
///      double r2 = g->Integral(0,1000);
///
///      //try with user directives computing more points
///      Int_t np = 1000;
///      double *x=new double[np];
///      double *w=new double[np];
///      g->CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
///      double r3 = g->IntegralFast(np,x,w,0,5);
///      double r4 = g->IntegralFast(np,x,w,0,1000);
///      double r5 = g->IntegralFast(np,x,w,0,10000);
///      double r6 = g->IntegralFast(np,x,w,0,100000);
///      printf("g->Integral(0,5)               = %g\n",r1);
///      printf("g->Integral(0,1000)            = %g\n",r2);
///      printf("g->IntegralFast(n,x,w,0,5)     = %g\n",r3);
///      printf("g->IntegralFast(n,x,w,0,1000)  = %g\n",r4);
///      printf("g->IntegralFast(n,x,w,0,10000) = %g\n",r5);
///      printf("g->IntegralFast(n,x,w,0,100000)= %g\n",r6);
///      delete [] x;
///      delete [] w;
///   }
/// ~~~
///
///   This example produces the following results:
///
/// ~~~ {.cpp}
///      g->Integral(0,5)               = 1.25331
///      g->Integral(0,1000)            = 1.25319
///      g->IntegralFast(n,x,w,0,5)     = 1.25331
///      g->IntegralFast(n,x,w,0,1000)  = 1.25331
///      g->IntegralFast(n,x,w,0,10000) = 1.25331
///      g->IntegralFast(n,x,w,0,100000)= 1.253
/// ~~~

Double_t TF1::IntegralOneDim(Double_t a, Double_t b,  Double_t epsrel, Double_t epsabs, Double_t &error)
{
   //Double_t *parameters = GetParameters();
   TF1_EvalWrapper wf1(this, 0, fgAbsValue);
   Double_t result = 0;
   Int_t status = 0;
   if (epsrel <= 0) epsrel = ROOT::Math::IntegratorOneDimOptions::DefaultRelTolerance();
   if (epsabs <= 0) epsabs = ROOT::Math::IntegratorOneDimOptions::DefaultAbsTolerance();
   if (ROOT::Math::IntegratorOneDimOptions::DefaultIntegratorType() == ROOT::Math::IntegrationOneDim::kGAUSS) {
      ROOT::Math::GaussIntegrator iod(epsabs, epsrel);
      iod.SetFunction(wf1);
      if (a != - TMath::Infinity() && b != TMath::Infinity())
         result =  iod.Integral(a, b);
      else if (a == - TMath::Infinity() && b != TMath::Infinity())
         result = iod.IntegralLow(b);
      else if (a != - TMath::Infinity() && b == TMath::Infinity())
         result = iod.IntegralUp(a);
      else if (a == - TMath::Infinity() && b == TMath::Infinity())
         result = iod.Integral();
      error = iod.Error();
      status = iod.Status();
   } else {
      ROOT::Math::IntegratorOneDim iod(wf1, ROOT::Math::IntegratorOneDimOptions::DefaultIntegratorType(), epsabs, epsrel);
      if (a != - TMath::Infinity() && b != TMath::Infinity())
         result =  iod.Integral(a, b);
      else if (a == - TMath::Infinity() && b != TMath::Infinity())
         result = iod.IntegralLow(b);
      else if (a != - TMath::Infinity() && b == TMath::Infinity())
         result = iod.IntegralUp(a);
      else if (a == - TMath::Infinity() && b == TMath::Infinity())
         result = iod.Integral();
      error = iod.Error();
      status = iod.Status();
   }
   if (status != 0) {
      std::string igName = ROOT::Math::IntegratorOneDim::GetName(ROOT::Math::IntegratorOneDimOptions::DefaultIntegratorType());
      Warning("IntegralOneDim", "Error found in integrating function %s in [%f,%f] using %s. Result = %f +/- %f  - status = %d", GetName(), a, b, igName.c_str(), result, error, status);
      TString msg("\t\tFunction Parameters = {");
      for (int ipar = 0; ipar < GetNpar(); ++ipar) {
         msg += TString::Format(" %s =  %f ", GetParName(ipar), GetParameter(ipar));
         if (ipar < GetNpar() - 1) msg += TString(",");
         else msg += TString("}");
      }
      Info("IntegralOneDim", "%s", msg.Data());
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Error on Integral of a parametric function between a and b
/// due to the parameter uncertainties and their covariance matrix from the fit.
/// In addition to the integral limits, this method takes as input a pointer to the fitted parameter values
/// and a pointer the covariance matrix from the fit. These pointers should be retrieved from the
/// previously performed fit using the TFitResult class.
/// Note that to get the TFitResult, te fit should be done using the fit option `S`.
/// Example:
/// ~~~~{.cpp}
/// TFitResultPtr r = histo->Fit(func, "S");
/// func->IntegralError(x1,x2,r->GetParams(), r->GetCovarianceMatrix()->GetMatrixArray() );
/// ~~~~
///
/// IMPORTANT NOTE1:
///
/// A null pointer to the parameter values vector and to the covariance matrix can be passed.
/// In this case, when the parameter values pointer is null, the parameter values stored in this
/// TF1 function object are used in the integral error computation.
/// When the poassed pointer to the covariance matrix is null, a covariance matrix from the last fit is retrieved
/// from a global fitter instance when it exists. Note that the global fitter instance
/// esists only when ROOT is not running with multi-threading enabled (ROOT::IsImplicitMTEnabled() == True).
/// When the ovariance matrix from the last fit cannot be retrieved, an error message is printed and a a zero value is
/// returned.
///
///
/// IMPORTANT NOTE2:
///
/// When no covariance matrix is passed and in the meantime a fit is done
/// using another function, the routine will signal an error and it will return zero only
/// when the number of fit parameter is different than the values stored in TF1 (TF1::GetNpar() ).
/// In the case that npar is the same, an incorrect result is returned.
///
/// IMPORTANT NOTE3:
///
/// The user must pass a pointer to the elements of the full covariance matrix
/// dimensioned with the right size (npar*npar), where npar is the total number of parameters (TF1::GetNpar()),
/// including also the fixed parameters. The covariance matrix must be retrieved from the TFitResult class as
/// shown above and not from TVirtualFitter::GetCovarianceMatrix() function.

Double_t TF1::IntegralError(Double_t a, Double_t b, const Double_t *params, const Double_t *covmat, Double_t epsilon)
{
   Double_t x1[1];
   Double_t x2[1];
   x1[0] = a, x2[0] = b;
   return ROOT::TF1Helper::IntegralError(this, 1, x1, x2, params, covmat, epsilon);
}

////////////////////////////////////////////////////////////////////////////////
/// Return Error on Integral of a parametric function with dimension larger than one
/// between a[] and b[]  due to the parameters uncertainties.
/// For a TF1 with dimension larger than 1 (for example a TF2 or TF3)
/// TF1::IntegralMultiple is used for the integral calculation
///
/// In addition to the integral limits, this method takes as input a pointer to the fitted parameter values
/// and a pointer the covariance matrix from the fit. These pointers should be retrieved from the
/// previously performed fit using the TFitResult class.
/// Note that to get the TFitResult, te fit should be done using the fit option `S`.
/// Example:
/// ~~~~{.cpp}
/// TFitResultPtr r = histo2d->Fit(func2, "S");
/// func2->IntegralError(a,b,r->GetParams(), r->GetCovarianceMatrix()->GetMatrixArray() );
/// ~~~~
///
/// IMPORTANT NOTE1:
///
/// A null pointer to the parameter values vector and to the covariance matrix can be passed.
/// In this case, when the parameter values pointer is null, the parameter values stored in this
/// TF1 function object are used in the integral error computation.
/// When the poassed pointer to the covariance matrix is null, a covariance matrix from the last fit is retrieved
/// from a global fitter instance when it exists. Note that the global fitter instance
/// esists only when ROOT is not running with multi-threading enabled (ROOT::IsImplicitMTEnabled() == True).
/// When the ovariance matrix from the last fit cannot be retrieved, an error message is printed and a a zero value is
/// returned.
///
///
/// IMPORTANT NOTE2:
///
/// When no covariance matrix is passed and in the meantime a fit is done
/// using another function, the routine will signal an error and it will return zero only
/// when the number of fit parameter is different than the values stored in TF1 (TF1::GetNpar() ).
/// In the case that npar is the same, an incorrect result is returned.
///
/// IMPORTANT NOTE3:
///
/// The user must pass a pointer to the elements of the full covariance matrix
/// dimensioned with the right size (npar*npar), where npar is the total number of parameters (TF1::GetNpar()),
/// including also the fixed parameters. The covariance matrix must be retrieved from the TFitResult class as
/// shown above and not from TVirtualFitter::GetCovarianceMatrix() function.

Double_t TF1::IntegralError(Int_t n, const Double_t *a, const Double_t *b, const Double_t *params, const  Double_t *covmat, Double_t epsilon)
{
   return ROOT::TF1Helper::IntegralError(this, n, a, b, params, covmat, epsilon);
}

#ifdef INTHEFUTURE
////////////////////////////////////////////////////////////////////////////////
/// Gauss-Legendre integral, see CalcGaussLegendreSamplingPoints

Double_t TF1::IntegralFast(const TGraph *g, Double_t a, Double_t b, Double_t *params)
{
   if (!g) return 0;
   return IntegralFast(g->GetN(), g->GetX(), g->GetY(), a, b, params);
}
#endif


////////////////////////////////////////////////////////////////////////////////
/// Gauss-Legendre integral, see CalcGaussLegendreSamplingPoints

Double_t TF1::IntegralFast(Int_t num, Double_t * /* x */, Double_t * /* w */, Double_t a, Double_t b, Double_t *params, Double_t epsilon)
{
   // Now x and w are not used!

   ROOT::Math::WrappedTF1 wf1(*this);
   if (params)
      wf1.SetParameters(params);
   ROOT::Math::GaussLegendreIntegrator gli(num, epsilon);
   gli.SetFunction(wf1);
   return gli.Integral(a, b);

}


////////////////////////////////////////////////////////////////////////////////
/// See more general prototype below.
/// This interface kept for back compatibility
/// It is recommended to use the other interface where one can specify also epsabs and the maximum number of
/// points

Double_t TF1::IntegralMultiple(Int_t n, const Double_t *a, const Double_t *b, Double_t epsrel, Double_t &relerr)
{
   Int_t nfnevl, ifail;
   UInt_t maxpts = TMath::Max(UInt_t(20 * TMath::Power(fNpx, GetNdim())), ROOT::Math::IntegratorMultiDimOptions::DefaultNCalls());
   Double_t result = IntegralMultiple(n, a, b, maxpts, epsrel, epsrel, relerr, nfnevl, ifail);
   if (ifail > 0) {
      Warning("IntegralMultiple", "failed code=%d, ", ifail);
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// This function computes, to an attempted specified accuracy, the value of
/// the integral
///
/// \param[in] n   Number of dimensions [2,15]
/// \param[in] a,b One-dimensional arrays of length >= N . On entry A[i],  and  B[i],
///   contain the lower and upper limits of integration, respectively.
/// \param[in] maxpts Maximum number of function evaluations to be allowed.
///   maxpts >= 2^n +2*n*(n+1) +1
///   if maxpts<minpts, maxpts is set to 10*minpts
/// \param[in] epsrel Specified relative accuracy.
/// \param[in] epsabs Specified absolute accuracy.
///   The integration algorithm will attempt to reach either the relative or the absolute accuracy.
///   In case the maximum function called is reached the algorithm will stop earlier without having reached
///   the desired accuracy
///
/// \param[out] relerr Contains, on exit, an estimation of the relative accuracy of the result.
/// \param[out] nfnevl number of function evaluations performed.
/// \param[out] ifail
/// \parblock
///        0 Normal exit. At least minpts and at most maxpts calls to the function were performed.
///
///        1 maxpts is too small for the specified accuracy eps. The result and relerr contain the values obtainable for the
///          specified value of maxpts.
///
///        3 n<2 or n>15
/// \endparblock
///
/// Method:
///
/// The default method used is the Genz-Mallik adaptive multidimensional algorithm
/// using the class ROOT::Math::AdaptiveIntegratorMultiDim (see the reference documentation of the class)
///
/// Other methods can be used by setting ROOT::Math::IntegratorMultiDimOptions::SetDefaultIntegrator()
/// to different integrators.
/// Other possible integrators are MC integrators based on the ROOT::Math::GSLMCIntegrator class
/// Possible methods are : Vegas, Miser or Plain
/// IN case of MC integration the accuracy is determined by the number of function calls, one should be
/// careful not to use a too large value of maxpts
///

Double_t TF1::IntegralMultiple(Int_t n, const Double_t *a, const Double_t *b, Int_t maxpts, Double_t epsrel, Double_t epsabs, Double_t &relerr, Int_t &nfnevl, Int_t &ifail)
{
   ROOT::Math::WrappedMultiFunction<TF1 &> wf1(*this, n);

   double result = 0;
   if (epsrel <= 0) epsrel = ROOT::Math::IntegratorMultiDimOptions::DefaultRelTolerance();
   if (epsabs <= 0) epsabs = ROOT::Math::IntegratorMultiDimOptions::DefaultAbsTolerance();
   if (ROOT::Math::IntegratorMultiDimOptions::DefaultIntegratorType() == ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
      ROOT::Math::AdaptiveIntegratorMultiDim aimd(wf1, epsabs, epsrel, maxpts);
      //aimd.SetMinPts(minpts); // use default minpts ( n^2 + 2 * n * (n+1) +1 )
      result = aimd.Integral(a, b);
      relerr = aimd.RelError();
      nfnevl = aimd.NEval();
      ifail =  aimd.Status();
   } else {
      // use default abs tolerance = relative tolerance
      ROOT::Math::IntegratorMultiDim imd(wf1, ROOT::Math::IntegratorMultiDimOptions::DefaultIntegratorType(), epsabs, epsrel, maxpts);
      result = imd.Integral(a, b);
      relerr = (result != 0) ? imd.Error() / std::abs(result) : imd.Error();
      nfnevl = 0;
      ifail = imd.Status();
   }


   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the function is valid

Bool_t TF1::IsValid() const
{
   if (fFormula) return fFormula->IsValid();
   if (fMethodCall) return fMethodCall->IsValid();
   // function built on compiled functors are always valid by definition
   // (checked at compiled time)
   // invalid is a TF1 where the functor is null pointer and has not been saved
   if (!fFunctor && fSave.empty()) return kFALSE;
   return kTRUE;
}


//______________________________________________________________________________


void TF1::Print(Option_t *option) const
{
   if (fType == EFType::kFormula) {
      printf("Formula based function:     %s \n", GetName());
      assert(fFormula);
      fFormula->Print(option);
   } else if (fType >  0) {
      if (fType == EFType::kInterpreted)
         printf("Interpreted based function: %s(double *x, double *p).  Ndim = %d, Npar = %d  \n", GetName(), GetNdim(),
                GetNpar());
      else if (fType == EFType::kCompositionFcn) {
         printf("Composition based function: %s. Ndim = %d, Npar = %d \n", GetName(), GetNdim(), GetNpar());
         if (!fComposition)
            printf("fComposition not found!\n"); // this would be bad
      } else {
         if (fFunctor)
            printf("Compiled based function: %s  based on a functor object.  Ndim = %d, Npar = %d\n", GetName(),
                   GetNdim(), GetNpar());
         else {
            printf("Function based on a list of points from a compiled based function: %s.  Ndim = %d, Npar = %d, Npx "
                   "= %zu\n",
                   GetName(), GetNdim(), GetNpar(), fSave.size());
            if (fSave.empty())
               Warning("Print", "Function %s is based on a list of points but list is empty", GetName());
         }
      }
      TString opt(option);
      opt.ToUpper();
      if (opt.Contains("V")) {
         // print list of parameters
         if (fNpar > 0) {
            printf("List of  Parameters: \n");
            for (int i = 0; i < fNpar; ++i)
               printf(" %20s =  %10f \n", GetParName(i), GetParameter(i));
         }
         if (!fSave.empty()) {
            // print list of saved points
            printf("List of  Saved points (N=%d): \n", int(fSave.size()));
            for (auto &x : fSave)
               printf("( %10f )  ", x);
            printf("\n");
         }
      }
   }
   if (fHistogram) {
      printf("Contained histogram\n");
      fHistogram->Print(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this function with its current attributes.
/// The function is going to be converted in an histogram and the corresponding
/// histogram is painted.
/// The painted histogram can be retrieved calling afterwards the method TF1::GetHistogram()

void TF1::Paint(Option_t *choptin)
{
   fgCurrent = this;

   char option[32];
   strlcpy(option,choptin,32);

   TString opt = option;
   opt.ToLower();

   Bool_t optSAME = kFALSE;
   if (opt.Contains("same")) {
      opt.ReplaceAll("same","");
      optSAME = kTRUE;
   }
   opt.ReplaceAll(' ', "");

   Double_t xmin = fXmin, xmax = fXmax, pmin = fXmin, pmax = fXmax;
   if (gPad) {
      pmin = gPad->PadtoX(gPad->GetUxmin());
      pmax = gPad->PadtoX(gPad->GetUxmax());
   }
   if (optSAME) {
      if (xmax < pmin) return;  // Completely outside.
      if (xmin > pmax) return;
      if (xmin < pmin) xmin = pmin;
      if (xmax > pmax) xmax = pmax;
   }

   // create an histogram using the function content (re-use it if already existing)
   fHistogram = DoCreateHistogram(xmin, xmax, kFALSE);

   char *l1 = strstr(option,"PFC"); // Automatic Fill Color
   char *l2 = strstr(option,"PLC"); // Automatic Line Color
   char *l3 = strstr(option,"PMC"); // Automatic Marker Color
   if (l1 || l2 || l3) {
      Int_t i = gPad->NextPaletteColor();
      if (l1) {memcpy(l1,"   ",3); fHistogram->SetFillColor(i);}
      if (l2) {memcpy(l2,"   ",3); fHistogram->SetLineColor(i);}
      if (l3) {memcpy(l3,"   ",3); fHistogram->SetMarkerColor(i);}
   }

   // set the optimal minimum and maximum
   Double_t minimum   = fHistogram->GetMinimumStored();
   Double_t maximum   = fHistogram->GetMaximumStored();
   if (minimum <= 0 && gPad && gPad->GetLogy()) minimum = -1111; // This can happen when switching from lin to log scale.
   if (gPad && gPad->GetUymin() < fHistogram->GetMinimum() &&
         !fHistogram->TestBit(TH1::kIsZoomed)) minimum = -1111; // This can happen after unzooming a fit.
   if (minimum == -1111) { // This can happen after unzooming.
      if (fHistogram->TestBit(TH1::kIsZoomed)) {
         minimum = fHistogram->GetYaxis()->GetXmin();
      } else {
         minimum = fMinimum;
         // Optimize the computation of the scale in Y in case the min/max of the
         // function oscillate around a constant value
         if (minimum == -1111) {
            Double_t hmin;
            if (optSAME && gPad) hmin = gPad->GetUymin();
            else         hmin = fHistogram->GetMinimum();
            if (hmin > 0) {
               Double_t hmax;
               Double_t hminpos = hmin;
               if (optSAME && gPad) hmax = gPad->GetUymax();
               else         hmax = fHistogram->GetMaximum();
               hmin -= 0.05 * (hmax - hmin);
               if (hmin < 0) hmin = 0;
               if (hmin <= 0 && gPad && gPad->GetLogy()) hmin = hminpos;
               minimum = hmin;
            }
         }
      }
      fHistogram->SetMinimum(minimum);
   }
   if (maximum == -1111) {
      if (fHistogram->TestBit(TH1::kIsZoomed)) {
         maximum = fHistogram->GetYaxis()->GetXmax();
      } else {
         maximum = fMaximum;
      }
      fHistogram->SetMaximum(maximum);
   }


   // Draw the histogram.
   if (!gPad) return;
   if (opt.Length() == 0) {
      if (optSAME) fHistogram->Paint("lfsame");
      else         fHistogram->Paint("lf");
   } else {
      fHistogram->Paint(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create histogram with bin content equal to function value
/// computed at the bin center
/// This histogram will be used to paint the function
/// A re-creation is forced and a new histogram is done if recreate=true

TH1   *TF1::DoCreateHistogram(Double_t xmin, Double_t  xmax, Bool_t recreate)
{
   Int_t i;
   Double_t xv[1];

   TH1 *histogram = 0;


   //  Create a temporary histogram and fill each channel with the function value
   //  Preserve axis titles
   TString xtitle = "";
   TString ytitle = "";
   char *semicol = (char *)strstr(GetTitle(), ";");
   if (semicol) {
      Int_t nxt = strlen(semicol);
      char *ctemp = new char[nxt];
      strlcpy(ctemp, semicol + 1, nxt);
      semicol = (char *)strstr(ctemp, ";");
      if (semicol) {
         *semicol = 0;
         ytitle = semicol + 1;
      }
      xtitle = ctemp;
      delete [] ctemp;
   }
   if (fHistogram) {
      // delete previous histograms if were done if done in different mode
      xtitle = fHistogram->GetXaxis()->GetTitle();
      ytitle = fHistogram->GetYaxis()->GetTitle();
      if (!gPad->GetLogx()  &&  fHistogram->TestBit(TH1::kLogX)) {
         delete fHistogram;
         fHistogram = 0;
         recreate = kTRUE;
      }
      if (gPad->GetLogx()  && !fHistogram->TestBit(TH1::kLogX)) {
         delete fHistogram;
         fHistogram = 0;
         recreate = kTRUE;
      }
   }

   if (fHistogram && !recreate) {
      histogram = fHistogram;
      fHistogram->GetXaxis()->SetLimits(xmin, xmax);
   } else {
      // If logx, we must bin in logx and not in x
      // otherwise in case of several decades, one gets wrong results.
      if (xmin > 0 && gPad && gPad->GetLogx()) {
         Double_t *xbins    = new Double_t[fNpx + 1];
         Double_t xlogmin = TMath::Log10(xmin);
         Double_t xlogmax = TMath::Log10(xmax);
         Double_t dlogx   = (xlogmax - xlogmin) / ((Double_t)fNpx);
         for (i = 0; i <= fNpx; i++) {
            xbins[i] = gPad->PadtoX(xlogmin + i * dlogx);
         }
         histogram = new TH1D("Func", GetTitle(), fNpx, xbins);
         histogram->SetBit(TH1::kLogX);
         delete [] xbins;
      } else {
         histogram = new TH1D("Func", GetTitle(), fNpx, xmin, xmax);
      }
      if (fMinimum != -1111) histogram->SetMinimum(fMinimum);
      if (fMaximum != -1111) histogram->SetMaximum(fMaximum);
      histogram->SetDirectory(0);
   }
   R__ASSERT(histogram);

   // Restore axis titles.
   histogram->GetXaxis()->SetTitle(xtitle.Data());
   histogram->GetYaxis()->SetTitle(ytitle.Data());
   Double_t *parameters = GetParameters();

   InitArgs(xv, parameters);
   for (i = 1; i <= fNpx; i++) {
      xv[0] = histogram->GetBinCenter(i);
      histogram->SetBinContent(i, EvalPar(xv, parameters));
   }

   // Copy Function attributes to histogram attributes.
   histogram->SetBit(TH1::kNoStats);
   histogram->SetLineColor(GetLineColor());
   histogram->SetLineStyle(GetLineStyle());
   histogram->SetLineWidth(GetLineWidth());
   histogram->SetFillColor(GetFillColor());
   histogram->SetFillStyle(GetFillStyle());
   histogram->SetMarkerColor(GetMarkerColor());
   histogram->SetMarkerStyle(GetMarkerStyle());
   histogram->SetMarkerSize(GetMarkerSize());

   // update saved histogram in case it was deleted or if it is the first time the method is called
   // for example when called from TF1::GetHistogram()
   if (!fHistogram) fHistogram = histogram;
   return histogram;

}


////////////////////////////////////////////////////////////////////////////////
/// Release parameter number ipar If used in a fit, the parameter
/// can vary freely. The parameter limits are reset to 0,0.

void TF1::ReleaseParameter(Int_t ipar)
{
   if (ipar < 0 || ipar > GetNpar() - 1) return;
   SetParLimits(ipar, 0, 0);
}


////////////////////////////////////////////////////////////////////////////////
/// Save values of function in array fSave

void TF1::Save(Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t)
{
   Double_t *parameters = GetParameters();
   //if (fSave != 0) {delete [] fSave; fSave = 0;}
   if (fParent && fParent->InheritsFrom(TH1::Class())) {
      //if parent is a histogram save the function at the center of the bins
      if ((xmin > 0 && xmax > 0) && TMath::Abs(TMath::Log10(xmax / xmin) > TMath::Log10(fNpx))) {
         TH1 *h = (TH1 *)fParent;
         Int_t bin1 = h->GetXaxis()->FindBin(xmin);
         Int_t bin2 = h->GetXaxis()->FindBin(xmax);
         int fNsave = bin2 - bin1 + 4;
         //fSave  = new Double_t[fNsave];
         fSave.resize(fNsave);
         Double_t xv[1];

         InitArgs(xv, parameters);
         for (Int_t i = bin1; i <= bin2; i++) {
            xv[0]    = h->GetXaxis()->GetBinCenter(i);
            fSave[i - bin1] = EvalPar(xv, parameters);
         }
         fSave[fNsave - 3] = xmin;
         fSave[fNsave - 2] = xmax;
         fSave[fNsave - 1] = xmax;
         return;
      }
   }
   int fNsave = fNpx + 3;
   if (fNsave <= 3) {
      return;
   }
   //fSave  = new Double_t[fNsave];
   fSave.resize(fNsave);
   Double_t dx = (xmax - xmin) / fNpx;
   if (dx <= 0) {
      dx = (fXmax - fXmin) / fNpx;
      fNsave--;
      xmin = fXmin + 0.5 * dx;
      xmax = fXmax - 0.5 * dx;
   }
   Double_t xv[1];
   InitArgs(xv, parameters);
   for (Int_t i = 0; i <= fNpx; i++) {
      xv[0]    = xmin + dx * i;
      fSave[i] = EvalPar(xv, parameters);
   }
   fSave[fNpx + 1] = xmin;
   fSave[fNpx + 2] = xmax;
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TF1::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   Int_t i;
   char quote = '"';

   // Save the function as C code independant from ROOT.
   if (strstr(option, "cc")) {
      out << "double " << GetName() << "(double xv) {" << std::endl;
      Double_t dx = (fXmax - fXmin) / (fNpx - 1);
      out << "   double x[" << fNpx << "] = {" << std::endl;
      out << "   ";
      Int_t n = 0;
      for (i = 0; i < fNpx; i++) {
         out << fXmin + dx *i ;
         if (i < fNpx - 1) out << ", ";
         if (n++ == 10) {
            out << std::endl;
            out << "   ";
            n = 0;
         }
      }
      out << std::endl;
      out << "   };" << std::endl;
      out << "   double y[" << fNpx << "] = {" << std::endl;
      out << "   ";
      n = 0;
      for (i = 0; i < fNpx; i++) {
         out << Eval(fXmin + dx * i);
         if (i < fNpx - 1) out << ", ";
         if (n++ == 10) {
            out << std::endl;
            out << "   ";
            n = 0;
         }
      }
      out << std::endl;
      out << "   };" << std::endl;
      out << "   if (xv<x[0]) return y[0];" << std::endl;
      out << "   if (xv>x[" << fNpx - 1 << "]) return y[" << fNpx - 1 << "];" << std::endl;
      out << "   int i, j=0;" << std::endl;
      out << "   for (i=1; i<" << fNpx << "; i++) { if (xv < x[i]) break; j++; }" << std::endl;
      out << "   return y[j] + (y[j + 1] - y[j]) / (x[j + 1] - x[j]) * (xv - x[j]);" << std::endl;
      out << "}" << std::endl;
      return;
   }

   out << "   " << std::endl;

   // Either f1Number is computed locally or set from outside
   static Int_t f1Number = 0;
   TString f1Name(GetName());
   const char *l = strstr(option, "#");
   if (l != 0) {
      sscanf(&l[1], "%d", &f1Number);
   } else {
      ++f1Number;
   }
   f1Name += f1Number;

   const char *addToGlobList = fParent ? ", TF1::EAddToList::kNo" : ", TF1::EAddToList::kDefault";

   if (!fType) {
      out << "   TF1 *" << f1Name.Data() << " = new TF1(" << quote << GetName() << quote << "," << quote << GetTitle() << quote << "," << fXmin << "," << fXmax <<  addToGlobList << ");" << std::endl;
      if (fNpx != 100) {
         out << "   " << f1Name.Data() << "->SetNpx(" << fNpx << ");" << std::endl;
      }
   } else {
      out << "   TF1 *" << f1Name.Data() << " = new TF1(" << quote << "*" << GetName() << quote << "," << fXmin << "," << fXmax << "," << GetNpar() << ");" << std::endl;
      out << "    //The original function : " << GetTitle() << " had originally been created by:" << std::endl;
      out << "    //TF1 *" << GetName() << " = new TF1(" << quote << GetName() << quote << "," << GetTitle() << "," << fXmin << "," << fXmax << "," << GetNpar();
      out << ", 1" << addToGlobList << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetRange(" << fXmin << "," << fXmax << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetName(" << quote << GetName() << quote << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;
      if (fNpx != 100) {
         out << "   " << f1Name.Data() << "->SetNpx(" << fNpx << ");" << std::endl;
      }
      Double_t dx = (fXmax - fXmin) / fNpx;
      Double_t xv[1];
      Double_t *parameters = GetParameters();
      InitArgs(xv, parameters);
      for (i = 0; i <= fNpx; i++) {
         xv[0]    = fXmin + dx * i;
         Double_t save = EvalPar(xv, parameters);
         out << "   " << f1Name.Data() << "->SetSavedPoint(" << i << "," << save << ");" << std::endl;
      }
      out << "   " << f1Name.Data() << "->SetSavedPoint(" << fNpx + 1 << "," << fXmin << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetSavedPoint(" << fNpx + 2 << "," << fXmax << ");" << std::endl;
   }

   if (TestBit(kNotDraw)) {
      out << "   " << f1Name.Data() << "->SetBit(TF1::kNotDraw);" << std::endl;
   }
   if (GetFillColor() != 0) {
      if (GetFillColor() > 228) {
         TColor::SaveColor(out, GetFillColor());
         out << "   " << f1Name.Data() << "->SetFillColor(ci);" << std::endl;
      } else
         out << "   " << f1Name.Data() << "->SetFillColor(" << GetFillColor() << ");" << std::endl;
   }
   if (GetFillStyle() != 1001) {
      out << "   " << f1Name.Data() << "->SetFillStyle(" << GetFillStyle() << ");" << std::endl;
   }
   if (GetMarkerColor() != 1) {
      if (GetMarkerColor() > 228) {
         TColor::SaveColor(out, GetMarkerColor());
         out << "   " << f1Name.Data() << "->SetMarkerColor(ci);" << std::endl;
      } else
         out << "   " << f1Name.Data() << "->SetMarkerColor(" << GetMarkerColor() << ");" << std::endl;
   }
   if (GetMarkerStyle() != 1) {
      out << "   " << f1Name.Data() << "->SetMarkerStyle(" << GetMarkerStyle() << ");" << std::endl;
   }
   if (GetMarkerSize() != 1) {
      out << "   " << f1Name.Data() << "->SetMarkerSize(" << GetMarkerSize() << ");" << std::endl;
   }
   if (GetLineColor() != 1) {
      if (GetLineColor() > 228) {
         TColor::SaveColor(out, GetLineColor());
         out << "   " << f1Name.Data() << "->SetLineColor(ci);" << std::endl;
      } else
         out << "   " << f1Name.Data() << "->SetLineColor(" << GetLineColor() << ");" << std::endl;
   }
   if (GetLineWidth() != 4) {
      out << "   " << f1Name.Data() << "->SetLineWidth(" << GetLineWidth() << ");" << std::endl;
   }
   if (GetLineStyle() != 1) {
      out << "   " << f1Name.Data() << "->SetLineStyle(" << GetLineStyle() << ");" << std::endl;
   }
   if (GetChisquare() != 0) {
      out << "   " << f1Name.Data() << "->SetChisquare(" << GetChisquare() << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetNDF(" << GetNDF() << ");" << std::endl;
   }

   if (GetXaxis()) GetXaxis()->SaveAttributes(out, f1Name.Data(), "->GetXaxis()");
   if (GetYaxis()) GetYaxis()->SaveAttributes(out, f1Name.Data(), "->GetYaxis()");

   Double_t parmin, parmax;
   for (i = 0; i < GetNpar(); i++) {
      out << "   " << f1Name.Data() << "->SetParameter(" << i << "," << GetParameter(i) << ");" << std::endl;
      out << "   " << f1Name.Data() << "->SetParError(" << i << "," << GetParError(i) << ");" << std::endl;
      GetParLimits(i, parmin, parmax);
      out << "   " << f1Name.Data() << "->SetParLimits(" << i << "," << parmin << "," << parmax << ");" << std::endl;
   }
   if (!strstr(option, "nodraw")) {
      out << "   " << f1Name.Data() << "->Draw("
          << quote << option << quote << ");" << std::endl;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Static function setting the current function.
/// the current function may be accessed in static C-like functions
/// when fitting or painting a function.

void TF1::SetCurrent(TF1 *f1)
{
   fgCurrent = f1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the result from the fit
/// parameter values, errors, chi2, etc...
/// Optionally a pointer to a vector (with size fNpar) of the parameter indices in the FitResult can be passed
/// This is useful in the case of a combined fit with different functions, and the FitResult contains the global result
/// By default it is assume that indpar = {0,1,2,....,fNpar-1}.

void TF1::SetFitResult(const ROOT::Fit::FitResult &result, const Int_t *indpar)
{
   Int_t npar = GetNpar();
   if (result.IsEmpty()) {
      Warning("SetFitResult", "Empty Fit result - nothing is set in TF1");
      return;
   }
   if (indpar == 0 && npar != (int) result.NPar()) {
      Error("SetFitResult", "Invalid Fit result passed - number of parameter is %d , different than TF1::GetNpar() = %d", npar, result.NPar());
      return;
   }
   if (result.Chi2() > 0)
      SetChisquare(result.Chi2());
   else
      SetChisquare(result.MinFcnValue());

   SetNDF(result.Ndf());
   SetNumberFitPoints(result.Ndf() + result.NFreeParameters());


   for (Int_t i = 0; i < npar; ++i) {
      Int_t ipar = (indpar != 0) ? indpar[i] : i;
      if (ipar < 0) continue;
      GetParameters()[i] = result.Parameter(ipar);
      // in case errors are not present do not set them
      if (ipar < (int) result.Errors().size())
         fParErrors[i] = result.Error(ipar);
   }
   //invalidate cached integral since parameters have changed
   Update();

}


////////////////////////////////////////////////////////////////////////////////
/// Set the maximum value along Y for this function
/// In case the function is already drawn, set also the maximum in the
/// helper histogram

void TF1::SetMaximum(Double_t maximum)
{
   fMaximum = maximum;
   if (fHistogram) fHistogram->SetMaximum(maximum);
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Set the minimum value along Y for this function
/// In case the function is already drawn, set also the minimum in the
/// helper histogram

void TF1::SetMinimum(Double_t minimum)
{
   fMinimum = minimum;
   if (fHistogram) fHistogram->SetMinimum(minimum);
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Set the number of degrees of freedom
/// ndf should be the number of points used in a fit - the number of free parameters

void TF1::SetNDF(Int_t ndf)
{
   fNDF = ndf;
}


////////////////////////////////////////////////////////////////////////////////
/// Set the number of points used to draw the function
///
/// The default number of points along x is 100 for 1-d functions and 30 for 2-d/3-d functions
/// You can increase this value to get a better resolution when drawing
/// pictures with sharp peaks or to get a better result when using TF1::GetRandom
/// the minimum number of points is 4, the maximum is 10000000 for 1-d and 10000 for 2-d/3-d functions

void TF1::SetNpx(Int_t npx)
{
   const Int_t minPx = 4;
   Int_t maxPx = 10000000;
   if (GetNdim() > 1) maxPx = 10000;
   if (npx >= minPx && npx <= maxPx) {
      fNpx = npx;
   } else {
      if (npx < minPx) fNpx = minPx;
      if (npx > maxPx) fNpx = maxPx;
      Warning("SetNpx", "Number of points must be >=%d && <= %d, fNpx set to %d", minPx, maxPx, fNpx);
   }
   Update();
}
////////////////////////////////////////////////////////////////////////////////
/// Set name of parameter number ipar

void TF1::SetParName(Int_t ipar, const char *name)
{
   if (fFormula) {
      if (ipar < 0 || ipar >= GetNpar()) return;
      fFormula->SetParName(ipar, name);
   } else
      fParams->SetParName(ipar, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set up to 10 parameter names

void TF1::SetParNames(const char *name0, const char *name1, const char *name2, const char *name3, const char *name4,
                      const char *name5, const char *name6, const char *name7, const char *name8, const char *name9, const char *name10)
{
   if (fFormula)
      fFormula->SetParNames(name0, name1, name2, name3, name4, name5, name6, name7, name8, name9, name10);
   else
      fParams->SetParNames(name0, name1, name2, name3, name4, name5, name6, name7, name8, name9, name10);
}
////////////////////////////////////////////////////////////////////////////////
/// Set error for parameter number ipar

void TF1::SetParError(Int_t ipar, Double_t error)
{
   if (ipar < 0 || ipar > GetNpar() - 1) return;
   fParErrors[ipar] = error;
}


////////////////////////////////////////////////////////////////////////////////
/// Set errors for all active parameters
/// when calling this function, the array errors must have at least fNpar values

void TF1::SetParErrors(const Double_t *errors)
{
   if (!errors) return;
   for (Int_t i = 0; i < GetNpar(); i++) fParErrors[i] = errors[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Set limits for parameter ipar.
///
/// The specified limits will be used in a fit operation
/// when the option "B" is specified (Bounds).
/// To fix a parameter, use TF1::FixParameter

void TF1::SetParLimits(Int_t ipar, Double_t parmin, Double_t parmax)
{
   Int_t npar = GetNpar();
   if (ipar < 0 || ipar > npar - 1) return;
   if (int(fParMin.size()) != npar) {
      fParMin.resize(npar);
   }
   if (int(fParMax.size()) != npar) {
      fParMax.resize(npar);
   }
   fParMin[ipar] = parmin;
   fParMax[ipar] = parmax;
}


////////////////////////////////////////////////////////////////////////////////
/// Initialize the upper and lower bounds to draw the function.
///
/// The function range is also used in an histogram fit operation
/// when the option "R" is specified.

void TF1::SetRange(Double_t xmin, Double_t xmax)
{
   fXmin = xmin;
   fXmax = xmax;
   if (fType == EFType::kCompositionFcn && fComposition) {
      fComposition->SetRange(xmin, xmax); // automatically updates sub-functions
   }
   Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Restore value of function saved at point

void TF1::SetSavedPoint(Int_t point, Double_t value)
{
   if (fSave.size() == 0) {
      fSave.resize(fNpx + 3);
   }
   if (point < 0 || point >= int(fSave.size())) return;
   fSave[point] = value;
}


////////////////////////////////////////////////////////////////////////////////
/// Set function title
///  if title has the form "fffffff;xxxx;yyyy", it is assumed that
///  the function title is "fffffff" and "xxxx" and "yyyy" are the
///  titles for the X and Y axis respectively.

void TF1::SetTitle(const char *title)
{
   if (!title) return;
   fTitle = title;
   if (!fHistogram) return;
   fHistogram->SetTitle(title);
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TF1::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      // process new version with new TFormula class which is contained in TF1
      //printf("reading TF1....- version  %d..\n",v);

      if (v > 7) {
         // new classes with new TFormula
         // need to register the objects
         b.ReadClassBuffer(TF1::Class(), this, v, R__s, R__c);
         if (!TestBit(kNotGlobal)) {
            R__LOCKGUARD(gROOTMutex);
            gROOT->GetListOfFunctions()->Add(this);
         }
         if (v >= 10 && v < 11)
            fComposition = std::unique_ptr<TF1AbsComposition>(fComposition_ptr);
         return;
      } else {
         ROOT::v5::TF1Data fold;
         //printf("Reading TF1 as v5::TF1Data- version %d \n",v);
         fold.Streamer(b, v, R__s, R__c, TF1::Class());
         // convert old TF1 to new one
         ((TF1v5Convert *)this)->Convert(fold);
      }
   }

   // Writing
   else {
      Int_t saved = 0;
      // save not-formula functions as array of points
      if (fType > 0 && fSave.empty() && fType != EFType::kCompositionFcn) {
         saved = 1;
         Save(fXmin, fXmax, 0, 0, 0, 0);
      }
      b.WriteClassBuffer(TF1::Class(), this);

      // clear vector contents
      if (saved) {
         fSave.clear();
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Called by functions such as SetRange, SetNpx, SetParameters
/// to force the deletion of the associated histogram or Integral

void TF1::Update()
{
   delete fHistogram;
   fHistogram = 0;
   if (!fIntegral.empty()) {
      fIntegral.clear();
      fAlpha.clear();
      fBeta.clear();
      fGamma.clear();
   }
   if (fNormalized) {
      // need to compute the integral of the not-normalized function
      fNormalized = false;
      fNormIntegral = Integral(fXmin, fXmax, 0.0);
      fNormalized = true;
   } else
      fNormIntegral = 0;

   // std::vector<double>x(fNdim);
   // if ((fType == 1) && !fFunctor->Empty())  (*fFunctor)x.data(), (Double_t*)fParams);
   if (fType == EFType::kCompositionFcn && fComposition) {
      // double-check that the parameters are correct
      fComposition->SetParameters(GetParameters());

      fComposition->Update(); // should not be necessary, but just to be safe
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to set the global flag to reject points
/// the fgRejectPoint global flag is tested by all fit functions
/// if TRUE the point is not included in the fit.
/// This flag can be set by a user in a fitting function.
/// The fgRejectPoint flag is reset by the TH1 and TGraph fitting functions.

void TF1::RejectPoint(Bool_t reject)
{
   fgRejectPoint = reject;
}


////////////////////////////////////////////////////////////////////////////////
/// See TF1::RejectPoint above

Bool_t TF1::RejectedPoint()
{
   return fgRejectPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Return nth moment of function between a and b
///
/// See TF1::Integral() for parameter definitions

Double_t TF1::Moment(Double_t n, Double_t a, Double_t b, const Double_t *params, Double_t epsilon)
{
   // wrapped function in interface for integral calculation
   // using abs value of integral

   TF1_EvalWrapper func(this, params, kTRUE, n);

   ROOT::Math::GaussIntegrator giod;

   giod.SetFunction(func);
   giod.SetRelTolerance(epsilon);

   Double_t norm =  giod.Integral(a, b);
   if (norm == 0) {
      Error("Moment", "Integral zero over range");
      return 0;
   }

   // calculate now integral of x^n f(x)
   // wrapped the member function EvalNum in  interface required by integrator using the functor class
   ROOT::Math::Functor1D xnfunc(&func, &TF1_EvalWrapper::EvalNMom);
   giod.SetFunction(xnfunc);

   Double_t res = giod.Integral(a, b) / norm;

   return res;
}


////////////////////////////////////////////////////////////////////////////////
/// Return nth central moment of function between a and b
/// (i.e the n-th moment around the mean value)
///
/// See TF1::Integral() for parameter definitions
///
/// \author Gene Van Buren <gene@bnl.gov>

Double_t TF1::CentralMoment(Double_t n, Double_t a, Double_t b, const Double_t *params, Double_t epsilon)
{
   TF1_EvalWrapper func(this, params, kTRUE, n);

   ROOT::Math::GaussIntegrator giod;

   giod.SetFunction(func);
   giod.SetRelTolerance(epsilon);

   Double_t norm =  giod.Integral(a, b);
   if (norm == 0) {
      Error("Moment", "Integral zero over range");
      return 0;
   }

   // calculate now integral of xf(x)
   // wrapped the member function EvalFirstMom in  interface required by integrator using the functor class
   ROOT::Math::Functor1D xfunc(&func, &TF1_EvalWrapper::EvalFirstMom);
   giod.SetFunction(xfunc);

   // estimate of mean value
   Double_t xbar = giod.Integral(a, b) / norm;

   // use different mean value in function wrapper
   func.fX0 = xbar;
   ROOT::Math::Functor1D xnfunc(&func, &TF1_EvalWrapper::EvalNMom);
   giod.SetFunction(xnfunc);

   Double_t res = giod.Integral(a, b) / norm;
   return res;
}


//______________________________________________________________________________
// some useful static utility functions to compute sampling points for IntegralFast
////////////////////////////////////////////////////////////////////////////////
/// Type safe interface (static method)
/// The number of sampling points are taken from the TGraph

#ifdef INTHEFUTURE
void TF1::CalcGaussLegendreSamplingPoints(TGraph *g, Double_t eps)
{
   if (!g) return;
   CalcGaussLegendreSamplingPoints(g->GetN(), g->GetX(), g->GetY(), eps);
}


////////////////////////////////////////////////////////////////////////////////
/// Type safe interface (static method)
/// A TGraph is created with new with num points and the pointer to the
/// graph is returned by the function. It is the responsibility of the
/// user to delete the object.
/// if num is invalid (<=0) NULL is returned

TGraph *TF1::CalcGaussLegendreSamplingPoints(Int_t num, Double_t eps)
{
   if (num <= 0)
      return 0;

   TGraph *g = new TGraph(num);
   CalcGaussLegendreSamplingPoints(g->GetN(), g->GetX(), g->GetY(), eps);
   return g;
}
#endif


////////////////////////////////////////////////////////////////////////////////
/// Type: unsafe but fast interface filling the arrays x and w (static method)
///
/// Given the number of sampling points this routine fills the arrays x and w
/// of length num, containing the abscissa and weight of the Gauss-Legendre
/// n-point quadrature formula.
///
/// Gauss-Legendre:
/** \f[
          W(x)=1  -1<x<1 \\
          (j+1)P_{j+1} = (2j+1)xP_j-jP_{j-1}
    \f]
**/
/// num is the number of sampling points (>0)
/// x and w are arrays of size num
/// eps is the relative precision
///
/// If num<=0 or eps<=0 no action is done.
///
/// Reference: Numerical Recipes in C, Second Edition

void TF1::CalcGaussLegendreSamplingPoints(Int_t num, Double_t *x, Double_t *w, Double_t eps)
{
   // This function is just kept like this for backward compatibility!

   ROOT::Math::GaussLegendreIntegrator gli(num, eps);
   gli.GetWeightVectors(x, w);


}


/** \class TF1Parameters
TF1 Parameters class
*/

////////////////////////////////////////////////////////////////////////////////
/// Returns the parameter number given a name
/// not very efficient but list of parameters is typically small
/// could use a map if needed

Int_t TF1Parameters::GetParNumber(const char *name) const
{
   for (unsigned int i = 0; i < fParNames.size(); ++i) {
      if (fParNames[i] == std::string(name)) return i;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter values

void  TF1Parameters::SetParameters(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4,
                                   Double_t p5, Double_t p6, Double_t p7, Double_t p8,
                                   Double_t p9, Double_t p10)
{
   unsigned int npar = fParameters.size();
   if (npar > 0) fParameters[0] = p0;
   if (npar > 1) fParameters[1] = p1;
   if (npar > 2) fParameters[2] = p2;
   if (npar > 3) fParameters[3] = p3;
   if (npar > 4) fParameters[4] = p4;
   if (npar > 5) fParameters[5] = p5;
   if (npar > 6) fParameters[6] = p6;
   if (npar > 7) fParameters[7] = p7;
   if (npar > 8) fParameters[8] = p8;
   if (npar > 9) fParameters[9] = p9;
   if (npar > 10) fParameters[10] = p10;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter names

void TF1Parameters::SetParNames(const char *name0, const char *name1, const char *name2, const char *name3,
                                const char *name4, const char *name5, const char *name6, const char *name7,
                                const char *name8, const char *name9, const char *name10)
{
   unsigned int npar = fParNames.size();
   if (npar > 0) fParNames[0] = name0;
   if (npar > 1) fParNames[1] = name1;
   if (npar > 2) fParNames[2] = name2;
   if (npar > 3) fParNames[3] = name3;
   if (npar > 4) fParNames[4] = name4;
   if (npar > 5) fParNames[5] = name5;
   if (npar > 6) fParNames[6] = name6;
   if (npar > 7) fParNames[7] = name7;
   if (npar > 8) fParNames[8] = name8;
   if (npar > 9) fParNames[9] = name9;
   if (npar > 10) fParNames[10] = name10;
}
