// @(#)root/hist:$Id$
// Author: Maciej Zimnoch 30/09/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#if __cplusplus >= 201103L
#define ROOT_CPLUSPLUS11 1
#endif

#include "TROOT.h"
#include "TClass.h"
#include "TMethod.h"
#include "TMath.h"
#include "TF1.h"
#include "TMethodCall.h"
#include <TBenchmark.h>
#include "TError.h"
#include "TInterpreter.h"
#include "TFormula.h"
#include "TMap.h"
#include "TRegexp.h"
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <functional>

using namespace std;

// #define __STDC_LIMIT_MACROS
// #define __STDC_CONSTANT_MACROS

// #include  "cling/Interpreter/Interpreter.h"
// #include  "cling/Interpreter/Value.h"
// #include  "cling/Interpreter/StoredValueRef.h"


#ifdef WIN32
#pragma optimize("",off)
#endif
#include "v5/TFormula.h"

ClassImp(TFormula);

/** \class TFormula  TFormula.h "inc/TFormula.h"
    \ingroup Hist
    The Formula class

    This is a new version of the TFormula class based on Cling.
    This class is not 100% backward compatible with the old TFormula class, which is still available in ROOT as
    `ROOT::v5::TFormula`. Some of the TFormula member funtions available in version 5, such as
    `Analyze` and `AnalyzeFunction` are not available in the new TFormula.
    On the other hand formula expressions which were valid in version 5 are still valid in TFormula version 6

    This class has been implemented during Google Summer of Code 2013 by Maciej Zimnoch.

    ### Example of valid expressions:

    - `sin(x)/x`
    - `[0]*sin(x) + [1]*exp(-[2]*x)`
    - `x + y**2`
    - `x^2 + y^2`
    - `[0]*pow([1],4)`
    - `2*pi*sqrt(x/y)`
    - `gaus(0)*expo(3)  + ypol3(5)*x`
    - `gausn(0)*expo(3) + ypol3(5)*x`

    In the last example above:

    - `gaus(0)` is a substitute for `[0]*exp(-0.5*((x-[1])/[2])**2)`
    and (0) means start numbering parameters at 0
    - `gausn(0)` is a substitute for `[0]*exp(-0.5*((x-[1])/[2])**2)/(sqrt(2*pi)*[2]))`
    and (0) means start numbering parameters at 0
    - `expo(3)` is a substitute for `exp([3]+[4]*x`
    - `pol3(5)` is a substitute for `par[5]+par[6]*x+par[7]*x**2+par[8]*x**3`
    (`PolN` stands for Polynomial of degree N)

    `TMath` functions can be part of the expression, eg:

    - `TMath::Landau(x)*sin(x)`
    - `TMath::Erf(x)`

    Formula may contain constans, eg:

    - `sqrt2`
    - `e`
    - `pi`
    - `ln10`
    - `infinity`

    and more.

    Comparisons operators are also supported `(&amp;&amp;, ||, ==, &lt;=, &gt;=, !)`

    Examples:

    `sin(x*(x&lt;0.5 || x&gt;1))`

    If the result of a comparison is TRUE, the result is 1, otherwise 0.

    Already predefined names can be given. For example, if the formula

    `TFormula old("old",sin(x*(x&lt;0.5 || x&gt;1)))`

    one can assign a name to the formula. By default the name of the object = title = formula itself.

    `TFormula new("new","x*old")`

    is equivalent to:

    `TFormula new("new","x*sin(x*(x&lt;0.5 || x&gt;1))")`

    The class supports unlimited numer of variables and parameters.
    By default the names which can be used for the variables are `x,y,z,t` or
    `x[0],x[1],x[2],x[3],....x[N]` for N-dimensionals formula.

    This class is not anymore the base class for the function classes `TF1`, but it has now
    adata member of TF1 which can be access via `TF1::GetFormula`.

    \class TFormulaFunction
    Helper class for TFormula

    \class TFormulaVariable
    Another helper class for TFormula

    \class TFormulaParamOrder
    Functor defining the parameter order
*/

// prefix used for function name passed to Cling
static const TString gNamePrefix = "TFormula__";

// static map of function pointers and expressions
//static std::unordered_map<std::string,  TInterpreter::CallFuncIFacePtr_t::Generic_t> gClingFunctions = std::unordered_map<TString,  TInterpreter::CallFuncIFacePtr_t::Generic_t>();
static std::unordered_map<std::string,  void *> gClingFunctions = std::unordered_map<std::string,  void * >();

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsOperator(const char c)
{
   // operator ":" must be handled separately
   char ops[] = { '+','^', '-','/','*','<','>','|','&','!','=','?'};
   Int_t opsLen = sizeof(ops)/sizeof(char);
   for(Int_t i = 0; i < opsLen; ++i)
      if(ops[i] == c)
         return true;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsBracket(const char c)
{
   // Note that square brackets do not count as brackets here!!!
   char brackets[] = { ')','(','{','}'};
   Int_t bracketsLen = sizeof(brackets)/sizeof(char);
   for(Int_t i = 0; i < bracketsLen; ++i)
      if(brackets[i] == c)
         return true;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsFunctionNameChar(const char c)
{
   return !IsBracket(c) && !IsOperator(c) && c != ',' && c != ' ';
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsDefaultVariableName(const TString &name)
{
   return name == "x" || name == "z" || name == "y" || name == "t";
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsScientificNotation(const TString & formula, int i)
{
   // check if the character at position i  is part of a scientific notation
   if ( (formula[i] == 'e' || formula[i] == 'E')  &&  (i > 0 && i <  formula.Length()-1) )  {
      // handle cases:  2e+3 2e-3 2e3 and 2.e+3
      if ( (isdigit(formula[i-1]) || formula[i-1] == '.') && ( isdigit(formula[i+1]) || formula[i+1] == '+' || formula[i+1] == '-' ) )
         return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::IsHexadecimal(const TString & formula, int i)
{
   // check if the character at position i  is part of a scientific notation
   if ( (formula[i] == 'x' || formula[i] == 'X')  &&  (i > 0 && i <  formula.Length()-1) && formula[i-1] == '0')  {
      if (isdigit(formula[i+1]) )
         return true;
      static char hex_values[12] = { 'a','A', 'b','B','c','C','d','D','e','E','f','F'};
      for (int jjj = 0; jjj < 12; ++jjj) {
         if (formula[i+1] == hex_values[jjj])
            return true;
      }
   }
   // else
   //    return false;
   //    // handle cases:  2e+3 2e-3 2e3 and 2.e+3
   //    if ( (isdigit(formula[i-1]) || formula[i-1] == '.') && ( isdigit(formula[i+1]) || formula[i+1] == '+' || formula[i+1] == '-' ) )
   //       return true;
   // }
   return false;
}
////////////////////////////////////////////////////////////////////////////
// check is given position is in a parameter name i.e. within "[ ]"
////
Bool_t TFormula::IsAParameterName(const TString & formula, int pos) {

   Bool_t foundOpenParenthesis = false;
   if (pos == 0 || pos == formula.Length()-1) return false; 
   for (int i = pos-1; i >=0; i--) {
      if (formula[i] == ']' ) return false;
      if (formula[i] == '[' ) {
         foundOpenParenthesis = true;
         break;
      }
   }
   if (!foundOpenParenthesis ) return false; 

   // search after the position
   for (int i = pos+1; i < formula.Length(); i++) {
      if (formula[i] == ']' ) return true;
   }
   return false; 
}


////////////////////////////////////////////////////////////////////////////////
bool TFormulaParamOrder::operator() (const TString& a, const TString& b) const {
   // implement comparison used to set parameter orders in TFormula
   // want p2 to be before p10

   // Returns true if (a < b), meaning a comes before b, and false if (a >= b)

   TRegexp numericPattern("p?[0-9]+");
   Ssiz_t *len = new Ssiz_t(); // buffer to store length of regex match

   int patternStart = numericPattern.Index(a, len);
   bool aNumeric = (patternStart == 0 && *len == a.Length());

   patternStart = numericPattern.Index(b, len);
   bool bNumeric = (patternStart == 0 && *len == b.Length());

   if (aNumeric && !bNumeric)
      return true; // assume a (numeric) is always before b (not numeric)
   else if (!aNumeric && bNumeric)
      return false; // b comes before a
   else if (!aNumeric && !bNumeric)
      return a < b;
   else {
      int aInt = (a[0] == 'p') ? TString(a(1, a.Length())).Atoi() : a.Atoi();
      int bInt = (b[0] == 'p') ? TString(b(1, b.Length())).Atoi() : b.Atoi();
      return aInt < bInt;
   }

}

////////////////////////////////////////////////////////////////////////////////
void TFormula::ReplaceAllNames(TString &formula, map<TString, TString> &substitutions)
{
   /// Apply the name substitutions to the formula, doing all replacements in one pass

   for (int i = 0; i < formula.Length(); i++) {
      // start of name
      // (a little subtle, since we want to match names like "{V0}" and "[0]")
      if (isalpha(formula[i]) || formula[i] == '{' || formula[i] == '[') {
         int j; // index to end of name
         for (j = i + 1;
              j < formula.Length() && (IsFunctionNameChar(formula[j]) // square brackets are function name chars
                                       || (formula[i] == '{' && formula[j] == '}'));
              j++)
            ;
         TString name = (TString)formula(i, j - i);

         // std::cout << "Looking for name: " << name << std::endl;

         // if we find the name, do the substitution
         if (substitutions.find(name) != substitutions.end()) {
            formula.Replace(i, name.Length(), "(" + substitutions[name] + ")");
            i += substitutions[name].Length() + 2 - 1; // +2 for parentheses
            // std::cout << "made substitution: " << name << " to " << substitutions[name] << std::endl;
         } else if (isalpha(formula[i])) {
            // if formula[i] is alpha, can skip to end of candidate name, otherwise, we'll just
            // move one character ahead and try again
            i += name.Length() - 1;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
TFormula::TFormula()
{
   fName = "";
   fTitle = "";
   fClingInput = "";
   fReadyToExecute = false;
   fClingInitialized = false;
   fAllParametersSetted = false;
   fMethod = 0;
   fNdim = 0;
   fNpar = 0;
   fNumber = 0;
   fClingName = "";
   fFormula = "";
   fLambdaPtr = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
static bool IsReservedName(const char* name){
   if (strlen(name)!=1) return false;
   for (auto const & specialName : {"x","y","z","t"}){
      if (strcmp(name,specialName)==0) return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
TFormula::~TFormula()
{

   // N.B. a memory leak may happen if user set bit after constructing the object,
   // Setting of bit should be done only internally
   if (!TestBit(TFormula::kNotGlobal) && gROOT ) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfFunctions()->Remove(this);
   }

   if (fMethod) {
      fMethod->Delete();
   }
   int nLinParts = fLinearParts.size();
   if (nLinParts > 0) {
      for (int i = 0; i < nLinParts; ++i) delete fLinearParts[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
TFormula::TFormula(const char *name, const char *formula, bool addToGlobList)   :
   TNamed(name,formula),
   fClingInput(formula),fFormula(formula)
{
   fReadyToExecute = false;
   fClingInitialized = false;
   fMethod = 0;
   fNdim = 0;
   fNpar = 0;
   fNumber = 0;
   fMethod = 0;
   fLambdaPtr = nullptr;

   FillDefaults();


   if (addToGlobList && gROOT) {
      TFormula *old = 0;
      R__LOCKGUARD(gROOTMutex);
      old = dynamic_cast<TFormula*> ( gROOT->GetListOfFunctions()->FindObject(name) );
      if (old)
         gROOT->GetListOfFunctions()->Remove(old);
      if (IsReservedName(name))
         Error("TFormula","The name %s is reserved as a TFormula variable name.\n",name);
      else
         gROOT->GetListOfFunctions()->Add(this);
   }
   SetBit(kNotGlobal,!addToGlobList);

   //fName = gNamePrefix + name;  // is this needed

   // do not process null formulas.
   if (!fFormula.IsNull() ) {
      PreProcessFormula(fFormula);

      PrepareFormula(fFormula);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a full compileable C++ expression

TFormula::TFormula(const char *name, const char *formula, int ndim, int npar, bool addToGlobList)   :
   TNamed(name,formula),
   fClingInput(formula),fFormula(formula)
{
   fReadyToExecute = false;
   fClingInitialized = false;
   fNpar = 0;
   fMethod = 0;
   fNumber = 0;
   fLambdaPtr = nullptr;
   fFuncPtr = nullptr;


   fNdim = ndim;
   for (int i = 0; i < npar; ++i) {
      DoAddParameter(TString::Format("p%d",i), 0, false);
   }
   fAllParametersSetted = true;
   assert (fNpar == npar);

   bool ret = InitLambdaExpression(formula);

   if (ret)  {

      SetBit(TFormula::kLambda);

      fReadyToExecute = true;

      if (addToGlobList && gROOT) {
         TFormula *old = 0;
         R__LOCKGUARD(gROOTMutex);
         old = dynamic_cast<TFormula*> ( gROOT->GetListOfFunctions()->FindObject(name) );
         if (old)
            gROOT->GetListOfFunctions()->Remove(old);
         if (IsReservedName(name))
            Error("TFormula","The name %s is reserved as a TFormula variable name.\n",name);
         else
            gROOT->GetListOfFunctions()->Add(this);
      }
      SetBit(kNotGlobal,!addToGlobList);
   }
   else
      Error("TFormula","Syntax error in building the lambda expression %s", formula );
}

////////////////////////////////////////////////////////////////////////////////
TFormula::TFormula(const TFormula &formula) :
   TNamed(formula.GetName(),formula.GetTitle())
{
   fReadyToExecute = false;
   fClingInitialized = false;
   fMethod = 0;
   fNdim = formula.GetNdim();
   fNpar = formula.GetNpar();
   fNumber = formula.GetNumber();
   fFormula = formula.GetExpFormula();   // returns fFormula in case of Lambda's
   fLambdaPtr = nullptr;
   fFuncPtr = nullptr;

   // case of function based on a C++  expression (lambda's) which is ready to be compiled
   if (formula.fLambdaPtr && formula.TestBit(TFormula::kLambda)) {

      fClingInput = fFormula;
      fParams = formula.fParams;
      fClingParameters = formula.fClingParameters;
      fAllParametersSetted = formula.fAllParametersSetted;

      bool ret = InitLambdaExpression(fFormula);

      if (ret)  {
         SetBit(TFormula::kLambda);
         fReadyToExecute = true;
      }
      else
         Error("TFormula","Syntax error in building the lambda expression %s", fFormula.Data() );

   }
   else {

      FillDefaults();

      PreProcessFormula(fFormula);
      PrepareFormula(fFormula);
   }


   if (!TestBit(TFormula::kNotGlobal) && gROOT ) {
      R__LOCKGUARD(gROOTMutex);
      TFormula *old = (TFormula*)gROOT->GetListOfFunctions()->FindObject(formula.GetName());
      if (old)
         gROOT->GetListOfFunctions()->Remove(old);

      if (IsReservedName(formula.GetName())) {
         Error("TFormula","The name %s is reserved as a TFormula variable name.\n",formula.GetName());
      } else
         gROOT->GetListOfFunctions()->Add(this);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// = operator.

TFormula& TFormula::operator=(const TFormula &rhs)
{

   if (this != &rhs) {
      rhs.Copy(*this);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
Bool_t TFormula::InitLambdaExpression(const char * formula) {

   std::string lambdaExpression = formula;

   // check if formula exist already in the map
   {
      R__LOCKGUARD(gROOTMutex);

      auto funcit = gClingFunctions.find(lambdaExpression);
      if (funcit != gClingFunctions.end() ) {
         fLambdaPtr = funcit->second;
         fClingInitialized = true;
         return true;
      }
   }

   // set the cling name using hash of the static formulae map
   auto hasher = gClingFunctions.hash_function();
   TString lambdaName = TString::Format("lambda__id%zu", hasher(lambdaExpression) );

   //lambdaExpression = TString::Format("[&](double * x, double *){ return %s ;}",formula);
   //TString lambdaName = TString::Format("mylambda_%s",GetName() );
   TString lineExpr = TString::Format("std::function<double(double*,double*)> %s = %s ;",lambdaName.Data(), lambdaExpression.c_str() );
   gInterpreter->ProcessLine(lineExpr);
   fLambdaPtr = (void*) gInterpreter->ProcessLine(TString(lambdaName)+TString(";"));  // add ; to avoid printing
   if (fLambdaPtr != nullptr) {
      R__LOCKGUARD(gROOTMutex);
      gClingFunctions.insert ( std::make_pair ( lambdaExpression, fLambdaPtr) );
      fClingInitialized = true;
      return true;
   }
   fClingInitialized = false;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Compile the given expression with Cling
/// backward compatibility method to be used in combination with the empty constructor
/// if no expression is given , the current stored formula (retrieved with GetExpFormula()) or the title  is used.
/// return 0 if the formula compilation is successfull

Int_t TFormula::Compile(const char *expression)
{
   TString formula = expression;
   if (formula.IsNull() ) {
      formula = fFormula;
      if (formula.IsNull() ) formula = GetTitle();
   }

   if (formula.IsNull() ) return -1;

   // do not re-process if it was done before
   if (IsValid() && formula == fFormula ) return 0;

   // clear if a formula was already existing
   if (!fFormula.IsNull() ) Clear();

   fFormula = formula;

   if (TestBit(TFormula::kLambda) ) {
      bool ret = InitLambdaExpression(fFormula);
      return (ret) ? 0 : 1;
   }

   if (fVars.empty() ) FillDefaults();
   // prepare the formula for Cling
   //printf("compile: processing formula %s\n",fFormula.Data() );
   PreProcessFormula(fFormula);
   // pass formula in CLing
   bool ret = PrepareFormula(fFormula);

   return (ret) ? 0 : 1;
}

////////////////////////////////////////////////////////////////////////////////
void TFormula::Copy(TObject &obj) const
{
   TNamed::Copy(obj);
   // need to copy also cling parameters
   TFormula & fnew = dynamic_cast<TFormula&>(obj);

   fnew.fClingParameters = fClingParameters;
   fnew.fClingVariables = fClingVariables;

   fnew.fFuncs = fFuncs;
   fnew.fVars = fVars;
   fnew.fParams = fParams;
   fnew.fConsts = fConsts;
   fnew.fFunctionsShortcuts = fFunctionsShortcuts;
   fnew.fFormula  = fFormula;
   fnew.fNdim = fNdim;
   fnew.fNpar = fNpar;
   fnew.fNumber = fNumber;
   fnew.SetParameters(GetParameters());
   // copy Linear parts (it is a vector of TFormula pointers) needs to be copied one by one
   // looping at all the elements
   // delete first previous elements
   int nLinParts = fnew.fLinearParts.size();
   if (nLinParts > 0) {
      for (int i = 0; i < nLinParts; ++i) delete fnew.fLinearParts[i];
      fnew.fLinearParts.clear();
   }
   // old size that needs to be copied
   nLinParts = fLinearParts.size();
   if (nLinParts > 0) {
      fnew.fLinearParts.reserve(nLinParts);
      for (int i = 0; i < nLinParts; ++i) {
         TFormula * linearNew = new TFormula();
         TFormula * linearOld = (TFormula*) fLinearParts[i];
         if (linearOld) {
            linearOld->Copy(*linearNew);
            fnew.fLinearParts.push_back(linearNew);
         }
         else
            Warning("Copy","Function %s - expr %s has a dummy linear part %d",GetName(),GetExpFormula().Data(),i);
      }
   }

   fnew.fClingInput = fClingInput;
   fnew.fReadyToExecute = fReadyToExecute;
   fnew.fClingInitialized = fClingInitialized;
   fnew.fAllParametersSetted = fAllParametersSetted;
   fnew.fClingName = fClingName;

   // case of function based on a C++  expression (lambda's) which is ready to be compiled
   if (fLambdaPtr && TestBit(TFormula::kLambda)) {

      bool ret = fnew.InitLambdaExpression(fnew.fFormula);
      if (ret)  {
         fnew.SetBit(TFormula::kLambda);
         fnew.fReadyToExecute = true;
      }
      else {
         Error("TFormula","Syntax error in building the lambda expression %s", fFormula.Data() );
         fnew.fReadyToExecute = false;
      }
   }
   else if (fMethod) {
      if (fnew.fMethod) delete fnew.fMethod;
      // use copy-constructor of TMethodCall
      TMethodCall *m = new TMethodCall(*fMethod);
      fnew.fMethod  = m;
   }

   fnew.fFuncPtr = fFuncPtr;

}

////////////////////////////////////////////////////////////////////////////////
/// Clear the formula setting expression to empty and reset the variables and
/// parameters containers.

void TFormula::Clear(Option_t * )
{
   fNdim = 0;
   fNpar = 0;
   fNumber = 0;
   fFormula = "";
   fClingName = "";


   if(fMethod) fMethod->Delete();
   fMethod = nullptr;

   fClingVariables.clear();
   fClingParameters.clear();
   fReadyToExecute = false;
   fClingInitialized = false;
   fAllParametersSetted = false;
   fFuncs.clear();
   fVars.clear();
   fParams.clear();
   fConsts.clear();
   fFunctionsShortcuts.clear();

   // delete linear parts
   int nLinParts = fLinearParts.size();
   if (nLinParts > 0) {
      for (int i = 0; i < nLinParts; ++i) delete fLinearParts[i];
   }
   fLinearParts.clear();

}

////////////////////////////////////////////////////////////////////////////////
/// Sets TMethodCall to function inside Cling environment.
/// TFormula uses it to execute function.
/// After call, TFormula should be ready to evaluate formula.

bool TFormula::PrepareEvalMethod()
{

   if(!fMethod)
      {
      fMethod = new TMethodCall();

      Bool_t hasParameters = (fNpar > 0);
      Bool_t hasVariables = (fNdim > 0);
      TString prototypeArguments = "";
      if (hasVariables) {
         prototypeArguments.Append("Double_t*");
      }
      if (hasVariables && hasParameters) {
         prototypeArguments.Append(",");
      }
      if (hasParameters) {
         prototypeArguments.Append("Double_t*");
      }
      // init method call using real function name (cling name) which is defined in ProcessFormula
      fMethod->InitWithPrototype(fClingName, prototypeArguments);
      if (!fMethod->IsValid()) {
         Error("Eval", "Can't find %s function prototype with arguments %s", fClingName.Data(),
               prototypeArguments.Data());
         return false;
      }

      // not needed anymore since we use the function pointer
      // if(hasParameters)
      // {
      //    Long_t args[2];
      //    args[0] = (Long_t)fClingVariables.data();
      //    args[1] = (Long_t)fClingParameters.data();
      //    fMethod->SetParamPtrs(args,2);
      // }
      // else
      // {
      //    Long_t args[1];
      //    args[0] = (Long_t)fClingVariables.data();
      //    fMethod->SetParamPtrs(args,1);
      // }

      CallFunc_t *callfunc = fMethod->GetCallFunc();
      TInterpreter::CallFuncIFacePtr_t faceptr = gCling->CallFunc_IFacePtr(callfunc);
      fFuncPtr = faceptr.fGeneric;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
///    Inputs formula, transfered to C++ code into Cling

void TFormula::InputFormulaIntoCling()
{

   if (!fClingInitialized && fReadyToExecute && fClingInput.Length() > 0) {
      gCling->Declare(fClingInput);
      fClingInitialized = PrepareEvalMethod();
   }
}

////////////////////////////////////////////////////////////////////////////////
///    Fill structures with default variables, constants and function shortcuts

void TFormula::FillDefaults()
{
   //#ifdef ROOT_CPLUSPLUS11

   const TString defvars[] = { "x","y","z","t"};
   const pair<TString, Double_t> defconsts[] = {{"pi", TMath::Pi()},
                                                {"sqrt2", TMath::Sqrt2()},
                                                {"infinity", TMath::Infinity()},
                                                {"e", TMath::E()},
                                                {"ln10", TMath::Ln10()},
                                                {"loge", TMath::LogE()},
                                                {"c", TMath::C()},
                                                {"g", TMath::G()},
                                                {"h", TMath::H()},
                                                {"k", TMath::K()},
                                                {"sigma", TMath::Sigma()},
                                                {"r", TMath::R()},
                                                {"eg", TMath::EulerGamma()},
                                                {"true", 1},
                                                {"false", 0}};
   // const pair<TString,Double_t> defconsts[] = { {"pi",TMath::Pi()}, {"sqrt2",TMath::Sqrt2()},
   //       {"infinity",TMath::Infinity()}, {"ln10",TMath::Ln10()},
   //       {"loge",TMath::LogE()}, {"true",1},{"false",0} };
   const pair<TString,TString> funShortcuts[] =
      { {"sin","TMath::Sin" },
        {"cos","TMath::Cos" }, {"exp","TMath::Exp"}, {"log","TMath::Log"}, {"log10","TMath::Log10"},
        {"tan","TMath::Tan"}, {"sinh","TMath::SinH"}, {"cosh","TMath::CosH"},
        {"tanh","TMath::TanH"}, {"asin","TMath::ASin"}, {"acos","TMath::ACos"},
        {"atan","TMath::ATan"}, {"atan2","TMath::ATan2"}, {"sqrt","TMath::Sqrt"},
        {"ceil","TMath::Ceil"}, {"floor","TMath::Floor"}, {"pow","TMath::Power"},
        {"binomial","TMath::Binomial"},{"abs","TMath::Abs"},
        {"min","TMath::Min"},{"max","TMath::Max"},{"sign","TMath::Sign" },
        {"sq","TMath::Sq"}
      };

   std::vector<TString> defvars2(10);
   for (int i = 0; i < 9; ++i)
      defvars2[i] = TString::Format("x[%d]",i);

   for (auto var : defvars) {
      int pos = fVars.size();
      fVars[var] = TFormulaVariable(var, 0, pos);
      fClingVariables.push_back(0);
   }
   // add also the variables definesd like x[0],x[1],x[2],...
   // support up to x[9] - if needed extend that to higher value
   // const int maxdim = 10;
   // for (int i = 0; i < maxdim;  ++i) {
   //    TString xvar = TString::Format("x[%d]",i);
   //    fVars[xvar] =  TFormulaVariable(xvar,0,i);
   //    fClingVariables.push_back(0);
   // }

   for (auto con : defconsts) {
      fConsts[con.first] = con.second;
   }
   for (auto fun : funShortcuts) {
      fFunctionsShortcuts[fun.first] = fun.second;
   }

   /*** - old code to support C++03
        #else

        TString  defvarsNames[] = {"x","y","z","t"};
        Int_t    defvarsLength = sizeof(defvarsNames)/sizeof(TString);

        TString  defconstsNames[] =
   {"pi","sqrt2","infinity","e","ln10","loge","c","g","h","k","sigma","r","eg","true","false"};
        Double_t defconstsValues[] =
   {TMath::Pi(),TMath::Sqrt2(),TMath::Infinity(),TMath::E(),TMath::Ln10(),TMath::LogE(),
        TMath::C(),TMath::G(),TMath::H(),TMath::K(),TMath::Sigma(),TMath::R(),TMath::EulerGamma(), 1, 0};
        Int_t    defconstsLength = sizeof(defconstsNames)/sizeof(TString);

        TString  funShortcutsNames[] =
   {"sin","cos","exp","log","tan","sinh","cosh","tanh","asin","acos","atan","atan2","sqrt",
        "ceil","floor","pow","binomial","abs"};
        TString  funShortcutsExtendedNames[] =
   {"TMath::Sin","TMath::Cos","TMath::Exp","TMath::Log","TMath::Tan","TMath::SinH",
        "TMath::CosH","TMath::TanH","TMath::ASin","TMath::ACos","TMath::ATan","TMath::ATan2",
        "TMath::Sqrt","TMath::Ceil","TMath::Floor","TMath::Power","TMath::Binomial","TMath::Abs"};
        Int_t    funShortcutsLength = sizeof(funShortcutsNames)/sizeof(TString);

        for(Int_t i = 0; i < defvarsLength; ++i)
        {
        TString var = defvarsNames[i];
        Double_t value = 0;
        unsigned int pos = fVars.size();
        fVars[var] = TFormulaVariable(var,value,pos);
        fClingVariables.push_back(value);
        }

        for(Int_t i = 0; i < defconstsLength; ++i)
        {
        fConsts[defconstsNames[i]] = defconstsValues[i];
        }
        for(Int_t i = 0; i < funShortcutsLength; ++i)
        {
        pair<TString,TString> fun(funShortcutsNames[i],funShortcutsExtendedNames[i]);
        fFunctionsShortcuts[fun.first] = fun.second;
        }

        #endif
   ***/
}

////////////////////////////////////////////////////////////////////////////////
///    Handling polN
///    If before 'pol' exist any name, this name will be treated as variable used in polynomial
///    eg.
///    varpol2(5) will be replaced with: [5] + [6]*var + [7]*var^2
///    Empty name is treated like variable x.

void TFormula::HandlePolN(TString &formula)
{
   Int_t polPos = formula.Index("pol");
   while(polPos != kNPOS  && !IsAParameterName(formula,polPos) )
      {

      Bool_t defaultVariable = false;
      TString variable;
      Int_t openingBracketPos = formula.Index('(', polPos);
      Bool_t defaultCounter = openingBracketPos == kNPOS;
      Bool_t defaultDegree = true;
      Int_t degree, counter;
      TString sdegree;
      if (!defaultCounter) {
         // veryfy first of opening parenthesis belongs to pol expression
         // character between 'pol' and '(' must all be digits
         sdegree = formula(polPos + 3, openingBracketPos - polPos - 3);
         if (!sdegree.IsDigit())
            defaultCounter = true;
      }
      if (!defaultCounter) {
         degree = sdegree.Atoi();
         counter = TString(formula(openingBracketPos + 1, formula.Index(')', polPos) - openingBracketPos)).Atoi();
      } else {
         Int_t temp = polPos + 3;
         while (temp < formula.Length() && isdigit(formula[temp])) {
            defaultDegree = false;
            temp++;
         }
         degree = TString(formula(polPos + 3, temp - polPos - 3)).Atoi();
         counter = 0;
      }

      TString replacement = TString::Format("[%d]", counter);
      if (polPos - 1 < 0 || !IsFunctionNameChar(formula[polPos - 1]) || formula[polPos - 1] == ':') {
         variable = "x";
         defaultVariable = true;
      } else {
         Int_t tmp = polPos - 1;
         while (tmp >= 0 && IsFunctionNameChar(formula[tmp]) && formula[tmp] != ':') {
            tmp--;
         }
         variable = formula(tmp + 1, polPos - (tmp + 1));
      }
      Int_t param = counter + 1;
      Int_t tmp = 1;
      while (tmp <= degree) {
         if (tmp > 1)
            replacement.Append(TString::Format("+[%d]*%s^%d", param, variable.Data(), tmp));
         else
            replacement.Append(TString::Format("+[%d]*%s", param, variable.Data()));
         param++;
         tmp++;
      }
      // add paranthesis before and after
      if (degree > 0) {
         replacement.Insert(0, '(');
         replacement.Append(')');
         }
         TString pattern;
         if (defaultCounter && !defaultDegree) {
            pattern = TString::Format("%spol%d", (defaultVariable ? "" : variable.Data()), degree);
         } else if (defaultCounter && defaultDegree) {
            pattern = TString::Format("%spol", (defaultVariable ? "" : variable.Data()));
         } else {
            pattern = TString::Format("%spol%d(%d)", (defaultVariable ? "" : variable.Data()), degree, counter);
         }

         if (!formula.Contains(pattern)) {
            Error("HandlePolN", "Error handling polynomial function - expression is %s - trying to replace %s with %s ",
                  formula.Data(), pattern.Data(), replacement.Data());
            break;
         }
         if (formula == pattern) {
            // case of single polynomial
            SetBit(kLinear, 1);
            fNumber = 300 + degree;
         }
         formula.ReplaceAll(pattern, replacement);
         polPos = formula.Index("pol");
      }
}

////////////////////////////////////////////////////////////////////////////////
///    Handling parametrized functions
///    Function can be normalized, and have different variable then x.
///    Variables should be placed in brackets after function name.
///    No brackets are treated like [x].
///    Normalized function has char 'n' after name, eg.
///    gausn[var](0) will be replaced with [0]*exp(-0.5*((var-[1])/[2])^2)/(sqrt(2*pi)*[2])
///
///    Adding function is easy, just follow these rules, and add to
///    `TFormula::FillParametrizedFunctions` defined further below:
///
///    - Key for function map is pair of name and dimension of function
///    - value of key is a pair function body and normalized function body
///    - {Vn} is a place where variable appear, n represents n-th variable from variable list.
///      Count starts from 0.
///    - [num] stands for parameter number.
///      If user pass to function argument 5, num will stand for (5 + num) parameter.
///

void TFormula::HandleParametrizedFunctions(TString &formula)
{
   // define all parametrized functions
   map< pair<TString,Int_t> ,pair<TString,TString> > functions;
   FillParametrizedFunctions(functions);

   map<TString,Int_t> functionsNumbers;
   functionsNumbers["gaus"] = 100;
   functionsNumbers["bigaus"] = 102;
   functionsNumbers["landau"] = 400;
   functionsNumbers["expo"] = 200;
   functionsNumbers["crystalball"] = 500;

   // replace old names xygaus -> gaus[x,y]
   formula.ReplaceAll("xygaus","gaus[x,y]");
   formula.ReplaceAll("xgaus","gaus[x]");
   formula.ReplaceAll("ygaus","gaus[y]");
   formula.ReplaceAll("zgaus","gaus[z]");
   formula.ReplaceAll("xexpo","expo[x]");
   formula.ReplaceAll("yexpo","expo[y]");
   formula.ReplaceAll("zexpo","expo[z]");
   formula.ReplaceAll("xylandau","landau[x,y]");
   formula.ReplaceAll("xyexpo","expo[x,y]");
   // at the moment pre-defined functions have no more than 3 dimensions
   const char * defaultVariableNames[] = { "x","y","z"};

   for (map<pair<TString, Int_t>, pair<TString, TString>>::iterator it = functions.begin(); it != functions.end();
        ++it) {

      TString funName = it->first.first;
      Int_t funDim = it->first.second;
      Int_t funPos = formula.Index(funName);

      // std::cout << formula << " ---- " << funName << "  " << funPos << std::endl;
      while (funPos != kNPOS && !IsAParameterName(formula, funPos)) {

         // should also check that function is not something else (e.g. exponential - parse the expo)
         Int_t lastFunPos = funPos + funName.Length();

         // check that first and last character is not alphanumeric
         Int_t iposBefore = funPos - 1;
         // std::cout << "looping on  funpos is " << funPos << " formula is " << formula << " function " << funName <<
         // std::endl;
         if (iposBefore >= 0) {
            assert(iposBefore < formula.Length());
            if (isalpha(formula[iposBefore])) {
               // std::cout << "previous character for function " << funName << " is " << formula[iposBefore] << "- skip
               // " << std::endl;
               funPos = formula.Index(funName, lastFunPos);
               continue;
            }
         }

         Bool_t isNormalized = false;
         if (lastFunPos < formula.Length()) {
            // check if function is normalized by looking at "n" character after function name (e.g. gausn)
            isNormalized = (formula[lastFunPos] == 'n');
            if (isNormalized)
               lastFunPos += 1;
            if (lastFunPos < formula.Length()) {
               char c = formula[lastFunPos];
               // check if also last character is not alphanumeric or is not an operator and not a parenthesis ( or [.
               // Parenthesis [] are used to express the variables
               if (isalnum(c) || (!IsOperator(c) && c != '(' && c != ')' && c != '[' && c != ']')) {
                  // std::cout << "last character for function " << funName << " is " << c << " skip .." <<  std::endl;
                  funPos = formula.Index(funName, lastFunPos);
                  continue;
               }
            }
         }

         if (isNormalized) {
            SetBit(kNormalized, 1);
         }
         std::vector<TString> variables;
         Int_t dim = 0;
         TString varList = "";
         Bool_t defaultVariables = false;

         // check if function has specified the [...] e.g. gaus[x,y]
         Int_t openingBracketPos = funPos + funName.Length() + (isNormalized ? 1 : 0);
         Int_t closingBracketPos = kNPOS;
         if (openingBracketPos > formula.Length() || formula[openingBracketPos] != '[') {
            dim = funDim;
            variables.resize(dim);
            for (Int_t idim = 0; idim < dim; ++idim)
               variables[idim] = defaultVariableNames[idim];
            defaultVariables = true;
         } else {
            // in case of [..] found, assume they specify all the variables. Use it to get function dimension
            closingBracketPos = formula.Index(']', openingBracketPos);
            varList = formula(openingBracketPos + 1, closingBracketPos - openingBracketPos - 1);
            dim = varList.CountChar(',') + 1;
            variables.resize(dim);
            Int_t Nvar = 0;
            TString varName = "";
            for (Int_t i = 0; i < varList.Length(); ++i) {
               if (IsFunctionNameChar(varList[i])) {
                  varName.Append(varList[i]);
               }
               if (varList[i] == ',') {
                  variables[Nvar] = varName;
                  varName = "";
                  Nvar++;
               }
            }
            if (varName != "") // we will miss last variable
            {
               variables[Nvar] = varName;
            }
         }
         // chech if dimension obtained from [...] is compatible with what is defined in existing pre-defined functions
         // std::cout << " Found dim = " << dim  << " and function dimension is " << funDim << std::endl;
         if (dim != funDim) {
            pair<TString, Int_t> key = make_pair(funName, dim);
            if (functions.find(key) == functions.end()) {
               Error("PreProcessFormula", "Dimension of function %s is detected to be of dimension %d and is not "
                                          "compatible with existing pre-defined function which has dim %d",
                     funName.Data(), dim, funDim);
               return;
            }
            // skip the particular function found - we might find later on the corresponding pre-defined function
            funPos = formula.Index(funName, lastFunPos);
            continue;
         }
         // look now for the (..) brackets to get the parameter counter (e.g. gaus(0) + gaus(3) )
         // need to start for a position
         Int_t openingParenthesisPos = (closingBracketPos == kNPOS) ? openingBracketPos : closingBracketPos + 1;
         bool defaultCounter = (openingParenthesisPos > formula.Length() || formula[openingParenthesisPos] != '(');

         // Int_t openingParenthesisPos = formula.Index('(',funPos);
         // Bool_t defaultCounter = (openingParenthesisPos == kNPOS);
         Int_t counter;
         if (defaultCounter) {
            counter = 0;
         } else {
            // Check whether this is just a number in parentheses. If not, leave
            // it to `HandleFunctionArguments` to be parsed

            TRegexp counterPattern("([0-9]+)");
            Ssiz_t *len = new Ssiz_t();
            if (counterPattern.Index(formula, len, openingParenthesisPos) == -1) {
               funPos = formula.Index(funName, funPos + 1);
               continue;
            } else {
               counter =
                  TString(formula(openingParenthesisPos + 1, formula.Index(')', funPos) - openingParenthesisPos - 1))
                     .Atoi();
            }
         }
         // std::cout << "openingParenthesisPos  " << openingParenthesisPos << " counter is " << counter <<  std::endl;

         TString body = (isNormalized ? it->second.second : it->second.first);
         if (isNormalized && body == "") {
            Error("PreprocessFormula", "%d dimension function %s has no normalized form.", it->first.second,
                  funName.Data());
            break;
         }
         for (int i = 0; i < body.Length(); ++i) {
            if (body[i] == '{') {
               // replace {Vn} with variable names
               i += 2; // skip '{' and 'V'
               Int_t num = TString(body(i, body.Index('}', i) - i)).Atoi();
               TString variable = variables[num];
               TString pattern = TString::Format("{V%d}", num);
               i -= 2; // restore original position
               body.Replace(i, pattern.Length(), variable, variable.Length());
               i += variable.Length() - 1; // update i to reflect change in body string
            } else if (body[i] == '[') {
               // update parameter counters in case of many functions (e.g. gaus(0)+gaus(3) )
               Int_t tmp = i;
               while (tmp < body.Length() && body[tmp] != ']') {
                  tmp++;
               }
               Int_t num = TString(body(i + 1, tmp - 1 - i)).Atoi();
               num += counter;
               TString replacement = TString::Format("%d", num);

               body.Replace(i + 1, tmp - 1 - i, replacement, replacement.Length());
               i += replacement.Length() + 1;
            }
         }
         TString pattern;
         if (defaultCounter && defaultVariables) {
            pattern = TString::Format("%s%s", funName.Data(), (isNormalized ? "n" : ""));
         }
         if (!defaultCounter && defaultVariables) {
            pattern = TString::Format("%s%s(%d)", funName.Data(), (isNormalized ? "n" : ""), counter);
         }
         if (defaultCounter && !defaultVariables) {
            pattern = TString::Format("%s%s[%s]", funName.Data(), (isNormalized ? "n" : ""), varList.Data());
         }
         if (!defaultCounter && !defaultVariables) {
            pattern =
               TString::Format("%s%s[%s](%d)", funName.Data(), (isNormalized ? "n" : ""), varList.Data(), counter);
         }
         TString replacement = body;

         // set the number (only in case a function exists without anything else
         if (fNumber == 0 && formula.Length() <= (pattern.Length() - funPos) + 1) { // leave 1 extra
            fNumber = functionsNumbers[funName] + 10 * (dim - 1);
         }

         // std::cout << " replace " << pattern << " with " << replacement << std::endl;

         formula.Replace(funPos, pattern.Length(), replacement, replacement.Length());

         funPos = formula.Index(funName);
      }
      // std::cout << " End loop of " << funName << " formula is now " << formula << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
///   Handling parameter ranges, in the form of [1..5]
void TFormula::HandleParamRanges(TString &formula)
{
   TRegexp rangePattern("\\[[0-9]+\\.\\.[0-9]+\\]");
   Ssiz_t *len = new Ssiz_t();
   int matchIdx = 0;
   while ((matchIdx = rangePattern.Index(formula, len, matchIdx)) != -1) {
      int startIdx = matchIdx + 1;
      int endIdx = formula.Index("..", startIdx) + 2; // +2 for ".."
      int startCnt = TString(formula(startIdx, formula.Length())).Atoi();
      int endCnt = TString(formula(endIdx, formula.Length())).Atoi();

      if (endCnt <= startCnt)
         Error("HandleParamRanges", "End parameter (%d) <= start parameter (%d) in parameter range", endCnt, startCnt);

      TString newString = "[";
      for (int cnt = startCnt; cnt < endCnt; cnt++)
         newString += TString::Format("%d],[", cnt);
      newString += TString::Format("%d]", endCnt);

      // std::cout << "newString generated by HandleParamRanges is " << newString << std::endl;
      formula.Replace(matchIdx, formula.Index("]", matchIdx) + 1 - matchIdx, newString);

      matchIdx += newString.Length();
   }

   // std::cout << "final formula is now " << formula << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
///   Handling user functions (and parametrized functions)
///   to take variables and optionally parameters as arguments
void TFormula::HandleFunctionArguments(TString &formula)
{
   // std::cout << "calling `HandleFunctionArguments` on " << formula << std::endl;

   // Define parametrized functions, in case we need to use them
   std::map<std::pair<TString, Int_t>, std::pair<TString, TString>> parFunctions;
   FillParametrizedFunctions(parFunctions);

   // loop through characters
   for (Int_t i = 0; i < formula.Length(); ++i) {
      // List of things to ignore (copied from `TFormula::ExtractFunctors`)

      // ignore things that start with square brackets
      if (formula[i] == '[') {
         while (formula[i] != ']')
            i++;
         continue;
      }
      // ignore strings
      if (formula[i] == '\"') {
         do
            i++;
         while (formula[i] != '\"');
         continue;
      }
      // ignore numbers (scientific notation)
      if (IsScientificNotation(formula, i))
         continue;
      // ignore x in hexadecimal number
      if (IsHexadecimal(formula, i)) {
         while (!IsOperator(formula[i]) && i < formula.Length())
            i++;
         continue;
      }

      // investigate possible start of function name
      if (isalpha(formula[i]) && !IsOperator(formula[i])) {
         // std::cout << "character : " << i << " " << formula[i] << " is not an operator and is alpha" << std::endl;

         int j; // index to end of name
         for (j = i; j < formula.Length() && IsFunctionNameChar(formula[j]); j++)
            ;
         TString name = (TString)formula(i, j - i);
         // std::cout << "parsed name " << name << std::endl;

         // Count arguments (careful about parentheses depth)
         // Make list of indices where each argument is separated
         int nArguments = 1;
         int depth = 1;
         std::vector<int> argSeparators;
         argSeparators.push_back(j); // opening parenthesis
         int k;                      // index for end of closing parenthesis
         for (k = j + 1; depth >= 1 && k < formula.Length(); k++) {
            if (formula[k] == ',' && depth == 1) {
               nArguments++;
               argSeparators.push_back(k);
            } else if (formula[k] == '(')
               depth++;
            else if (formula[k] == ')')
               depth--;
         }
         argSeparators.push_back(k - 1); // closing parenthesis

         // retrieve `f` (code copied from ExtractFunctors)
         TObject *obj = 0;
         {
            R__LOCKGUARD(gROOTMutex);
            obj = gROOT->GetListOfFunctions()->FindObject(name);
         }
         TFormula *f = dynamic_cast<TFormula *>(obj);
         if (!f) {
            // maybe object is a TF1
            TF1 *f1 = dynamic_cast<TF1 *>(obj);
            if (f1)
               f = f1->GetFormula();
         }
         // `f` should be found by now, if it was a user-defined function.
         // The other possibility we need to consider is that this is a
         // parametrized function (else case below)

         bool nameRecognized = (f != NULL);

         // Get ndim, npar, and replacementFormula of function
         int ndim;
         int npar;
         TString replacementFormula;
         if (f) {
            ndim = f->GetNdim();
            npar = f->GetNpar();
            replacementFormula = f->GetExpFormula();
         } else {
            // otherwise, try to match default parametrized functions

            for (auto keyval : parFunctions) {
               // (name, ndim)
               pair<TString, Int_t> name_ndim = keyval.first;
               // (formula without normalization, formula with normalization)
               pair<TString, TString> formulaPair = keyval.second;

               // match names like gaus, gausn, breitwigner
               if (name == name_ndim.first)
                  replacementFormula = formulaPair.first;
               else if (name == name_ndim.first + "n" && formulaPair.second != "")
                  replacementFormula = formulaPair.second;
               else
                  continue;

               // set ndim
               ndim = name_ndim.second;

               // go through replacementFormula to find the number of parameters
               npar = 0;
               int idx = 0;
               while ((idx = replacementFormula.Index('[', idx)) != kNPOS) {
                  npar = max(npar, 1 + TString(replacementFormula(idx + 1, replacementFormula.Length())).Atoi());
                  idx = replacementFormula.Index(']', idx);
                  if (idx == kNPOS)
                     Error("HandleFunctionArguments", "Square brackets not matching in formula %s",
                           (const char *)replacementFormula);
               }
               // npar should be set correctly now

               // break if number of arguments is good (note: `gaus`, has two
               // definitions with different numbers of arguments, but it works
               // out so that it should be unambiguous)
               if (nArguments == ndim + npar || nArguments == npar || nArguments == ndim) {
                  nameRecognized = true;
                  break;
               }
            }
         }
         if (nameRecognized && ndim > 4)
            Error("HandleFunctionArguments", "Number of dimensions %d greater than 4. Cannot parse formula.", ndim);

         // if we have "recognizedName(...", then apply substitutions
         if (nameRecognized && j < formula.Length() && formula[j] == '(') {
            // std::cout << "naive replacement formula: " << replacementFormula << std::endl;
            // std::cout << "formula: " << formula << std::endl;

            // map to rename each argument in `replacementFormula`
            map<TString, TString> argSubstitutions;

            const char *defaultVariableNames[] = {"x", "y", "z", "t"};

            // check nArguments and add to argSubstitutions map as appropriate
            bool canReplace = false;
            if (nArguments == ndim + npar) {
               // loop through all variables and parameters, filling in argSubstitutions
               for (int argNr = 0; argNr < nArguments; argNr++) {

                  // Get new name (for either variable or parameter)
                  TString newName =
                     TString(formula(argSeparators[argNr] + 1, argSeparators[argNr + 1] - argSeparators[argNr] - 1));
                  PreProcessFormula(newName); // so that nesting works

                  // Get old name(s)
                  // and add to argSubstitutions map as appropriate
                  if (argNr < ndim) { // variable
                     TString oldName = (f) ? TString::Format("x[%d]", argNr) : TString::Format("{V%d}", argNr);
                     argSubstitutions[oldName] = newName;

                     if (f)
                        argSubstitutions[defaultVariableNames[argNr]] = newName;

                  } else { // parameter
                     int parNr = argNr - ndim;
                     TString oldName =
                        (f) ? TString::Format("[%s]", f->GetParName(parNr)) : TString::Format("[%d]", parNr);
                     argSubstitutions[oldName] = newName;

                     // If the name stays the same, keep the old value of the parameter
                     if (f && oldName == newName)
                        DoAddParameter(f->GetParName(parNr), f->GetParameter(parNr), false);
                  }
               }

               canReplace = true;
            } else if (nArguments == npar) {
               // Try to assume variables are implicit (need all arguments to be
               // parameters)

               // loop to check if all arguments are parameters
               bool varsImplicit = true;
               for (int argNr = 0; argNr < nArguments && varsImplicit; argNr++) {
                  int openIdx = argSeparators[argNr] + 1;
                  int closeIdx = argSeparators[argNr + 1] - 1;

                  // check brackets on either end
                  if (formula[openIdx] != '[' || formula[closeIdx] != ']' || closeIdx <= openIdx + 1)
                     varsImplicit = false;

                  // check that the middle is a single function-name
                  for (int idx = openIdx + 1; idx < closeIdx && varsImplicit; idx++)
                     if (!IsFunctionNameChar(formula[idx]))
                        varsImplicit = false;

                  if (!varsImplicit)
                     Warning("HandleFunctionArguments",
                             "Argument %d is not a parameter. Cannot assume variables are implicit.", argNr);
               }

               // loop to replace parameter names
               if (varsImplicit) {
                  // if parametrized function, still need to replace parameter names
                  if (!f) {
                     for (int dim = 0; dim < ndim; dim++) {
                        argSubstitutions[TString::Format("{V%d}", dim)] = defaultVariableNames[dim];
                     }
                  }

                  for (int argNr = 0; argNr < nArguments; argNr++) {
                     TString oldName =
                        (f) ? TString::Format("[%s]", f->GetParName(argNr)) : TString::Format("[%d]", argNr);
                     TString newName =
                        TString(formula(argSeparators[argNr] + 1, argSeparators[argNr + 1] - argSeparators[argNr] - 1));

                     // preprocess the formula so that nesting works
                     PreProcessFormula(newName);
                     argSubstitutions[oldName] = newName;

                     // If the name stays the same, keep the old value of the parameter
                     if (f && oldName == newName)
                        DoAddParameter(f->GetParName(argNr), f->GetParameter(argNr), false);
                  }

                  canReplace = true;
               }
            }
            if (!canReplace && nArguments == ndim) {
               // Treat parameters as implicit

               // loop to replace variable names
               for (int argNr = 0; argNr < nArguments; argNr++) {
                  TString oldName = (f) ? TString::Format("x[%d]", argNr) : TString::Format("{V%d}", argNr);
                  TString newName =
                     TString(formula(argSeparators[argNr] + 1, argSeparators[argNr + 1] - argSeparators[argNr] - 1));

                  // preprocess so nesting works
                  PreProcessFormula(newName);
                  argSubstitutions[oldName] = newName;

                  if (f) // x, y, z are not used in parametrized function definitions
                     argSubstitutions[defaultVariableNames[argNr]] = newName;
               }

               if (f) {
                  // keep old values of the parameters
                  for (int parNr = 0; parNr < npar; parNr++)
                     DoAddParameter(f->GetParName(parNr), f->GetParameter(parNr), false);
               }

               canReplace = true;
            }

            if (canReplace)
               ReplaceAllNames(replacementFormula, argSubstitutions);
            // std::cout << "after replacement, replacementFormula is " << replacementFormula << std::endl;

            if (canReplace) {
               // std::cout << "about to replace position " << i << " length " << k-i << " in formula : " << formula <<
               // std::endl;
               formula.Replace(i, k - i, replacementFormula);
               i += replacementFormula.Length() - 1; // skip to end of replacement
               // std::cout << "new formula is : " << formula << std::endl;
            } else {
               Warning("HandleFunctionArguments", "Unable to make replacement. Number of parameters doesn't work : "
                                                  "%d arguments, %d dimensions, %d parameters",
                       nArguments, ndim, npar);
               i = j;
            }

         } else {
            i = j; // skip to end of candidate "name"
         }
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
///    Handling exponentiation
///    Can handle multiple carets, eg.
///    2^3^4 will be treated like 2^(3^4)

void TFormula::HandleExponentiation(TString &formula)
{
   Int_t caretPos = formula.Last('^');
   while(caretPos != kNPOS  && !IsAParameterName(formula,caretPos) )
      {

      TString right, left;
      Int_t temp = caretPos;
      temp--;
      // get the expression in ( ) which has the operator^ applied
      if (formula[temp] == ')') {
         Int_t depth = 1;
         temp--;
         while (depth != 0 && temp > 0) {
            if (formula[temp] == ')')
               depth++;
            if (formula[temp] == '(')
               depth--;
            temp--;
         }
         if (depth == 0)
            temp++;
      }
      // this in case of someting like sin(x+2)^2
      do {
         temp--; // go down one
         // handle scientific notation cases (1.e-2 ^ 3 )
         if (temp >= 2 && IsScientificNotation(formula, temp - 1))
            temp -= 3;
      } while (temp >= 0 && !IsOperator(formula[temp]) && !IsBracket(formula[temp]));

      assert(temp + 1 >= 0);
      Int_t leftPos = temp + 1;
      left = formula(leftPos, caretPos - leftPos);
      // std::cout << "left to replace is " << left << std::endl;

      // look now at the expression after the ^ operator
      temp = caretPos;
      temp++;
      if (temp >= formula.Length()) {
         Error("HandleExponentiation", "Invalid position of operator ^");
         return;
         }
         if (formula[temp] == '(') {
            Int_t depth = 1;
            temp++;
            while (depth != 0 && temp < formula.Length()) {
               if (formula[temp] == ')')
                  depth--;
               if (formula[temp] == '(')
                  depth++;
               temp++;
            }
            temp--;
         } else {
            // handle case  first character is operator - or + continue
            if (formula[temp] == '-' || formula[temp] == '+')
               temp++;
            // handle cases x^-2 or x^+2
            // need to handle also cases x^sin(x+y)
            Int_t depth = 0;
            // stop right expression if is an operator or if is a ")" from a zero depth
            while (temp < formula.Length() && ((depth > 0) || !IsOperator(formula[temp]))) {
               temp++;
               // handle scientific notation cases (1.e-2 ^ 3 )
               if (temp >= 2 && IsScientificNotation(formula, temp))
                  temp += 2;
               // for internal parenthesis
               if (temp < formula.Length() && formula[temp] == '(')
                  depth++;
               if (temp < formula.Length() && formula[temp] == ')') {
                  if (depth > 0)
                     depth--;
                  else
                     break; // case of end of a previously started expression e.g. sin(x^2)
               }
            }
         }
         right = formula(caretPos + 1, (temp - 1) - caretPos);
         // std::cout << "right to replace is " << right << std::endl;

         TString pattern = TString::Format("%s^%s", left.Data(), right.Data());
         TString replacement = TString::Format("pow(%s,%s)", left.Data(), right.Data());

         // std::cout << "pattern : " << pattern << std::endl;
         // std::cout << "replacement : " << replacement << std::endl;
         formula.Replace(leftPos, pattern.Length(), replacement, replacement.Length());

         caretPos = formula.Last('^');
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle linear functions defined with the operator ++.

void TFormula::HandleLinear(TString &formula)
{
   // Handle Linear functions identified with "@" operator
   Int_t linPos = formula.Index("@");
   if (linPos == kNPOS ) return;  // function is not linear
   Int_t nofLinParts = formula.CountChar((int)'@');
   assert(nofLinParts > 0);
   fLinearParts.reserve(nofLinParts + 1);
   Int_t Nlinear = 0;
   bool first = true;
   while(linPos != kNPOS && !IsAParameterName(formula,linPos))
      {
      SetBit(kLinear, 1);
      // analyze left part only the first time
      Int_t temp = 0;
      TString left;
      if (first) {
         temp = linPos - 1;
         while (temp >= 0 && formula[temp] != '@') {
            temp--;
         }
         left = formula(temp + 1, linPos - (temp + 1));
      }
      temp = linPos + 1;
      while (temp < formula.Length() && formula[temp] != '@') {
         temp++;
      }
      TString right = formula(linPos + 1, temp - (linPos + 1));

      TString pattern =
         (first) ? TString::Format("%s@%s", left.Data(), right.Data()) : TString::Format("@%s", right.Data());
      TString replacement =
         (first) ? TString::Format("([%d]*(%s))+([%d]*(%s))", Nlinear, left.Data(), Nlinear + 1, right.Data())
                 : TString::Format("+([%d]*(%s))", Nlinear, right.Data());
      Nlinear += (first) ? 2 : 1;

      formula.ReplaceAll(pattern, replacement);
      if (first) {
         TFormula *lin1 = new TFormula("__linear1", left, false);
         fLinearParts.push_back(lin1);
      }
      TFormula *lin2 = new TFormula("__linear2", right, false);
      fLinearParts.push_back(lin2);

      linPos = formula.Index("@");
      first = false;
   }
}

////////////////////////////////////////////////////////////////////////////////
///    Preprocessing of formula
///    Replace all ** by ^, and removes spaces.
///    Handle also parametrized functions like polN,gaus,expo,landau
///    and exponentiation.
///    Similar functionality should be added here.

void TFormula::PreProcessFormula(TString &formula)
{
   formula.ReplaceAll("**","^");
   formula.ReplaceAll("++","@");  // for linear functions
   formula.ReplaceAll(" ","");
   HandlePolN(formula);
   HandleParametrizedFunctions(formula);
   HandleParamRanges(formula);
   HandleFunctionArguments(formula);
   HandleExponentiation(formula);
   // "++" wil be dealt with Handle Linear
   HandleLinear(formula);
   // special case for "--" and "++"
   // ("++" needs to be written with whitespace that is removed before but then we re-add it again
   formula.ReplaceAll("--","- -");
   formula.ReplaceAll("++","+ +");
}

////////////////////////////////////////////////////////////////////////////////
/// prepare the formula to be executed
/// normally is called with fFormula

Bool_t TFormula::PrepareFormula(TString &formula)
{
   fFuncs.clear();
   fReadyToExecute = false;
   ExtractFunctors(formula);

   // update the expression with the new formula
   fFormula = formula;
   // save formula to parse variable and parameters for Cling
   fClingInput = formula;
   // replace all { and }
   fFormula.ReplaceAll("{","");
   fFormula.ReplaceAll("}","");

   //std::cout << "functors are extracted formula is " << std::endl;
   //std::cout << fFormula << std::endl << std::endl;

   fFuncs.sort();
   fFuncs.unique();

   // use inputFormula for Cling
   ProcessFormula(fClingInput);

   // for pre-defined functions (need after processing)
   if (fNumber != 0) SetPredefinedParamNames();

   return fReadyToExecute && fClingInitialized;
}

////////////////////////////////////////////////////////////////////////////////
///    Extracts functors from formula, and put them in fFuncs.
///    Simple grammar:
///  -  <function>  := name(arg1,arg2...)
///  -  <variable>  := name
///  -  <parameter> := [number]
///  -  <name>      := String containing lower and upper letters, numbers, underscores
///  -  <number>    := Integer number
///    Operators are omitted.

void TFormula::ExtractFunctors(TString &formula)
{
   // std::cout << "Commencing ExtractFunctors on " << formula << std::endl;

   TString name = "";
   TString body = "";
   // printf("formula is : %s \n",formula.Data() );
   for(Int_t i = 0 ; i < formula.Length(); ++i )
      {

      // std::cout << "loop on character : " << i << " " << formula[i] << std::endl;
      // case of parameters
      if (formula[i] == '[') {
         Int_t tmp = i;
         i++;
         TString param = "";
         while (formula[i] != ']' && i < formula.Length()) {
            param.Append(formula[i++]);
         }
         i++;
         // rename parameter name XX to pXX
         // std::cout << "examine parameters " << param << std::endl;
         int paramIndex = -1;
         if (param.IsDigit()) {
            paramIndex = param.Atoi();
            param.Insert(0, 'p'); // needed for the replacement
            if (paramIndex >= fNpar || fParams.find(param) == fParams.end()) {
               // add all parameters up to given index found
               for (int idx = 0; idx <= paramIndex; ++idx) {
                  TString pname = TString::Format("p%d", idx);
                  if (fParams.find(pname) == fParams.end())
                     DoAddParameter(pname, 0, false);
               }
            }
         } else {
            // handle whitespace characters in parname
            param.ReplaceAll("\\s", " ");

            // only add if parameter does not already exist, because maybe
            // `HandleFunctionArguments` already assigned a default value to the
            // parameter
            if (fParams.find(param) == fParams.end() || GetParNumber(param) < 0 ||
                (unsigned)GetParNumber(param) >= fClingParameters.size()) {
               // std::cout << "Setting parameter " << param << " to 0" << std::endl;
               DoAddParameter(param, 0, false);
            }
         }
         TString replacement = TString::Format("{[%s]}", param.Data());
         formula.Replace(tmp, i - tmp, replacement, replacement.Length());
         fFuncs.push_back(TFormulaFunction(param));
         // printf("found parameter %s \n",param.Data() );
         continue;
            }
            // case of strings
            if (formula[i] == '\"') {
               // look for next instance of "\"
               do {
                  i++;
               } while (formula[i] != '\"');
            }
            // case of e or E for numbers in exponential notaton (e.g. 2.2e-3)
            if (IsScientificNotation(formula, i))
               continue;
            // case of x for hexadecimal numbers
            if (IsHexadecimal(formula, i)) {
               // find position of operator
               // do not check cases if character is not only a to f, but accept anything
               while (!IsOperator(formula[i]) && i < formula.Length()) {
                  i++;
               }
               continue;
         }

         // std::cout << "investigating character : " << i << " " << formula[i] << " of formula " << formula <<
         // std::endl;
         // look for variable and function names. They  start in C++ with alphanumeric characters
         if (isalpha(formula[i]) &&
             !IsOperator(formula[i])) // not really needed to check if operator (if isalpha is not an operator)
         {
            // std::cout << "character : " << i << " " << formula[i] << " is not an operator and is alpha " <<
            // std::endl;

            while (IsFunctionNameChar(formula[i]) && i < formula.Length()) {
               // need special case for separting operator  ":" from scope operator "::"
               if (formula[i] == ':' && ((i + 1) < formula.Length())) {
                  if (formula[i + 1] == ':') {
                     // case of :: (scopeOperator)
                     name.Append("::");
                     i += 2;
                     continue;
                  } else
                     break;
               }

               name.Append(formula[i++]);
            }
            // printf(" build a name %s \n",name.Data() );
            if (formula[i] == '(') {
               i++;
               if (formula[i] == ')') {
                  fFuncs.push_back(TFormulaFunction(name, body, 0));
                  name = body = "";
                  continue;
               }
               Int_t depth = 1;
               Int_t args = 1; // we will miss first argument
               while (depth != 0 && i < formula.Length()) {
                  switch (formula[i]) {
                  case '(': depth++; break;
                  case ')': depth--; break;
                  case ',':
                     if (depth == 1)
                        args++;
                     break;
                  }
                  if (depth != 0) // we don't want last ')' inside body
                  {
                     body.Append(formula[i++]);
                  }
               }
               Int_t originalBodyLen = body.Length();
               ExtractFunctors(body);
               formula.Replace(i - originalBodyLen, originalBodyLen, body, body.Length());
               i += body.Length() - originalBodyLen;
               fFuncs.push_back(TFormulaFunction(name, body, args));
            } else {

               // std::cout << "check if character : " << i << " " << formula[i] << " from name " << name << "  is a
               // function " << std::endl;

               // check if function is provided by gROOT
               TObject *obj = 0;
               {
                  R__LOCKGUARD(gROOTMutex);
                  obj = gROOT->GetListOfFunctions()->FindObject(name);
               }
               TFormula *f = dynamic_cast<TFormula *>(obj);
               if (!f) {
                  // maybe object is a TF1
                  TF1 *f1 = dynamic_cast<TF1 *>(obj);
                  if (f1)
                     f = f1->GetFormula();
               }
               if (f) {
                  // Replacing user formula the old way (as opposed to 'HandleFunctionArguments')
                  // Note this is only for replacing functions that do
                  // not specify variables and/or parameters in brackets
                  // (the other case is done by `HandleFunctionArguments`)

                  TString replacementFormula = f->GetExpFormula();

                  // analyze expression string
                  // std::cout << "formula to replace for " << f->GetName() << " is " << replacementFormula <<
                  // std::endl;
                  PreProcessFormula(replacementFormula);
                  // we need to define different parameters if we use the unnamed default parameters ([0])
                  // I need to replace all the terms in the functor for backward compatibility of the case
                  // f1("[0]*x") f2("[0]*x") f1+f2 - it is weird but it is better to support
                  // std::cout << "current number of parameter is " << fNpar << std::endl;
                  int nparOffset = 0;
                  // if (fParams.find("0") != fParams.end() ) {
                  // do in any case if parameters are existing
                  std::vector<TString> newNames;
                  if (fNpar > 0) {
                     nparOffset = fNpar;
                     newNames.resize(f->GetNpar());
                     // start from higher number to avoid overlap
                     for (int jpar = f->GetNpar() - 1; jpar >= 0; --jpar) {
                        // parameters name have a "p" added in front
                        TString pj = TString(f->GetParName(jpar));
                        if (pj[0] == 'p' && TString(pj(1, pj.Length())).IsDigit()) {
                           TString oldName = TString::Format("[%s]", f->GetParName(jpar));
                           TString newName = TString::Format("[p%d]", nparOffset + jpar);
                           // std::cout << "replace - parameter " << f->GetParName(jpar) << " with " <<  newName <<
                           // std::endl;
                           replacementFormula.ReplaceAll(oldName, newName);
                           newNames[jpar] = newName;
                        } else
                           newNames[jpar] = f->GetParName(jpar);
                     }
                     // std::cout << "after replacing params " << replacementFormula << std::endl;
                  }
                  ExtractFunctors(replacementFormula);
                  // std::cout << "after re-extracting functors " << replacementFormula << std::endl;

                  // set parameter value from replacement formula
                  for (int jpar = 0; jpar < f->GetNpar(); ++jpar) {
                     if (nparOffset > 0) {
                        // parameter have an offset- so take this into accound
                        assert((int)newNames.size() == f->GetNpar());
                        SetParameter(newNames[jpar], f->GetParameter(jpar));
                     } else
                        // names are the same between current formula and replaced one
                        SetParameter(f->GetParName(jpar), f->GetParameter(jpar));
                  }
                  // need to add parenthesis at begin and end of replacementFormula
                  replacementFormula.Insert(0, '(');
                  replacementFormula.Insert(replacementFormula.Length(), ')');
                  formula.Replace(i - name.Length(), name.Length(), replacementFormula, replacementFormula.Length());
                  // move forward the index i of the main loop
                  i += replacementFormula.Length() - name.Length();

                  // we have extracted all the functor for "fname"
                  // std::cout << "We have extracted all the functors for fname" << std::endl;
                  // std::cout << " i = " << i << " f[i] = " << formula[i] << " - " << formula << std::endl;
                  name = "";

                  continue;
               }

               // add now functor in
               TString replacement = TString::Format("{%s}", name.Data());
               formula.Replace(i - name.Length(), name.Length(), replacement, replacement.Length());
               i += 2;
               fFuncs.push_back(TFormulaFunction(name));
            }
            }
            name = body = "";
      }
}

////////////////////////////////////////////////////////////////////////////////
///    Iterates through functors in fFuncs and performs the appropriate action.
///    If functor has 0 arguments (has only name) can be:
///     - variable
///       * will be replaced with x[num], where x is an array containing value of this variable under num.
///     - pre-defined formula
///       * will be replaced with formulas body
///     - constant
///       * will be replaced with constant value
///     - parameter
///       * will be replaced with p[num], where p is an array containing value of this parameter under num.
///    If has arguments it can be :
///     - function shortcut, eg. sin
///       * will be replaced with fullname of function, eg. sin -> TMath::Sin
///     - function from cling environment, eg. TMath::BreitWigner(x,y,z)
///       * first check if function exists, and has same number of arguments, then accept it and set as found.
///    If all functors after iteration are matched with corresponding action,
///    it inputs C++ code of formula into cling, and sets flag that formula is ready to evaluate.

void TFormula::ProcessFormula(TString &formula)
{
   // std::cout << "Begin: formula is " << formula << " list of functors " << fFuncs.size() << std::endl;

   for (list<TFormulaFunction>::iterator funcsIt = fFuncs.begin(); funcsIt != fFuncs.end(); ++funcsIt) {
      TFormulaFunction &fun = *funcsIt;

      // std::cout << "fun is " << fun.GetName() << std::endl;

      if (fun.fFound)
         continue;
      if (fun.IsFuncCall()) {
         map<TString, TString>::iterator it = fFunctionsShortcuts.find(fun.GetName());
         if (it != fFunctionsShortcuts.end()) {
            TString shortcut = it->first;
            TString full = it->second;
            // std::cout << " functor " << fun.GetName() << " found - replace " <<  shortcut << " with " << full << " in
            // " << formula << std::endl;
            // replace all functors
            Ssiz_t index = formula.Index(shortcut, 0);
            while (index != kNPOS) {
               // check that function is not in a namespace and is not in other characters
               // std::cout << "analyzing " << shortcut << " in " << formula << std::endl;
               Ssiz_t i2 = index + shortcut.Length();
               if ((index > 0) && (isalpha(formula[index - 1]) || formula[index - 1] == ':')) {
                  index = formula.Index(shortcut, i2);
                  continue;
               }
               if (i2 < formula.Length() && formula[i2] != '(') {
                  index = formula.Index(shortcut, i2);
                  continue;
               }
               // now replace the string
               formula.Replace(index, shortcut.Length(), full);
               Ssiz_t inext = index + full.Length();
               index = formula.Index(shortcut, inext);
               fun.fFound = true;
            }
         }
         if (fun.fName.Contains("::")) // add support for nested namespaces
         {
            // look for last occurence of "::"
            std::string name(fun.fName);
            size_t index = name.rfind("::");
            assert(index != std::string::npos);
            TString className = fun.fName(0, fun.fName(0, index).Length());
            TString functionName = fun.fName(index + 2, fun.fName.Length());

            Bool_t silent = true;
            TClass *tclass = TClass::GetClass(className, silent);
            // std::cout << "looking for class " << className << std::endl;
            const TList *methodList = tclass->GetListOfAllPublicMethods();
            TIter next(methodList);
            TMethod *p;
            while ((p = (TMethod *)next())) {
               if (strcmp(p->GetName(), functionName.Data()) == 0 &&
                   (fun.GetNargs() <= p->GetNargs() && fun.GetNargs() >= p->GetNargs() - p->GetNargsOpt())) {
                  fun.fFound = true;
                  break;
               }
            }
         }
         if (!fun.fFound) {
            // try to look into all the global functions in gROOT
            TFunction *f;
            {
               R__LOCKGUARD(gROOTMutex);
               f = (TFunction *)gROOT->GetListOfGlobalFunctions(true)->FindObject(fun.fName);
            }
            // if found a function with matching arguments
            if (f && fun.GetNargs() <= f->GetNargs() && fun.GetNargs() >= f->GetNargs() - f->GetNargsOpt()) {
               fun.fFound = true;
            }
         }

         if (!fun.fFound) {
            // ignore not found functions
            if (gDebug)
               Info("TFormula", "Could not find %s function with %d argument(s)", fun.GetName(), fun.GetNargs());
            fun.fFound = false;
         }
      } else {
         TFormula *old = 0;
         {
            R__LOCKGUARD(gROOTMutex);
            old = (TFormula *)gROOT->GetListOfFunctions()->FindObject(gNamePrefix + fun.fName);
         }
         if (old) {
            // we should not go here (this analysis is done before in ExtractFunctors)
            assert(false);
            fun.fFound = true;
            TString pattern = TString::Format("{%s}", fun.GetName());
            TString replacement = old->GetExpFormula();
            PreProcessFormula(replacement);
            ExtractFunctors(replacement);
            formula.ReplaceAll(pattern, replacement);
            continue;
         }
         // looking for default variables defined in fVars

         map<TString, TFormulaVariable>::iterator varsIt = fVars.find(fun.GetName());
         if (varsIt != fVars.end()) {

            TString name = (*varsIt).second.GetName();
            Double_t value = (*varsIt).second.fValue;

            AddVariable(name, value); // this set the cling variable
            if (!fVars[name].fFound) {

               fVars[name].fFound = true;
               int varDim = (*varsIt).second.fArrayPos; // variable dimenions (0 for x, 1 for y, 2, for z)
               if (varDim >= fNdim) {
                  fNdim = varDim + 1;

                  // we need to be sure that all other variables are added with position less
                  for (auto &v : fVars) {
                     if (v.second.fArrayPos < varDim && !v.second.fFound) {
                        AddVariable(v.first, v.second.fValue);
                        v.second.fFound = true;
                     }
                  }
               }
            }
            // remove the "{.. }" added around the variable
            TString pattern = TString::Format("{%s}", name.Data());
            TString replacement = TString::Format("x[%d]", (*varsIt).second.fArrayPos);
            formula.ReplaceAll(pattern, replacement);

            // std::cout << "Found an observable for " << fun.GetName()  << std::endl;

            fun.fFound = true;
            continue;
         }
         // check for observables defined as x[0],x[1],....
         // maybe could use a regular expression here
         // only in case match with defined variables is not successfull
         TString funname = fun.GetName();
         if (funname.Contains("x[") && funname.Contains("]")) {
            TString sdigit = funname(2, funname.Index("]"));
            int digit = sdigit.Atoi();
            if (digit >= fNdim) {
               fNdim = digit + 1;
               // we need to add the variables in fVars all of them before x[n]
               for (int j = 0; j < fNdim; ++j) {
                  TString vname = TString::Format("x[%d]", j);
                  if (fVars.find(vname) == fVars.end()) {
                     fVars[vname] = TFormulaVariable(vname, 0, j);
                     fVars[vname].fFound = true;
                     AddVariable(vname, 0.);
                  }
               }
            }
            // std::cout << "Found matching observable for " << funname  << std::endl;
            fun.fFound = true;
            // remove the "{.. }" added around the variable
            TString pattern = TString::Format("{%s}", funname.Data());
            formula.ReplaceAll(pattern, funname);
            continue;
         }
         //}

         auto paramsIt = fParams.find(fun.GetName());
         if (paramsIt != fParams.end()) {
            // TString name = (*paramsIt).second.GetName();
            TString pattern = TString::Format("{[%s]}", fun.GetName());
            // std::cout << "pattern is " << pattern << std::endl;
            if (formula.Index(pattern) != kNPOS) {
               // TString replacement = TString::Format("p[%d]",(*paramsIt).second.fArrayPos);
               TString replacement = TString::Format("p[%d]", (*paramsIt).second);
               // std::cout << "replace pattern  " << pattern << " with " << replacement << std::endl;
               formula.ReplaceAll(pattern, replacement);
            }
            fun.fFound = true;
            continue;
         } else {
            // std::cout << "functor  " << fun.GetName() << " is not a parameter " << std::endl;
         }

         // looking for constants (needs to be done after looking at the parameters)
         map<TString, Double_t>::iterator constIt = fConsts.find(fun.GetName());
         if (constIt != fConsts.end()) {
            TString pattern = TString::Format("{%s}", fun.GetName());
            TString value = TString::Format("%lf", (*constIt).second);
            formula.ReplaceAll(pattern, value);
            fun.fFound = true;
            // std::cout << "constant with name " << fun.GetName() << " is found " << std::endl;
            continue;
         }

         fun.fFound = false;
      }
      }
   //std::cout << "End: formula is " << formula << std::endl;

   // ignore case of functors have been matched - try to pass it to Cling
   if(!fReadyToExecute)
      {
      fReadyToExecute = true;
      Bool_t hasVariables = (fNdim > 0);
      Bool_t hasParameters = (fNpar > 0);
      if (!hasParameters) {
         fAllParametersSetted = true;
      }
      // assume a function without variables is always 1-dimensional
      if (hasParameters && !hasVariables) {
         fNdim = 1;
         AddVariable("x", 0);
         hasVariables = true;
      }
      Bool_t hasBoth = hasVariables && hasParameters;
      Bool_t inputIntoCling = (formula.Length() > 0);
      if (inputIntoCling) {

         // save copy of inputFormula in a std::strig for the unordered map
         // and also formula is same as FClingInput typically and it will be modified
         std::string inputFormula = std::string(formula);

         // valid input formula - try to put into Cling
         TString argumentsPrototype = TString::Format("%s%s%s", (hasVariables ? "Double_t *x" : ""),
                                                      (hasBoth ? "," : ""), (hasParameters ? "Double_t *p" : ""));

         // set the name for Cling using the hash_function
         fClingName = gNamePrefix;

         // check if formula exist already in the map
         R__LOCKGUARD(gROOTMutex);

         auto funcit = gClingFunctions.find(inputFormula);

         if (funcit != gClingFunctions.end()) {
            fFuncPtr = (TInterpreter::CallFuncIFacePtr_t::Generic_t)funcit->second;
            fClingInitialized = true;
            inputIntoCling = false;
         }

         // set the cling name using hash of the static formulae map
         auto hasher = gClingFunctions.hash_function();
         fClingName = TString::Format("%s__id%zu", gNamePrefix.Data(), hasher(inputFormula));

         fClingInput = TString::Format("Double_t %s(%s){ return %s ; }", fClingName.Data(), argumentsPrototype.Data(),
                                       inputFormula.c_str());

         // this is not needed (maybe can be re-added in case of recompilation of identical expressions
         // // check in case of a change if need to re-initialize
         // if (fClingInitialized) {
         //    if (oldClingInput == fClingInput)
         //       inputIntoCling = false;
         //    else
         //       fClingInitialized = false;
         // }

         if (inputIntoCling) {
            InputFormulaIntoCling();
            if (fClingInitialized) {
               // if Cling has been succesfully initialized
               // dave function ptr in the static map
               R__LOCKGUARD(gROOTMutex);
               gClingFunctions.insert(std::make_pair(inputFormula, (void *)fFuncPtr));
            }

         } else {
            fAllParametersSetted = true;
            fClingInitialized = true;
         }
         }
      }


   // IN case of a Cling Error check components wich are not found in Cling
   // check that all formula components arematched otherwise emit an error
   if (!fClingInitialized) {
      Bool_t allFunctorsMatched = true;
      for(list<TFormulaFunction>::iterator it = fFuncs.begin(); it != fFuncs.end(); it++)
         {
         if (!it->fFound) {
            allFunctorsMatched = false;
            if (it->GetNargs() == 0)
               Error("ProcessFormula", "\"%s\" has not been matched in the formula expression", it->GetName());
            else
               Error("ProcessFormula", "Could not find %s function with %d argument(s)", it->GetName(), it->GetNargs());
         }
         }
      if (!allFunctorsMatched) {
         Error("ProcessFormula","Formula \"%s\" is invalid !", GetExpFormula().Data() );
         fReadyToExecute = false;
      }
   }

   // clean up un-used default variables in case formula is valid
   if (fClingInitialized && fReadyToExecute) {
      auto itvar = fVars.begin();
      do {
         if (!itvar->second.fFound) {
            // std::cout << "Erase variable " << itvar->first << std::endl;
            itvar = fVars.erase(itvar);
         } else
            itvar++;
         }
      while( itvar != fVars.end() );
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Fill map with parametrized functions

void TFormula::FillParametrizedFunctions(map<pair<TString, Int_t>, pair<TString, TString>> &functions)
{
   // map< pair<TString,Int_t> ,pair<TString,TString> > functions;
   functions.insert(
      make_pair(make_pair("gaus", 1), make_pair("[0]*exp(-0.5*(({V0}-[1])/[2])*(({V0}-[1])/[2]))",
                                                "[0]*exp(-0.5*(({V0}-[1])/[2])*(({V0}-[1])/[2]))/(sqrt(2*pi)*[2])")));
   functions.insert(make_pair(make_pair("landau", 1), make_pair("[0]*TMath::Landau({V0},[1],[2],false)",
                                                                "[0]*TMath::Landau({V0},[1],[2],true)")));
   functions.insert(make_pair(make_pair("expo", 1), make_pair("exp([0]+[1]*{V0})", "")));
   functions.insert(
      make_pair(make_pair("crystalball", 1), make_pair("[0]*ROOT::Math::crystalball_function({V0},[3],[4],[2],[1])",
                                                       "[0]*ROOT::Math::crystalball_pdf({V0},[3],[4],[2],[1])")));
   functions.insert(
      make_pair(make_pair("breitwigner", 1), make_pair("[0]*ROOT::Math::breitwigner_pdf({V0},[2],[1])",
                                                       "[0]*ROOT::Math::breitwigner_pdf({V0},[2],[4],[1])")));
   // chebyshev polynomial
   functions.insert(make_pair(make_pair("cheb0", 1), make_pair("ROOT::Math::Chebyshev0({V0},[0])", "")));
   functions.insert(make_pair(make_pair("cheb1", 1), make_pair("ROOT::Math::Chebyshev1({V0},[0],[1])", "")));
   functions.insert(make_pair(make_pair("cheb2", 1), make_pair("ROOT::Math::Chebyshev2({V0},[0],[1],[2])", "")));
   functions.insert(make_pair(make_pair("cheb3", 1), make_pair("ROOT::Math::Chebyshev3({V0},[0],[1],[2],[3])", "")));
   functions.insert(
      make_pair(make_pair("cheb4", 1), make_pair("ROOT::Math::Chebyshev4({V0},[0],[1],[2],[3],[4])", "")));
   functions.insert(
      make_pair(make_pair("cheb5", 1), make_pair("ROOT::Math::Chebyshev5({V0},[0],[1],[2],[3],[4],[5])", "")));
   functions.insert(
      make_pair(make_pair("cheb6", 1), make_pair("ROOT::Math::Chebyshev6({V0},[0],[1],[2],[3],[4],[5],[6])", "")));
   functions.insert(
      make_pair(make_pair("cheb7", 1), make_pair("ROOT::Math::Chebyshev7({V0},[0],[1],[2],[3],[4],[5],[6],[7])", "")));
   functions.insert(make_pair(make_pair("cheb8", 1),
                              make_pair("ROOT::Math::Chebyshev8({V0},[0],[1],[2],[3],[4],[5],[6],[7],[8])", "")));
   functions.insert(make_pair(make_pair("cheb9", 1),
                              make_pair("ROOT::Math::Chebyshev9({V0},[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])", "")));
   functions.insert(
      make_pair(make_pair("cheb10", 1),
                make_pair("ROOT::Math::Chebyshev10({V0},[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10])", "")));
   // 2-dimensional functions
   functions.insert(
      make_pair(make_pair("gaus", 2), make_pair("[0]*exp(-0.5*(({V0}-[1])/[2])^2 - 0.5*(({V1}-[3])/[4])^2)", "")));
   functions.insert(
      make_pair(make_pair("landau", 2),
                make_pair("[0]*TMath::Landau({V0},[1],[2],false)*TMath::Landau({V1},[3],[4],false)", "")));
   functions.insert(make_pair(make_pair("expo", 2), make_pair("exp([0]+[1]*{V0})", "exp([0]+[1]*{V0}+[2]*{V1})")));
   // gaussian with correlations
   functions.insert(
      make_pair(make_pair("bigaus", 2), make_pair("[0]*ROOT::Math::bigaussian_pdf({V0},{V1},[2],[4],[5],[1],[3])",
                                                  "[0]*ROOT::Math::bigaussian_pdf({V0},{V1},[2],[4],[5],[1],[3])")));
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter names only in case of pre-defined functions.

void TFormula::SetPredefinedParamNames() {

   if (fNumber == 0) return;

   if (fNumber == 100) { // Gaussian
      SetParName(0,"Constant");
      SetParName(1,"Mean");
      SetParName(2,"Sigma");
      return;
   }
   if (fNumber == 110) {
      SetParName(0,"Constant");
      SetParName(1,"MeanX");
      SetParName(2,"SigmaX");
      SetParName(3,"MeanY");
      SetParName(4,"SigmaY");
      return;
   }
   if (fNumber == 112) {  // bigaus
      SetParName(0,"Constant");
      SetParName(1,"MeanX");
      SetParName(2,"SigmaX");
      SetParName(3,"MeanY");
      SetParName(4,"SigmaY");
      SetParName(5,"Rho");
      return;
   }
   if (fNumber == 200) { // exponential
      SetParName(0,"Constant");
      SetParName(1,"Slope");
      return;
   }
   if (fNumber == 400) { // landau
      SetParName(0,"Constant");
      SetParName(1,"MPV");
      SetParName(2,"Sigma");
      return;
   }
   if (fNumber == 500) { // crystal-ball
      SetParName(0,"Constant");
      SetParName(1,"Mean");
      SetParName(2,"Sigma");
      SetParName(3,"Alpha");
      SetParName(4,"N");
      return;
   }
   if (fNumber == 600) { // breit-wigner
      SetParName(0,"Constant");
      SetParName(1,"Mean");
      SetParName(2,"Gamma");
      return;
   }
   // if formula is a polynomial (or chebyshev), set parameter names
   // not needed anymore (p0 is assigned by default)
   // if (fNumber == (300+fNpar-1) ) {
   //    for (int i = 0; i < fNpar; i++) SetParName(i,TString::Format("p%d",i));
   //    return;
   // }

   // // general case if parameters are digits (XX) change to pXX
   // auto paramMap = fParams;  // need to copy the map because SetParName is going to modify it
   // for ( auto & p : paramMap) {
   //    if (p.first.IsDigit() )
   //        SetParName(p.second,TString::Format("p%s",p.first.Data()));
   // }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return linear part.

const TObject* TFormula::GetLinearPart(Int_t i) const
{
   if (!fLinearParts.empty()) {
      int n = fLinearParts.size();
      if (i < 0 || i >= n ) {
         Error("GetLinearPart","Formula %s has only %d linear parts - requested %d",GetName(),n,i);
         return nullptr;
      }
      return fLinearParts[i];
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
///    Adds variable to known variables, and reprocess formula.

void TFormula::AddVariable(const TString &name, double value)
{
   if(fVars.find(name) != fVars.end() )
      {
      TFormulaVariable &var = fVars[name];
      var.fValue = value;

      // If the position is not defined in the Cling vectors, make space for it
      // but normally is variable is defined in fVars a slot should be also present in fClingVariables
      if (var.fArrayPos < 0) {
         var.fArrayPos = fVars.size();
      }
      if (var.fArrayPos >= (int)fClingVariables.size()) {
         fClingVariables.resize(var.fArrayPos + 1);
      }
      fClingVariables[var.fArrayPos] = value;
   } else {
      TFormulaVariable var(name, value, fVars.size());
      fVars[name] = var;
      fClingVariables.push_back(value);
      if (!fFormula.IsNull()) {
         // printf("process formula again - %s \n",fClingInput.Data() );
         ProcessFormula(fClingInput);
      }
      }

}

////////////////////////////////////////////////////////////////////////////////
///    Adds multiple variables.
///    First argument is an array of pairs<TString,Double>, where
///    first argument is name of variable,
///    second argument represents value.
///    size - number of variables passed in first argument

void TFormula::AddVariables(const TString *vars, const Int_t size)
{
   Bool_t anyNewVar = false;
   for(Int_t i = 0 ; i < size; ++i)
      {

      const TString &vname = vars[i];

      TFormulaVariable &var = fVars[vname];
      if (var.fArrayPos < 0) {

         var.fName = vname;
         var.fArrayPos = fVars.size();
         anyNewVar = true;
         var.fValue = 0;
         if (var.fArrayPos >= (int)fClingVariables.capacity()) {
            Int_t multiplier = 2;
            if (fFuncs.size() > 100) {
               multiplier = TMath::Floor(TMath::Log10(fFuncs.size()) * 10);
            }
            fClingVariables.reserve(multiplier * fClingVariables.capacity());
         }
         fClingVariables.push_back(0.0);
            }
            // else
            // {
            //    var.fValue = v.second;
            //    fClingVariables[var.fArrayPos] = v.second;
            // }
      }
      if (anyNewVar && !fFormula.IsNull()) {
         ProcessFormula(fClingInput);
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the name of the formula. We need to allow the list of function to
/// properly handle the hashes.

void TFormula::SetName(const char* name)
{
   if (IsReservedName(name)) {
      Error("SetName", "The name \'%s\' is reserved as a TFormula variable name.\n"
                       "\tThis function will not be renamed.",
            name);
   } else {
      // Here we need to remove and re-add to keep the hashes consistent with
      // the underlying names.
      auto listOfFunctions = gROOT->GetListOfFunctions();
      TObject* thisAsFunctionInList = nullptr;
      R__LOCKGUARD(gROOTMutex);
      if (listOfFunctions){
         thisAsFunctionInList = listOfFunctions->FindObject(this);
         if (thisAsFunctionInList) listOfFunctions->Remove(thisAsFunctionInList);
      }
      TNamed::SetName(name);
      if (thisAsFunctionInList) listOfFunctions->Add(thisAsFunctionInList);
   }
}

////////////////////////////////////////////////////////////////////////////////
///
///    Sets multiple variables.
///    First argument is an array of pairs<TString,Double>, where
///    first argument is name of variable,
///    second argument represents value.
///    size - number of variables passed in first argument

void TFormula::SetVariables(const pair<TString,Double_t> *vars, const Int_t size)
{
   for(Int_t i = 0; i < size; ++i)
      {
      pair<TString, Double_t> v = vars[i];
      if (fVars.find(v.first) != fVars.end()) {
         fVars[v.first].fValue = v.second;
         fClingVariables[fVars[v.first].fArrayPos] = v.second;
      } else {
         Error("SetVariables", "Variable %s is not defined.", v.first.Data());
      }
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns variable value.

Double_t TFormula::GetVariable(const char *name) const
{
   TString sname(name);
   if (fVars.find(sname) == fVars.end()) {
      Error("GetVariable", "Variable %s is not defined.", sname.Data());
      return -1;
   }
   return fVars.find(sname)->second.fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns variable number (positon in array) given its name.

Int_t TFormula::GetVarNumber(const char *name) const
{
   TString sname(name);
   if (fVars.find(sname) == fVars.end()) {
      Error("GetVarNumber", "Variable %s is not defined.", sname.Data());
      return -1;
   }
   return fVars.find(sname)->second.fArrayPos;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns variable name given its position in the array.

TString TFormula::GetVarName(Int_t ivar) const
{
   if (ivar < 0 || ivar >= fNdim) return "";

   // need to loop on the map to find corresponding variable
   for ( auto & v : fVars) {
      if (v.second.fArrayPos == ivar) return v.first;
   }
   Error("GetVarName","Variable with index %d not found !!",ivar);
   //return TString::Format("x%d",ivar);
   return TString();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets variable value.

void TFormula::SetVariable(const TString &name, Double_t value)
{
   if (fVars.find(name) == fVars.end()) {
      Error("SetVariable", "Variable %s is not defined.", name.Data());
      return;
   }
   fVars[name].fValue = value;
   fClingVariables[fVars[name].fArrayPos] = value;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds parameter to known parameters.
/// User should use SetParameter, because parameters are added during initialization part,
/// and after that adding new will be pointless.

void TFormula::DoAddParameter(const TString &name, Double_t value, Bool_t processFormula)
{
   //std::cout << "adding parameter " << name << std::endl;

   // if parameter is already defined in fParams - just set the new value
   if(fParams.find(name) != fParams.end() )
      {
      int ipos = fParams[name];
      // TFormulaVariable & par = fParams[name];
      // par.fValue = value;
      if (ipos < 0) {
         ipos = fParams.size();
         fParams[name] = ipos;
      }
      //
      if (ipos >= (int)fClingParameters.size()) {
         if (ipos >= (int)fClingParameters.capacity())
            fClingParameters.reserve(TMath::Max(int(fParams.size()), ipos + 1));
         fClingParameters.insert(fClingParameters.end(), ipos + 1 - fClingParameters.size(), 0.0);
      }
      fClingParameters[ipos] = value;
   } else {
      // new parameter defined
      fNpar++;
      // TFormulaVariable(name,value,fParams.size());
      int pos = fParams.size();
      // fParams.insert(std::make_pair<TString,TFormulaVariable>(name,TFormulaVariable(name,value,pos)));
      auto ret = fParams.insert(std::make_pair(name, pos));
      // map returns a std::pair<iterator, bool>
      // use map the order for defult position of parameters in the vector
      // (i.e use the alphabetic order)
      if (ret.second) {
         // a new element is inserted
         if (ret.first == fParams.begin())
            pos = 0;
         else {
            auto previous = (ret.first);
            --previous;
            pos = previous->second + 1;
         }

         if (pos < (int)fClingParameters.size())
            fClingParameters.insert(fClingParameters.begin() + pos, value);
         else {
            // this should not happen
            if (pos > (int)fClingParameters.size())
               Warning("inserting parameter %s at pos %d when vector size is  %d \n", name.Data(), pos,
                       (int)fClingParameters.size());

            if (pos >= (int)fClingParameters.capacity())
               fClingParameters.reserve(TMath::Max(int(fParams.size()), pos + 1));
            fClingParameters.insert(fClingParameters.end(), pos + 1 - fClingParameters.size(), 0.0);
            fClingParameters[pos] = value;
         }

         // need to adjust all other positions
         for (auto it = ret.first; it != fParams.end(); ++it) {
            it->second = pos;
            pos++;
         }

         // for (auto & p : fParams)
         //    std::cout << "Parameter " << p.first << " position " << p.second << " value " <<
         //    fClingParameters[p.second] << std::endl;
         // printf("inserted parameters size params %d size cling %d \n",fParams.size(), fClingParameters.size() );
      }
      if (processFormula) {
         // replace first in input parameter name with [name]
         fClingInput.ReplaceAll(name, TString::Format("[%s]", name.Data()));
         ProcessFormula(fClingInput);
      }
      }

}

////////////////////////////////////////////////////////////////////////////////
/// Return parameter index given a name (return -1 for not existing parameters)
/// non need to print an error

Int_t TFormula::GetParNumber(const char * name) const {
   auto it = fParams.find(name);
   if (it == fParams.end()) {
      return -1;
   }
   return it->second;

}

////////////////////////////////////////////////////////////////////////////////
/// Returns parameter value given by string.

Double_t TFormula::GetParameter(const char * name) const
{
   int i = GetParNumber(name);
   if (i  == -1) {
      Error("GetParameter","Parameter %s is not defined.",name);
      return TMath::QuietNaN();
   }

   return GetParameter( GetParNumber(name) );
}

////////////////////////////////////////////////////////////////////////////////
/// Return parameter value given by integer.

Double_t TFormula::GetParameter(Int_t param) const
{
   //TString name = TString::Format("%d",param);
   if(param >=0 && param < (int) fClingParameters.size())
      return fClingParameters[param];
   Error("GetParameter","wrong index used - use GetParameter(name)");
   return TMath::QuietNaN();
}

////////////////////////////////////////////////////////////////////////////////
///    Return parameter name given by integer.

const char * TFormula::GetParName(Int_t ipar) const
{
   if (ipar < 0 || ipar >= fNpar) return "";

   // need to loop on the map to find corresponding parameter
   for ( auto & p : fParams) {
      if (p.second == ipar) return p.first.Data();
   }
   Error("GetParName","Parameter with index %d not found !!",ipar);
   //return TString::Format("p%d",ipar);
   return TString();
}

////////////////////////////////////////////////////////////////////////////////
Double_t* TFormula::GetParameters() const
{
   if(!fClingParameters.empty())
      return const_cast<Double_t*>(&fClingParameters[0]);
   return 0;
}

void TFormula::GetParameters(Double_t *params) const
{
   for (Int_t i = 0; i < fNpar; ++i) {
      if (Int_t(fClingParameters.size()) > i)
         params[i] = fClingParameters[i];
      else
         params[i] = -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets parameter value.

void TFormula::SetParameter(const char *name, Double_t value)
{
   SetParameter( GetParNumber(name), value);

   // do we need this ???
#ifdef OLDPARAMS
   if (fParams.find(name) == fParams.end()) {
      Error("SetParameter", "Parameter %s is not defined.", name.Data());
      return;
   }
   fParams[name].fValue = value;
   fParams[name].fFound = true;
   fClingParameters[fParams[name].fArrayPos] = value;
   fAllParametersSetted = true;
   for(map<TString,TFormulaVariable>::iterator it = fParams.begin(); it != fParams.end(); it++)
      {
      if (!it->second.fFound) {
         fAllParametersSetted = false;
         break;
      }
      }
#endif
}

#ifdef OLDPARAMS

////////////////////////////////////////////////////////////////////////////////
/// Set multiple parameters.
/// First argument is an array of pairs<TString,Double>, where
/// first argument is name of parameter,
/// second argument represents value.
/// size - number of params passed in first argument

void TFormula::SetParameters(const pair<TString,Double_t> *params,const Int_t size)
{
   for(Int_t i = 0 ; i < size ; ++i)
      {
      pair<TString, Double_t> p = params[i];
      if (fParams.find(p.first) == fParams.end()) {
         Error("SetParameters", "Parameter %s is not defined", p.first.Data());
         continue;
      }
      fParams[p.first].fValue = p.second;
      fParams[p.first].fFound = true;
      fClingParameters[fParams[p.first].fArrayPos] = p.second;
      }
   fAllParametersSetted = true;
   for(map<TString,TFormulaVariable>::iterator it = fParams.begin(); it != fParams.end(); it++)
      {
      if (!it->second.fFound) {
         fAllParametersSetted = false;
         break;
      }
      }
}
#endif

////////////////////////////////////////////////////////////////////////////////
void TFormula::DoSetParameters(const Double_t *params, Int_t size)
{
   if(!params || size < 0 || size > fNpar) return;
   // reset vector of cling parameters
   if (size != (int) fClingParameters.size() ) {
      Warning("SetParameters","size is not same of cling parameter size %d - %d",size,int(fClingParameters.size()) );
      for (Int_t i = 0; i < size; ++i) {
         TString name = TString::Format("%d", i);
         SetParameter(name, params[i]);
      }
      return;
   }
   fAllParametersSetted = true;
   std::copy(params, params+size, fClingParameters.begin() );
}

////////////////////////////////////////////////////////////////////////////////
/// Set a vector of parameters value.
/// Order in the vector is by default the aphabetic order given to the parameters
/// apart if the users has defined explicitly the parameter names

void TFormula::SetParameters(const Double_t *params)
{
   DoSetParameters(params,fNpar);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a list of parameters.
/// The order is by default the aphabetic order given to the parameters
/// apart if the users has defined explicitly the parameter names

void TFormula::SetParameters(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5, Double_t p6,
                             Double_t p7, Double_t p8, Double_t p9, Double_t p10)
{
   if(fNpar >= 1) SetParameter(0,p0);
   if(fNpar >= 2) SetParameter(1,p1);
   if(fNpar >= 3) SetParameter(2,p2);
   if(fNpar >= 4) SetParameter(3,p3);
   if(fNpar >= 5) SetParameter(4,p4);
   if(fNpar >= 6) SetParameter(5,p5);
   if(fNpar >= 7) SetParameter(6,p6);
   if(fNpar >= 8) SetParameter(7,p7);
   if(fNpar >= 9) SetParameter(8,p8);
   if(fNpar >= 10) SetParameter(9,p9);
   if(fNpar >= 11) SetParameter(10,p10);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a parameter given a parameter index
/// The parameter index is by default the aphabetic order given to the parameters
/// apart if the users has defined explicitly the parameter names

void TFormula::SetParameter(Int_t param, Double_t value)
{
   if (param < 0 || param >= fNpar) return;
   assert(int(fClingParameters.size()) == fNpar);
   fClingParameters[param] = value;
   // TString name = TString::Format("%d",param);
   // SetParameter(name,value);
}

////////////////////////////////////////////////////////////////////////////////
void TFormula::SetParNames(const char *name0, const char *name1, const char *name2, const char *name3,
                           const char *name4, const char *name5, const char *name6, const char *name7,
                           const char *name8, const char *name9, const char *name10)
{
   if (fNpar >= 1)
      SetParName(0, name0);
   if (fNpar >= 2)
      SetParName(1, name1);
   if (fNpar >= 3)
      SetParName(2, name2);
   if (fNpar >= 4)
      SetParName(3, name3);
   if (fNpar >= 5)
      SetParName(4, name4);
   if (fNpar >= 6)
      SetParName(5, name5);
   if (fNpar >= 7)
      SetParName(6, name6);
   if (fNpar >= 8)
      SetParName(7, name7);
   if (fNpar >= 9)
      SetParName(8, name8);
   if (fNpar >= 10)
      SetParName(9, name9);
   if (fNpar >= 11)
      SetParName(10, name10);
}

////////////////////////////////////////////////////////////////////////////////
void TFormula::SetParName(Int_t ipar, const char * name)
{

   if (ipar < 0 || ipar > fNpar) {
      Error("SetParName","Wrong Parameter index %d ",ipar);
      return;
   }
   TString oldName;
   // find parameter with given index
   for ( auto &it : fParams) {
      if (it.second  == ipar) {
         oldName =  it.first;
         fParams.erase(oldName);
         fParams.insert(std::make_pair(name, ipar) );
         break;
      }
   }
   if (oldName.IsNull() ) {
      Error("SetParName","Parameter %d is not existing.",ipar);
      return;
   }

   //replace also parameter name in formula expression
   ReplaceParamName(fFormula, oldName, name);

}

////////////////////////////////////////////////////////////////////////////////
/// Replace in Formula expression the parameter name.

void TFormula::ReplaceParamName(TString & formula, const TString & oldName, const TString & name){
   if (!formula.IsNull() ) {
      bool found = false;
      for(list<TFormulaFunction>::iterator it = fFuncs.begin(); it != fFuncs.end(); ++it)
         {
         if (oldName == it->GetName()) {
            found = true;
            it->fName = name;
            break;
         }
         }
         if (!found) {
            Error("SetParName", "Parameter %s is not defined.", oldName.Data());
            return;
         }
         // change whitespace to \\s avoid problems in parsing
         TString newName = name;
         newName.ReplaceAll(" ", "\\s");
         TString pattern = TString::Format("[%s]", oldName.Data());
         TString replacement = TString::Format("[%s]", newName.Data());
         formula.ReplaceAll(pattern, replacement);
   }

}

////////////////////////////////////////////////////////////////////////////////
Double_t TFormula::EvalPar(const Double_t *x,const Double_t *params) const
{
   return DoEval(x, params);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets first 4  variables (e.g. x, y, z, t) and evaluate formula.

Double_t TFormula::Eval(Double_t x, Double_t y, Double_t z, Double_t t) const
{
   double xxx[4] = {x,y,z,t};
   return DoEval(xxx);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets first 3  variables (e.g. x, y, z) and evaluate formula.

Double_t TFormula::Eval(Double_t x, Double_t y , Double_t z) const
{
   double xxx[3] = {x,y,z};
   return DoEval(xxx);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets first 2  variables (e.g. x and y) and evaluate formula.

Double_t TFormula::Eval(Double_t x, Double_t y) const
{
   double xxx[2] = {x,y};
   return DoEval(xxx);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets first variable (e.g. x) and evaluate formula.

Double_t TFormula::Eval(Double_t x) const
{
   //double xxx[1] = {x};
   double * xxx = &x;
   return DoEval(xxx);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate formula.
/// If formula is not ready to execute(missing parameters/variables),
/// print these which are not known.
/// If parameter has default value, and has not been setted, appropriate warning is shown.

Double_t TFormula::DoEval(const double * x, const double * params) const
{
   if(!fReadyToExecute)
      {
      Error("Eval", "Formula is invalid and not ready to execute ");
      for (auto it = fFuncs.begin(); it != fFuncs.end(); ++it) {
         TFormulaFunction fun = *it;
         if (!fun.fFound) {
            printf("%s is unknown.\n", fun.GetName());
         }
      }
      return TMath::QuietNaN();
      }
   if (fLambdaPtr && TestBit(TFormula::kLambda)) {// case of lambda functions
      std::function<double(double *, double *)> & fptr = * ( (std::function<double(double *, double *)> *) fLambdaPtr);
      assert(x);
      //double * v = (x) ? const_cast<double*>(x) : const_cast<double*>(fClingVariables.data());
      double * v = const_cast<double*>(x);
      double * p = (params) ? const_cast<double*>(params) : const_cast<double*>(fClingParameters.data());
      return fptr(v, p);
   }
   // this is needed when reading from a file
   if (!fClingInitialized) {
      Error("Eval","Formula is invalid or not properly initialized - try calling TFormula::Compile");
      return TMath::QuietNaN();
#ifdef EVAL_IS_NOT_CONST
      // need to replace in cling the name of the pointer of this object
      TString oldClingName = fClingName;
      fClingName.Replace(fClingName.Index("_0x")+1,fClingName.Length(), TString::Format("%p",this) );
      fClingInput.ReplaceAll(oldClingName, fClingName);
      InputFormulaIntoCling();
#endif
   }

   Double_t result = 0;
   void* args[2];
   double * vars = (x) ? const_cast<double*>(x) : const_cast<double*>(fClingVariables.data());
   args[0] = &vars;
   if (fNpar <= 0)
      (*fFuncPtr)(0, 1, args, &result);
   else {
      double * pars = (params) ? const_cast<double*>(params) : const_cast<double*>(fClingParameters.data());
      args[1] = &pars;
      (*fFuncPtr)(0, 2, args, &result);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the expression formula.
///
///  - If option = "P" replace the parameter names with their values
///  - If option = "CLING" return the actual expression used to build the function  passed to cling
///  - If option = "CLINGP" replace in the CLING expression the parameter with their values

TString TFormula::GetExpFormula(Option_t *option) const
{
   TString opt(option);
   if (opt.IsNull() || TestBit(TFormula::kLambda) ) return fFormula;
   opt.ToUpper();

   //  if (opt.Contains("N") ) {
   //    TString formula = fFormula;
   //    ReplaceParName(formula, ....)
   // }

   if (opt.Contains("CLING") ) {
      std::string clingFunc = fClingInput.Data();
      std::size_t found = clingFunc.find("return");
      std::size_t found2 = clingFunc.rfind(";");
      if (found == std::string::npos || found2 == std::string::npos) {
         Error("GetExpFormula","Invalid Cling expression - return default formula expression");
         return fFormula;
      }
      TString clingFormula = fClingInput(found+7,found2-found-7);
      // to be implemented
      if (!opt.Contains("P")) return clingFormula;
      // replace all "p[" with "[parname"
      int i = 0;
      while (i < clingFormula.Length()-2 ) {
         // look for p[number
         if (clingFormula[i] == 'p' && clingFormula[i+1] == '[' && isdigit(clingFormula[i+2]) ) {
            int j = i+3;
            while ( isdigit(clingFormula[j]) ) { j++;}
            if (clingFormula[j] != ']') {
               Error("GetExpFormula","Parameters not found - invalid expression - return default cling formula");
               return clingFormula;
            }
            TString parNumbName = clingFormula(i+2,j-i-2);
            int parNumber = parNumbName.Atoi();
            assert(parNumber < fNpar);
            TString replacement = TString::Format("%f",GetParameter(parNumber));
            clingFormula.Replace(i,j-i+1, replacement );
            i += replacement.Length();
         }
         i++;
      }
      return clingFormula;
   }
   if (opt.Contains("P") ) {
      // replace parameter names with their values
      TString expFormula = fFormula;
      int i = 0;
      while (i < expFormula.Length()-2 ) {
         // look for [parName]
         if (expFormula[i] == '[') {
            int j = i+1;
            while ( expFormula[j] != ']' ) { j++;}
            if (expFormula[j] != ']') {
               Error("GetExpFormula","Parameter names not found - invalid expression - return default formula");
               return expFormula;
            }
            TString parName = expFormula(i+1,j-i-1);
            TString replacement = TString::Format("%g",GetParameter(parName));
            expFormula.Replace(i,j-i+1, replacement );
            i += replacement.Length();
         }
         i++;
      }
      return expFormula;
   }
   Warning("GetExpFormula","Invalid option - return defult formula expression");
   return fFormula;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the formula and its attributes.

void TFormula::Print(Option_t *option) const
{
   printf(" %20s : %s Ndim= %d, Npar= %d, Number= %d \n",GetName(),GetTitle(), fNdim,fNpar,fNumber);
   printf(" Formula expression: \n");
   printf("\t%s \n",fFormula.Data() );
   TString opt(option);
   opt.ToUpper();
   // do an evaluation as a cross-check
   //if (fReadyToExecute) Eval();

   if (opt.Contains("V") ) {
      if (fNdim > 0 && !TestBit(TFormula::kLambda)) {
         printf("List of  Variables: \n");
         assert(int(fClingVariables.size()) >= fNdim);
         for ( int ivar = 0; ivar < fNdim ; ++ivar) {
            printf("Var%4d %20s =  %10f \n",ivar,GetVarName(ivar).Data(), fClingVariables[ivar]);
         }
      }
      if (fNpar > 0) {
         printf("List of  Parameters: \n");
         if ( int(fClingParameters.size()) < fNpar)
            Error("Print","Number of stored parameters in vector %zu in map %zu is different than fNpar %d",fClingParameters.size(), fParams.size(), fNpar);
         assert(int(fClingParameters.size()) >= fNpar);
         // print with order passed to Cling function
         for ( int ipar = 0; ipar < fNpar ; ++ipar) {
            printf("Par%4d %20s =  %10f \n",ipar,GetParName(ipar), fClingParameters[ipar] );
         }
      }
      printf("Expression passed to Cling:\n");
      printf("\t%s\n",fClingInput.Data() );
   }
   if(!fReadyToExecute)
      {
      Warning("Print", "Formula is not ready to execute. Missing parameters/variables");
      for (list<TFormulaFunction>::const_iterator it = fFuncs.begin(); it != fFuncs.end(); ++it) {
         TFormulaFunction fun = *it;
         if (!fun.fFound) {
            printf("%s is unknown.\n", fun.GetName());
         }
      }
      }
      if (!fAllParametersSetted) {
         // we can skip this
         // Info("Print","Not all parameters are setted.");
         // for(map<TString,TFormulaVariable>::const_iterator it = fParams.begin(); it != fParams.end(); ++it)
         // {
         //    pair<TString,TFormulaVariable> param = *it;
         //    if(!param.second.fFound)
         //    {
         //       printf("%s has default value %lf\n",param.first.Data(),param.second.GetInitialValue());
         //    }
         // }
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TFormula::Streamer(TBuffer &b)
{
   if (b.IsReading() ) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      //std::cout << "version " << v << std::endl;
      if (v <= 8 && v > 3 && v != 6) {
         // old TFormula class
         ROOT::v5::TFormula * fold = new ROOT::v5::TFormula();
         // read old TFormula class
         fold->Streamer(b, v,  R__s, R__c, TFormula::Class());
         //std::cout << "read old tformula class " << std::endl;
         TFormula fnew(fold->GetName(), fold->GetExpFormula() );

         *this = fnew;

         //          printf("copying content in a new TFormula \n");
         SetParameters(fold->GetParameters() );
         if (!fReadyToExecute ) {
            Error("Streamer","Old formula read from file is NOT valid");
            Print("v");
         }
         delete fold;
         return;
      }
      else if (v > 8) {
         // new TFormula class
         b.ReadClassBuffer(TFormula::Class(), this, v, R__s, R__c);

         //std::cout << "reading npar = " << GetNpar() << std::endl;

         // initialize the formula
         // need to set size of fClingVariables which is transient
         //fClingVariables.resize(fNdim);

         // case of formula contains only parameters
         if (fFormula.IsNull() ) return;


         // store parameter values, names and order
         std::vector<double> parValues = fClingParameters;
         auto paramMap = fParams;
         fNpar = fParams.size();

         if (!TestBit(TFormula::kLambda) ) {

            //std::cout << "Streamer::Reading preprocess the formula " << fFormula << " ndim = " << fNdim << " npar = " << fNpar << std::endl;
            // for ( auto &p : fParams)
            //    std::cout << "parameter " << p.first << " index " << p.second << std::endl;

            fClingParameters.clear();  // need to be reset before re-initializing it

            FillDefaults();


            PreProcessFormula(fFormula);

            //std::cout << "Streamer::after pre-process the formula " << fFormula << " ndim = " << fNdim << " npar = " << fNpar << std::endl;

            PrepareFormula(fFormula);

            //std::cout << "Streamer::after prepared " << fClingInput << " ndim = " << fNdim << " npar = " << fNpar << std::endl;


            // restore parameter values
            if (fNpar != (int) parValues.size() ) {
               Error("Streamer","number of parameters computed (%d) is not same as the stored parameters (%d)",fNpar,int(parValues.size()) );
               Print("v");
            }
         }
         else {
            // case of lamda expressions
            bool ret = InitLambdaExpression(fFormula);
            if (ret) {
               fReadyToExecute  = true;
               fClingInitialized  = true;
            }
         }
         assert(fNpar == (int) parValues.size() );
         std::copy( parValues.begin(), parValues.end(), fClingParameters.begin() );
         // restore parameter names and order
         if (fParams.size() != paramMap.size() ) {
            Warning("Streamer","number of parameters list found (%zu) is not same as the stored one (%zu) - use re-created list",fParams.size(),paramMap.size()) ;
            //Print("v");
         }
         else
            //assert(fParams.size() == paramMap.size() );
            fParams = paramMap;

         // input formula into Cling
         // need to replace in cling the name of the pointer of this object
         // TString oldClingName = fClingName;
         // fClingName.Replace(fClingName.Index("_0x")+1,fClingName.Length(), TString::Format("%p",this) );
         // fClingInput.ReplaceAll(oldClingName, fClingName);
         // InputFormulaIntoCling();

         if (!TestBit(kNotGlobal)) {
            R__LOCKGUARD(gROOTMutex);
            gROOT->GetListOfFunctions()->Add(this);
         }
         if (!fReadyToExecute ) {
            Error("Streamer","Formula read from file is NOT ready to execute");
            Print("v");
         }
         //std::cout << "reading 2 npar = " << GetNpar() << std::endl;

         return;
      }
      else {
         Error("Streamer","Reading version %d is not supported",v);
         return;
      }
   }
   else {
      // case of writing
      b.WriteClassBuffer(TFormula::Class(), this);
      // std::cout << "writing npar = " << GetNpar() << std::endl;
   }
}
