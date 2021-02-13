// @(#)root/hist:$Id$
// Author: Nicolas Brun   19/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cmath>

#include "TROOT.h"
#include "TClass.h"
#include "TBuffer.h"
#include "v5/TFormula.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFunction.h"
#include "TMethodCall.h"
#include "TObjString.h"
#include "TError.h"
#include "v5/TFormulaPrimitive.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "strlcpy.h"
#include "snprintf.h"

#ifdef WIN32
#pragma optimize("",off)
#endif

static Int_t gMAXOP=1000,gMAXPAR=1000,gMAXCONST=1000;
const Int_t  gMAXSTRINGFOUND = 10;
const UInt_t kOptimizationError = BIT(19);

ClassImp(ROOT::v5::TFormula);

namespace ROOT {

   namespace v5 {

/** \class TFormula  TFormula.h "inc/v5/TFormula.h"
     \ingroup Hist
The FORMULA class (ROOT version 5)

 Example of valid expressions:

    -  `sin(x)/x`
    -  `[0]*sin(x) + [1]*exp(-[2]*x)`
    -  `x + y**2`
    -  `x^2 + y^2`
    -  `[0]*pow([1],4)`
    -  `2*pi*sqrt(x/y)`
    -  `gaus(0)*expo(3)  + ypol3(5)*x`
    -  `gausn(0)*expo(3) + ypol3(5)*x`

 In the last example above:

    gaus(0) is a substitute for `[0]*exp(-0.5*((x-[1])/[2])**2)`
       and (0) means start numbering parameters at 0

    gausn(0) is a substitute for `[0]*exp(-0.5*((x-[1])/[2])**2)/(sqrt(2*pi)*[2]))`
       and (0) means start numbering parameters at 0

    expo(3) is a substitute for `exp([3]+[4]*x)`

    pol3(5) is a substitute for `par[5]+par[6]*x+par[7]*x**2+par[8]*x**3`
        (here Pol3 stands for Polynomial of degree 3)

  TMath functions can be part of the expression, eg:

    -  `TMath::Landau(x)*sin(x)`
    -  `TMath::Erf(x)`

  Comparisons operators are also supported (&&, ||, ==, <=, >=, !)
  Examples:

        sin(x*(x<0.5 || x>1))

  If the result of a comparison is TRUE, the result is 1, otherwise 0.

  Already predefined names can be given. For example, if the formula

    TFormula old(sin(x*(x<0.5 || x>1))) one can assign a name to the formula. By default
    the name of the object = title = formula itself.
    old.SetName("old").
    then, old can be reused in a new expression.
    TFormula new("x*old") is equivalent to:
    TFormula new("x*sin(x*(x<0.5 || x>1))")

  Up to 4 dimensions are supported (indicated by x, y, z, t)
  An expression may have 0 parameters or a list of parameters
  indicated by the sequence [par_number]

  A graph showing the logic to compile and analyze a formula
  is shown in TFormula::Compile and TFormula::Analyze.
  Once a formula has been compiled, it can be evaluated for a given
  set of parameters. see graph in TFormula::EvalPar.

  This class is the base class for the function classes TF1,TF2 and TF3.
  It is also used by the ntuple selection mechanism TNtupleFormula.

  In version 7 of TFormula, the usage of fOper has been changed
  to improve the performance of TFormula::EvalPar.
  Conceptually, fOper was changed from a simple array of Int_t
  to an array of composite values.
  For example a 'ylandau(5)' operation used to be encoded as 4105;
  it is now encoded as (klandau >> kTFOperShift) + 5
  Any class inheriting from TFormula and using directly fOper (which
  is now a private data member), needs to be updated to take this
  in consideration.  The member functions recommended to set and
  access fOper are:  SetAction, GetAction, GetActionParam
  For more performant access to the information, see the implementation
  TFormula::EvalPar

  ### CHANGING DEFAULT SETTINGS

  When creating complex formula , it may be necessary to increase
  some default parameters. see static function TFormula::SetMaxima

  ### WHY TFormula CANNOT ACCEPT A CLASS MEMBER FUNCTION ?

  This is a frequently asked question.
  C++ is a strongly typed language. There is no way for TFormula (without
  recompiling this class) to know about all possible user defined data types.
  This also apply to the case of a static class function.
  Because TMath is a special and frequent case, TFormula is aware
  of all TMath functions.
*/

////////////////////////////////////////////////////////////////////////////////
/// Formula default constructor.

TFormula::TFormula(): TNamed()
{
   fNdim   = 0;
   fNpar   = 0;
   fNoper  = 0;
   fNconst = 0;
   fNumber = 0;
   fExpr   = 0;
   fOper   = 0;
   fConst  = 0;
   fParams = 0;
   fNstring= 0;
   fNames  = 0;
   fNval   = 0;
   //
   //MI change
   fNOperOptimized = 0;
   fExprOptimized  = 0;
   fOperOptimized  = 0;
   fOperOffset     = 0;
   fPredefined     = 0;
   fOptimal = (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalParOld;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Formula constructor.

TFormula::TFormula(const char *name,const char *expression) :
   TNamed(name,expression)
{
   fNdim   = 0;
   fNpar   = 0;
   fNoper  = 0;
   fNconst = 0;
   fNumber = 0;
   fExpr   = 0;
   fOper   = 0;
   fConst  = 0;
   fParams = 0;
   fNstring= 0;
   fNames  = 0;
   fNval   = 0;
   //
   //MI change
   fNOperOptimized = 0;
   fExprOptimized  = 0;
   fOperOptimized  = 0;
   fOperOffset     = 0;
   fPredefined     = 0;
   fOptimal = (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalParOld;

   if (!expression || !*expression) {
      Error("TFormula", "expression may not be 0 or have 0 length");
      return;
   }

   //eliminate blanks in expression
   Int_t i,j,nch;
   nch = strlen(expression);
   char *expr = new char[nch+1];
   j = 0;
   for (i=0;i<nch;i++) {
      if (expression[i] == ' ') continue;
      if (i > 0 && (expression[i] == '*') && (expression[i-1] == '*')) {
         expr[j-1] = '^';
         continue;
      }
      expr[j] = expression[i]; j++;
   }
   expr[j] = 0;
   Bool_t gausNorm   = kFALSE;
   Bool_t landauNorm = kFALSE;
   Bool_t linear = kFALSE;

   if (j) {
      TString chaine = expr;
      //special case for functions for linear fitting
      if (chaine.Contains("++"))
         linear = kTRUE;
      // special case for normalized gaus
      if (chaine.Contains("gausn")) {
         gausNorm = kTRUE;
         TString tmp = chaine;
         tmp.ReplaceAll("gausn","");
         tmp.ReplaceAll("landaun","");
         if ( tmp.Contains("gaus")  )
            Warning("TFormula","Cannot use both gaus and gausn - gaus will be treated as gausn");
         if ( tmp.Contains("landau")  )
            Warning("TFormula","Cannot use both gausn and landau - landau will be treated as landaun");
      }
      // special case for normalized landau
      if (chaine.Contains("landaun")) {
         landauNorm = kTRUE;
         TString tmp = chaine;
         tmp.ReplaceAll("landaun","");
         tmp.ReplaceAll("gausn","");
         if ( tmp.Contains("gaus")  ) {
            Warning("TFormula","Cannot use both gaus and landaun - gaus will be treated as gausn");
         }
         if ( tmp.Contains("landau") )
            Warning("TFormula","Cannot use both landau and landaun - landau will be treated as landaun");
      }
      // need to the replacement here for the error message before
      if (gausNorm)
         chaine.ReplaceAll("gausn","gaus");
      if (landauNorm)
         chaine.ReplaceAll("landaun","landau");

      SetTitle(chaine.Data());
   }
   delete [] expr;

   if (linear)    SetBit(kLinear);

   if (Compile()) return;

   if (gausNorm)   SetBit(kNormalized);
   if (landauNorm) SetBit(kNormalized);

   // Store formula in linked list of formula in ROOT


   if (strcmp(name,"x")==0 || strcmp(name,"y")==0 ||
       strcmp(name,"z")==0 || strcmp(name,"t")==0 )
   {
      Error("TFormula","The name \'%s\' is reserved as a TFormula variable name.\n"
         "\tThis function will not be registered in the list of functions",name);
   } else {
      R__LOCKGUARD(gROOTMutex);
      TFormula *old = (TFormula*)gROOT->GetListOfFunctions()->FindObject(name);
      if (old) {
         gROOT->GetListOfFunctions()->Remove(old);
      }
      gROOT->GetListOfFunctions()->Add(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TFormula::TFormula(const TFormula &formula) : TNamed()
{
   fNdim   = 0;
   fNpar   = 0;
   fNoper  = 0;
   fNconst = 0;
   fNumber = 0;
   fExpr   = 0;
   fOper   = 0;
   fConst  = 0;
   fParams = 0;
   fNstring= 0;
   fNames  = 0;
   fNval   = 0;
   fNOperOptimized = 0;
   fPredefined     = 0;
   fOperOffset     = 0;
   fExprOptimized  = 0;
   fOperOptimized  = 0;
   fOptimal = (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalParOld;

   ((TFormula&)formula).TFormula::Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TFormula& TFormula::operator=(const TFormula &rhs)
{
   if (this != &rhs) {
      rhs.Copy(*this);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Formula default destructor.

TFormula::~TFormula()
{
   if (gROOT) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfFunctions()->Remove(this);
   }

   ClearFormula();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the chain as function call.
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///     TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///     MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.

Bool_t TFormula::AnalyzeFunction(TString &chaine, Int_t &err, Int_t offset)
{
   int i;

   // We have to decompose the chain is 3 potential components:
   //   namespace::functionName( args )

   Ssiz_t argStart = chaine.First('(');
   if (argStart<0) return false;

   TString functionName = chaine(0,argStart);

   // This does not support template yet (where the scope operator might be in
   // one of the template arguments
   Ssiz_t scopeEnd = functionName.Last(':');
   TString spaceName;
   if (scopeEnd>0 && functionName[scopeEnd-1]==':') {
      spaceName = functionName(0,scopeEnd-1);
      functionName.Remove(0,scopeEnd+1);
   }

   // Now we need to count and decompose the actual arguments, we could also check the type
   // of the arguments
   if (chaine[chaine.Length()-1] != ')') {
      Error("AnalyzeFunction","We thought we had a function but we dont (in %s)\n",chaine.Data());
   }

   TString args = chaine(argStart+1,chaine.Length()-2-argStart);
   TObjArray argArr;
   argArr.SetOwner(kTRUE);
   //fprintf(stderr,"args are '%s'\n",args.Data());

   Bool_t inString = false;
   int paran = 0;
   int brack = 0;
   int prevComma = 0;
   int nargs = 0;
   for(i=0; i<args.Length(); i++) {
      if (args[i]=='"') inString = !inString;
      if (inString) continue;

      Bool_t foundArg = false;
      switch(args[i]) {

         case '(': paran++; break;
         case ')': paran--; break;
         case '[': brack++; break;
         case ']': brack--; break;

         case ',': if (paran==0 && brack==0) { foundArg = true; } break;
      }
      if ((i+1)==args.Length()) {
         foundArg = true; i++;
      }
      if (foundArg) {
         TString arg = args(prevComma,i-prevComma);

         // Here we could
         //   a) check the type
         //fprintf(stderr,"found #%d arg %s\n",nargs,arg.Data());

         // We register the arg for later usage
         argArr.Add(new TObjString(arg));
         nargs++;

         prevComma = i+1;
      };
   }

   if (nargs>999) {
      err = 7;
      return false;
   }

   // Now we need to lookup the function and check its arguments.
   TClass *ns = (spaceName.Length()) ? TClass::GetClass(spaceName) : 0;
   ClassInfo_t *cinfo = 0;
   if (ns) {
      cinfo = ns->GetClassInfo();
   } else {
      cinfo = gInterpreter->ClassInfo_Factory();
   }

   // ROOT does yet have a complete TType class, but TCling does,
   // so let's use that for now.
   static TypeInfo_t *const doubletype { gInterpreter->TypeInfo_Factory("double") };

   std::vector<TypeInfo_t*> proto(nargs,doubletype);

   CallFunc_t *callfunc = gInterpreter->CallFunc_Factory();
   Long_t func_offset;
   gInterpreter->CallFunc_SetFuncProto(callfunc,cinfo,functionName,proto,false,&func_offset,ROOT::kConversionMatch);

   TMethodCall *method = new TMethodCall(ns,callfunc,func_offset);

   if (!ns) gInterpreter->ClassInfo_Delete(cinfo);
   gInterpreter->CallFunc_Delete(callfunc);

   if (method->IsValid()) {
      if (method->ReturnType() == TMethodCall::kOther) {
         /*
           Error("Compile",
               "TFormula can only call interpreted and compiled function that returns a numerical type %s returns a %s\n",
               method->GetMethodName(), method->GetMethod()->GetReturnTypeName());
         */
         err=29;

      } else {

         // Analyze the arguments
         TIter next(&argArr);
         TObjString *objstr;
         while ( (objstr=(TObjString*)next()) ) {
            Analyze(objstr->String(),err,offset);
         }

         fFunctions.Add(method);
         fExpr[fNoper] = method->GetMethod()->GetPrototype();
         SetAction(fNoper, kFunctionCall, fFunctions.GetLast()*1000 + nargs);
         fNoper++;
         return true;
      }
   }

   delete method;
   //
   // MI change - extended space of functions
   // not forward compatible change
   //
   TString cbase(chaine);
   Int_t args_paran = cbase.First("(");
   if (args_paran>0){
      cbase[args_paran]=0;
   }

   ROOT::v5::TFormulaPrimitive *prim = ROOT::v5::TFormulaPrimitive::FindFormula(cbase, args_paran>0 ? cbase.Data() + args_paran + 1 : (const char*)0);
   if (prim &&   (!IsA()->GetBaseClass("TTreeFormula"))) {
      // TO BE DONE ALSO IN TTREFORMULA - temporary fix MI
      // Analyze the arguments
      TIter next(&argArr);
      TObjString *objstr;
      while ( (objstr=(TObjString*)next()) ) {
         Analyze(objstr->String(),err,offset); if (err) return kFALSE;
      }
      if (nargs!=prim->fNArguments) {
         Error("Compile",        "%s requires %d arguments",
            prim->GetName(), prim->fNArguments);
         return kFALSE;
      }
      fExpr[fNoper] = prim->GetName();
      if (prim->fType==10){
         SetAction(fNoper, kFD1);
      }
      if (prim->fType==110){
         SetAction(fNoper, kFD2);
      }
      if (prim->fType==1110){
         SetAction(fNoper, kFD3);
      }
      if (prim->fType==-1){
         SetAction(fNoper, kFDM);
         if (fNpar<prim->fNParameters) fNpar+=prim->fNParameters;
      }

      fNoper++;
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Analyze a sub-expression in one formula.
///
///   Expressions in one formula are recursively analyzed.
///   Result of analysis is stored in the object tables.
///
///   ### Table of function codes and errors
///
/// ~~~ {.cpp}
///   * functions :
///
///     +           1                   pow          20
///     -           2                   sq           21
///     *           3                   sqrt         22
///     /           4                   strstr       23
///     %           5                   min          24
///                                     max          25
///                                     log          30
///     cos         10                  exp          31
///     sin         11                  log10        32
///     tan         12
///     acos        13                  abs          41
///     asin        14                  sign         42
///     atan        15                  int          43
///     atan2       16
///     fmod        17                  rndm         50
///
///     cosh        70                  acosh        73
///     sinh        71                  asinh        74
///     tanh        72                  atanh        75
///
///     expo       100                  gaus        110     gausn  (see note below)
///     expo(0)    100 0                gaus(0)     110 0   gausn(0)
///     expo(1)    100 1                gaus(1)     110 1   gausn(1)
///     xexpo      100 x                xgaus       110 x   xgausn
///     yexpo      101 x                ygaus       111 x   ygausn
///     zexpo      102 x                zgaus       112 x   zgausn
///     xyexpo     105 x                xygaus      115 x   xygausn
///     yexpo(5)   102 5                ygaus(5)    111 5   ygausn(5)
///     xyexpo(2)  105 2                xygaus(2)   115 2   xygausn(2)
///
///     landau      120 x   landaun (see note below)
///     landau(0)   120 0   landaun(0)
///     landau(1)   120 1   landaun(1)
///     xlandau     120 x   xlandaun
///     ylandau     121 x   ylandaun
///     zlandau     122 x   zlandaun
///     xylandau    125 x   xylandaun
///     ylandau(5)  121 5   ylandaun(5)
///     xylandau(2) 125 2   xylandaun(2)
///
///     pol0        130 x               pol1        130 1xx
///     pol0(0)     130 0               pol1(0)     130 100
///     pol0(1)     130 1               pol1(1)     130 101
///     xpol0       130 x               xpol1       130 101
///     ypol0       131 x               ypol1       131 101
///     zpol0       132 x               zpol1       132 1xx
///     ypol0(5)    131 5               ypol1(5)    131 105
///
///     pi          40
///
///     &&          60                  <            64
///     ||          61                  >            65
///     ==          62                  <=           66
///     !=          63                  =>           67
///     !           68
///     ==(string)  76                  &            78
///     !=(string)  77                  |            79
///     <<(shift)   80                  >>(shift)    81
///     ? :         82
///
///   * constants (kConstants) :
///
///    c0  141 1      c1  141 2  etc..
///
///   * strings (kStringConst):
///
///    sX  143 x
///
///   * variables (kFormulaVar) :
///
///     x    144 0      y    144 1      z    144 2      t    144 3
///
///   * parameters :
///
///     [1]        140 1
///     [2]        140 2
///     etc.
/// ~~~
///
///   ### Special cases for normalized gaussian or landau distributions
///
///   the expression "gaus" is a substitute for
///
///       [0]*exp(-0.5*((x-[1])/[2])**2)
///
///   to obtain a standard normalized gaussian, use "gausn" instead of "gaus"
///   the expression "gausn" is a substitute for
///
///       [0]*exp(-0.5*((x-[1])/[2])**2)/(sqrt(2*pi)*[2]))
///
///   WARNING: gaus and gausn are mutually exclusive in the same expression.
///
///   In the same way the expression "landau" is a substitute for
///
///       [0]*TMath::Landau(x,[1],[2],kFALSE)
///
///   to obtain a standard normalized landau, use "landaun" instead of "landau"
///   the expression "landaun" is a substitute for
///
///       [0]*TMath::Landau(x,[1],[2],kTRUE)
///
///   WARNING: landau and landaun are mutually exclusive in the same expression.
///
///   ### Boolean optimization (kBoolOptmize) :
///
///     Those pseudo operation are used to implement lazy evaluation of
///     && and ||.  When the left hand of the expression if false
///     (respectively true), the evaluation of the right is entirely skipped
///     (since it would not change the value of the expression).
///
///     &&   142 11 (one operation on right) 142 21 (2 operations on right)
///     ||   142 12 (one operation on right) 142 22 (2 operations on right)
///
///   * functions calls (kFunctionCall) :
///
///    f0 145  0  f1 145  1  etc..
///
///   ### Errors :
///
///       1   : Division By Zero
///       2   : Invalid Floating Point Operation
///       4   : Empty String
///       5   : invalid syntax
///       6   : Too many operators
///       7   : Too many parameters
///       10  : z specified but not x and y
///       11  : z and y specified but not x
///       12  : y specified but not x
///       13  : z and x specified but not y
///       20  : non integer value for parameter number
///       21  : atan2 requires two arguments
///       22  : pow requires two arguments
///       23  : degree of polynomial not specified
///       24  : Degree of polynomial must be positive
///       25  : Degree of polynomial must be less than 20
///       26  : Unknown name
///       27  : Too many constants in expression
///       28  : strstr requires two arguments
///       29  : interpreted or compiled function have to return a numerical type
///       30  : Bad numerical expression
///       31  : Part of the variable exist but some of it is not accessible or useable
///       40  : '(' is expected
///       41  : ')' is expected
///       42  : '[' is expected
///       43  : ']' is expected
///
/// \image html TFormula_analyze.png
///
///  ### Special functions
///
///  By default, the formula is assigned fNumber=0. However, the following
///  formula built with simple functions are assigned  fNumber:
///
///       "gaus"      100  (or gausn)
///       "xygaus"    110
///       "expo"      200
///       "polN"      300+N
///       "landau"    400
///       "xylandau"  410
///
///  Note that expressions like gaus(0), expo(1) will force fNumber=0
///
///  ### Warning when deriving a class from TFormula
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///       TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///       MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.

void TFormula::Analyze(const char *schain, Int_t &err, Int_t offset)
{


   Int_t valeur,find,n,i,j,k,lchain,nomb,virgule,inter,nest;
   valeur=find=n=i=j=k=lchain=nomb=virgule=inter=nest = 0;
   Int_t compt,compt2,compt3,compt4;
   Bool_t inString;
   Double_t vafConst;
   ULong_t vafConst2;
   Bool_t parenthese;
   TString s,chaine_error,chaine1ST;
   TString s1,s2,s3,ctemp;

   TString chaine = schain;
   const TFormula *oldformula;
   Int_t modulo,plus,puiss10,puiss10bis,moins,multi,divi,puiss,et,ou,petit,grand,egal,diff,peteg,grdeg,etx,oux,rshift,lshift,tercond,terelse;
   char t;
   TString slash("/"), escapedSlash("\\/");
   Int_t inter2 = 0;
   SetNumber(0);
   Int_t actionCode,actionParam;
   Int_t err_hint = 0;

   // Verify correct matching of parenthesis and remove unnecessary parenthesis.
   lchain = chaine.Length();
   //if (chaine(lchain-2,2) == "^2") chaine = "sq(" + chaine(0,lchain-2) + ")";
   parenthese = kTRUE;
   lchain = chaine.Length();
   while (parenthese && lchain>0 && err==0){
      compt  = 0;
      compt2 = 0;
      inString = false;
      lchain = chaine.Length();
      if (lchain==0) err=4;
      else {
         for (i=1; i<=lchain; ++i) {
            if (chaine(i-1,1) == "\"") inString = !inString;
            if (!inString) {
               if (chaine(i-1,1) == "[") compt2++;
               if (chaine(i-1,1) == "]") compt2--;
               if (chaine(i-1,1) == "(") compt++;
               if (chaine(i-1,1) == ")") compt--;
            }
            if (compt < 0) err = 40; // more open parentheses than close parentheses
            if (compt2< 0) err = 42; // more ] than [
            if (compt==0 && (i!=lchain || lchain==1)) parenthese = kFALSE;
            // if (lchain<3 && chaine(0,1)!="(" && chaine(lchain-1,1)!=")") parenthese = kFALSE;
         }
         if (compt > 0) err = 41; // more ( than )
         if (compt2> 0) err = 43; // more [ than ]
         if (parenthese) chaine = chaine(1,lchain-2);
      }
   } // while parenthesis

   if (lchain==0) err=4; // empty string
   modulo=plus=moins=multi=divi=puiss=et=ou=petit=grand=egal=diff=peteg=grdeg=etx=oux=rshift=lshift=tercond=terelse=0;

   // Look for simple operators

   if (err==0) {
      compt = compt2 = compt3 = compt4 = 0;puiss10=0;puiss10bis = 0;
      inString = false;
      j = lchain;
      Bool_t isdecimal = 1; // indicates whether the left part is decimal.

      for (i=1;i<=lchain; i++) {

         puiss10=puiss10bis=0;
         if (i>2) {
            t = chaine[i-3];
            isdecimal = isdecimal && (strchr("0123456789.",t)!=0);
            if (isdecimal) {
               if ( chaine[i-2] == 'e' || chaine[i-2] == 'E' ) puiss10 = 1;
            } else if ( strchr("+-/[]()&|><=!*/%^\\",t) ) {
               isdecimal = 1; // reset after delimiter
            }
         }
         if (j>2) {
            if (chaine[j-2] == 'e' || chaine[j-2] == 'E') {
               Bool_t isrightdecimal = 1;

               for(k=j-3; k>=0 && isrightdecimal; --k) {
                  t = chaine[k];
                  isrightdecimal = isrightdecimal && (strchr("0123456789.",t)!=0);
                  if (!isrightdecimal) {
                     if (strchr("+-/[]()&|><=!*/%^\\",t)!=0) {
                        puiss10bis = 1;
                     }
                  }
               }
               if (k<0 && isrightdecimal)  puiss10bis = 1;
            }
         }
         if (puiss10 && (i<=lchain)) {
            t = chaine[i];
            puiss10 = (strchr("0123456789.",t)!=0);
         }
         if (puiss10bis && (j<=lchain)) {
            t = chaine[j];
            puiss10bis = (strchr("0123456789.",t)!=0);
         }

         if (chaine(i-1,1) == "\"") inString = !inString;
         if (inString) continue;
         if (chaine(i-1,1) == "[") compt2++;
         if (chaine(i-1,1) == "]") compt2--;
         if (chaine(i-1,1) == "(") compt++;
         if (chaine(i-1,1) == ")") compt--;
         if (chaine(j-1,1) == "[") compt3++;
         if (chaine(j-1,1) == "]") compt3--;
         if (chaine(j-1,1) == "(") compt4++;
         if (chaine(j-1,1) == ")") compt4--;
         if (chaine(i-1,2)=="&&" && !inString && compt==0 && compt2==0 && et==0) {et=i;puiss=0;}
         if (chaine(i-1,2)=="||" && compt==0 && compt2==0 && ou==0) {puiss10=0; ou=i;}
         if (chaine(i-1,1)=="&"  && compt==0 && compt2==0 && etx==0) {etx=i;puiss=0;}
         if (chaine(i-1,1)=="|"  && compt==0 && compt2==0 && oux==0) {puiss10=0; oux=i;}
         if (chaine(i-1,2)==">>" && compt==0 && compt2==0 && rshift==0) {puiss10=0; rshift=i;}
         if (chaine(i-1,1)==">"  && compt==0 && compt2==0 && rshift==0 && grand==0)
            {puiss10=0; grand=i;}
         if (chaine(i-1,2)=="<<" && compt==0 && compt2==0 && lshift==0) {puiss10=0; lshift=i;}
         if (chaine(i-1,1)=="<"  && compt==0 && compt2==0 && lshift==0 && petit==0)
            {puiss10=0; petit=i;
            // Check whether or not we have a template names! (actually this can
            // only happen in TTreeFormula.
            for(int ip = i,depth=0; ip < lchain; ++ip) {
               char c = chaine(ip);
               // The characters allowed in the template parameter are alpha-numerical characters,
               // underscores, comma, <, > and scope operator.
               if (isalnum(c) || c=='_' || c==',') continue;
               if (c==':' && chaine(ip+1)==':') { ++ip; continue; }
               if (c=='<') { ++depth; continue; }
               if (c=='>') {
                  if (depth) { --depth; continue; }
                  else {
                     // We reach the end of the template parameter.
                     petit = 0;
                     i = ip+1;
                     break;
                  }
               }
               // Character not authorized within a template parameter
               break;
            }
            if (petit==0) {
               // We found a template parameter and modified i
               continue; // the for(int i ,...)
            }
         }
         if ((chaine(i-1,2)=="<=" || chaine(i-1,2)=="=<") && compt==0 && compt2==0
            && peteg==0) {peteg=i; puiss10=0; petit=0;}
         if ((chaine(i-1,2)=="=>" || chaine(i-1,2)==">=") && compt==0 && compt2==0
            && grdeg==0) {puiss10=0; grdeg=i; grand=0;}
         if (chaine(i-1,2) == "==" && compt == 0 && compt2 == 0 && egal == 0) {puiss10=0; egal=i;}
         if (chaine(i-1,2) == "!=" && compt == 0 && compt2 == 0 && diff == 0) {puiss10=0; diff=i;}
         if (i>1 && chaine(i-1,1) == "+" && compt == 0 && compt2 == 0 && puiss10==0) plus=i;
         if (chaine(j-1,1) == "-" && chaine(j-2,1) != "*" && chaine(j-2,1) != "/"
            && chaine(j-2,1)!="^" && compt3==0 && compt4==0 && moins==0 && puiss10bis==0) moins=j;
         if (chaine(i-1,1)=="%" && compt==0 && compt2==0 && modulo==0) {puiss10=0; modulo=i;}
         if (chaine(i-1,1)=="*" && compt==0 && compt2==0 && multi==0)  {puiss10=0; multi=i;}
         if (chaine(j-1,1)=="/" && chaine(j-2,1)!="\\"
            && compt4==0 && compt3==0 && divi==0)
         {
            puiss10=0; divi=j;
         }
         if (chaine(j-1)=='^' && compt4==0 && compt3==0 && puiss==0) {puiss10=0; puiss=j;}
         if (chaine(i-1)=='?' && compt == 0 && compt2 == 0 && tercond == 0) {puiss10=0; tercond=i;}
         if (chaine(i-1)==':' && tercond && compt == 0 && compt2 == 0 && terelse == 0) {
            if (i>2 && chaine(i-2)!=':' && chaine(i)!=':') {
               puiss10=0; terelse=i;
            }
         }

         j--;
      }

   // If operator found, analyze left and right part of the statement

      actionParam = 0;
      if (tercond && terelse) {
         if (tercond == 1 || terelse == lchain || tercond == (terelse-1) ) {
            err = 5;
            chaine_error = "?:";
         } else {
            // Condition
            ctemp = chaine(0,tercond-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;

            fExpr[fNoper] = "?: condition jump";
            actionCode = kJumpIf;
            actionParam = 0;
            SetAction(fNoper,actionCode, actionParam);
            Int_t optloc = fNoper++;

            // Expression executed if condition is true.
            ctemp = chaine(tercond,terelse-tercond-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            actionParam = fNoper; // We want to skip the next instruction ('else jump'), so we set the param to the current cursor and the next instruction will be skip by the ++i in the eval loop
            SetAction(optloc, actionCode, actionParam);

            fExpr[fNoper] = "?: else jump";
            actionCode = kJump;
            actionParam = 0;
            // Set jump target.
            SetAction(fNoper,actionCode, actionParam);
            optloc = fNoper++;

            // Expression executed if condition is false.
            ctemp = chaine(terelse,lchain-terelse);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            // Set jump target.
            actionParam = fNoper - 1; // We need to not skip the next instruction, so we compensate for the ++i in the eval loop
            SetAction(optloc, actionCode, actionParam);

            if (IsString(optloc-1) != IsString(fNoper-1)) {
               err = 45;
               chaine_error = "?:";
            }
         }
      } else if (ou != 0) {    //check for ||
         if (ou==1 || ou==lchain-1) {
            err=5;
            chaine_error="||";
         }
         else {
            ctemp = chaine(0,ou-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;

            fExpr[fNoper] = "|| checkpoint";
            actionCode = kBoolOptimize;
            actionParam = 2;
            SetAction(fNoper,actionCode, actionParam);
            Int_t optloc = fNoper++;

            ctemp = chaine(ou+1,lchain-ou-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "||";
            actionCode = kOr;
            SetAction(fNoper,actionCode, 0);

            SetAction( optloc, GetAction(optloc), GetActionParam(optloc) + (fNoper-optloc) * 10);
            fNoper++;
            if (!CheckOperands(optloc-1,fNoper-1,err)) return;
         }
      } else if (et!=0) {
         if (et==1 || et==lchain-1) {
            err=5;
            chaine_error="&&";
         }
         else {
            ctemp = chaine(0,et-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;

            fExpr[fNoper] = "&& checkpoint";
            actionCode = kBoolOptimize;
            actionParam = 1;
            SetAction(fNoper,actionCode,actionParam);

            Int_t optloc = fNoper++;

            ctemp = chaine(et+1,lchain-et-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "&&";
            actionCode = kAnd;
            SetAction(fNoper,actionCode,0);

            SetAction(optloc, GetAction(optloc), GetActionParam(optloc) + (fNoper-optloc) * 10);
            fNoper++;
            if (!CheckOperands(optloc-1,fNoper-1,err)) return;
         }
      } else if (oux!=0) {
         if (oux==1 || oux==lchain) {
            err=5;
            chaine_error="|";
         }
         else {
            ctemp = chaine(0,oux-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(oux,lchain-oux);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "|";
            actionCode = kBitOr;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (etx!=0) {
         if (etx==1 || etx==lchain) {
            err=5;
            chaine_error="&";
         }
         else {
            ctemp = chaine(0,etx-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(etx,lchain-etx);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "&";
            actionCode = kBitAnd;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (petit != 0) {
         if (petit==1 || petit==lchain) {
            err=5;
            chaine_error="<";
         }
         else {
            ctemp = chaine(0,petit-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(petit,lchain-petit);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "<";
            actionCode = kLess;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (grand != 0) {
         if (grand==1 || grand==lchain) {
            err=5;
            chaine_error=">";
         }
         else {
            ctemp = chaine(0,grand-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(grand,lchain-grand);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = ">";
            actionCode = kGreater;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (peteg != 0) {
         if (peteg==1 || peteg==lchain-1) {
            err=5;
            chaine_error="<=";
         }
         else {
            ctemp = chaine(0,peteg-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            ctemp = chaine(peteg+1,lchain-peteg-1);
            UInt_t leftopr = fNoper-1;
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "<=";
            actionCode = kLessThan;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (grdeg != 0) {
         if (grdeg==1 || grdeg==lchain-1) {
            err=5;
            chaine_error="=>";
         }
         else {
            ctemp = chaine(0,grdeg-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(grdeg+1,lchain-grdeg-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = ">=";
            actionCode = kGreaterThan;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else if (egal != 0) {
         if (egal==1 || egal==lchain-1) {
            err=5;
            chaine_error="==";
         }
         else {
            ctemp = chaine(0,egal-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            Int_t optloc = fNoper-1;

            ctemp = chaine(egal+1,lchain-egal-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "==";
            actionCode = kEqual;

            Bool_t isstring = IsString(fNoper-1);
            if (IsString(optloc) != isstring) {
               err = 45;
               chaine_error = "==";
            } else if (isstring) {
               actionCode = kStringEqual;
            }
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
         }
      } else if (diff != 0) {
         if (diff==1 || diff==lchain-1) {
            err=5;
            chaine_error = "!=";
         }
         else {
            ctemp = chaine(0,diff-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            Int_t optloc = fNoper-1;

            ctemp = chaine(diff+1,lchain-diff-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "!=";
            actionCode = kNotEqual;

            Bool_t isstring = IsString(fNoper-1);
            if (IsString(optloc) != isstring) {
               err = 45;
               chaine_error = "!=";
            } else if (isstring) {
               actionCode = kStringNotEqual;
            }
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
         }
      } else if (plus != 0) {
         if (plus==lchain) {
            err=5;
            chaine_error = "+";
         }
         else {
            ctemp = chaine(0,plus-1);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            UInt_t leftopr = fNoper-1;
            ctemp = chaine(plus,lchain-plus);
            Analyze(ctemp.Data(),err,offset); if (err) return;
            fExpr[fNoper] = "+";
            actionCode = kAdd;
            SetAction(fNoper,actionCode,actionParam);
            fNoper++;
            if (!CheckOperands(leftopr,fNoper-1,err)) return;
         }
      } else {
         if (moins != 0) {
            if (moins == 1) {
               ctemp = chaine(moins,lchain-moins);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               fExpr[fNoper] = "-";
               actionCode = kSignInv;
               SetAction(fNoper,actionCode,actionParam);
               ++fNoper;
               if (!CheckOperands(fNoper-1,err)) return;
            } else {
               if (moins == lchain) {
                  err=5;
                  chaine_error = "-";
               } else {
                  ctemp = chaine(0,moins-1);
                  Analyze(ctemp.Data(),err,offset); if (err) return;
                  UInt_t leftopr = fNoper-1;
                  ctemp = chaine(moins,lchain-moins);
                  Analyze(ctemp.Data(),err,offset); if (err) return;
                  fExpr[fNoper] = "-";
                  actionCode = kSubstract;
                  SetAction(fNoper,actionCode,actionParam);
                  fNoper++;
                  if (!CheckOperands(leftopr,fNoper-1,err)) return;
               }
            }
         } else if (modulo != 0) {
            if (modulo == 1 || modulo == lchain) {
               err=5;
               chaine_error="%";
            } else {
               ctemp = chaine(0,modulo-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               UInt_t leftopr = fNoper-1;
               ctemp = chaine(modulo,lchain-modulo);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               fExpr[fNoper] = "%";
               actionCode = kModulo;
               SetAction(fNoper,actionCode,actionParam);
               fNoper++;
               if (!CheckOperands(leftopr,fNoper-1,err)) return;
            }
         } else if (rshift != 0) {
            if (rshift == 1 || rshift == lchain) {
               err=5;
               chaine_error=">>";
            } else {
               ctemp = chaine(0,rshift-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               UInt_t leftopr = fNoper-1;
               ctemp = chaine(rshift+1,lchain-rshift-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               fExpr[fNoper] = ">>";
               actionCode = kRightShift;
               SetAction(fNoper,actionCode,actionParam);
               fNoper++;
               if (!CheckOperands(leftopr,fNoper-1,err)) return;
            }
         } else if (lshift != 0) {
            if (lshift == 1 || lshift == lchain) {
               err=5;
               chaine_error=">>";
            } else {
               ctemp = chaine(0,lshift-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               UInt_t leftopr = fNoper-1;
               ctemp = chaine(lshift+1,lchain-lshift-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               fExpr[fNoper] = ">>";
               actionCode = kLeftShift;
               SetAction(fNoper,actionCode,actionParam);
               fNoper++;
               if (!CheckOperands(leftopr,fNoper-1,err)) return;
            }
         } else {
            if (multi != 0) {
               if (multi == 1 || multi == lchain) {
               err=5;
               chaine_error="*";
            }
            else {
               ctemp = chaine(0,multi-1);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               UInt_t leftopr = fNoper-1;
               ctemp = chaine(multi,lchain-multi);
               Analyze(ctemp.Data(),err,offset); if (err) return;
               fExpr[fNoper] = "*";
               actionCode = kMultiply;
               SetAction(fNoper,actionCode,actionParam);
               fNoper++;
               if (!CheckOperands(leftopr,fNoper-1,err)) return;
            }
         } else {
            if (divi != 0) {
               if (divi == 1 || divi == lchain) {
                  err=5;
                  chaine_error = "/";
               }
               else {
                  ctemp = chaine(0,divi-1);
                  Analyze(ctemp.Data(),err,offset); if (err) return;
                  UInt_t leftopr = fNoper-1;
                  ctemp = chaine(divi,lchain-divi);
                  Analyze(ctemp.Data(),err,offset); if (err) return;
                  fExpr[fNoper] = "/";
                  actionCode = kDivide;
                  SetAction(fNoper,actionCode,actionParam);
                  fNoper++;
                  if (!CheckOperands(leftopr,fNoper-1,err)) return;
               }
            } else {
               if (puiss != 0) {
                  if (puiss == 1 || puiss == lchain) {
                     err = 5;
                     chaine_error = "**";
                  }
                  else {
                     if (chaine(lchain-2,2) == "^2") {
                        ctemp = "sq(" + chaine(0,lchain-2) + ")";
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                     } else {
                        ctemp = chaine(0,puiss-1);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        UInt_t leftopr = fNoper-1;
                        ctemp = chaine(puiss,lchain-puiss);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "^";
                        actionCode = kpow;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(leftopr,fNoper-1,err)) return;
                     }
                  }
               } else {

                  find=0;

   // Check for a numerical expression
                  {
                     Bool_t hasDot = kFALSE;
                     Bool_t isHexa = kFALSE;
                     Bool_t hasExpo= kFALSE;
                     if ((chaine(0,2)=="0x")||(chaine(0,2)=="0X")) isHexa=kTRUE;
                     for (j=0; j<chaine.Length() && err==0; j++) {
                        t=chaine[j];
                        if (!isHexa) {
                           if (j>0 && (chaine(j,1)=="e" || chaine(j,2)=="e+" || chaine(j,2)=="e-" || chaine(j,1)=="E" || chaine(j,2)=="E+" || chaine(j,2)=="E-")) {
                              if (hasExpo) {
                                 err=26;
                                 chaine_error=chaine;
                              }
                              hasExpo = kTRUE;
                              // The previous implementation allowed a '.' in the exponent.
                              // That information was ignored (by sscanf), we now make it an error
                              // hasDot = kFALSE;
                              hasDot = kTRUE;  // forbid any additional '.'
                              if (chaine(j,2)=="e+" || chaine(j,2)=="e-" || chaine(j,2)=="E+" || chaine(j,2)=="E-") j++;
                           }
                           else {
                              if (chaine(j,1) == "." && !hasDot) hasDot = kTRUE; // accept only one '.' in the number
                              else {
                                 // The previous implementation was allowing ANYTHING after the '.' and thus
                                 // made legal code like '2.3 and fpx' and was just silently ignoring the
                                 // 'and fpx'.
                                 if (!strchr("0123456789",t) && (chaine(j,1)!="+" || j!=0)) {
                                    err = 30;
                                    chaine_error=chaine;
                                 }
                              }
                           }
                        }
                        else {
                           if (!strchr("0123456789abcdefABCDEF",t) && (j>1)) {
                              err = 30;
                              chaine_error=chaine;
                           }
                        }
                     }
                     if (fNconst >= gMAXCONST) err = 27;
                     if (!err) {
                        if (!isHexa) {if (sscanf((const char*)chaine,"%lg",&vafConst) > 0) err = 0; else err =1;}
                        else {if (sscanf((const char*)chaine,"%lx",&vafConst2) > 0) err = 0; else err=1;
                        vafConst = (Double_t) vafConst2;}
                        fExpr[fNoper] = chaine;
                        k = -1;
                        for (j=0;j<fNconst;j++) {
                           if (vafConst == fConst[j] ) k= j;
                        }
                        if ( k < 0) {  k = fNconst; fNconst++; fConst[k] = vafConst; }
                        actionCode = kConstant;
                        actionParam = k;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                     }
                     if (err==30) err=0;
                     else find = kTRUE;
                  }

                  // Look for an already defined expression

                  if (find==0) {
                     {
                        R__LOCKGUARD(gROOTMutex);
                        oldformula = (const TFormula*)gROOT->GetListOfFunctions()->FindObject((const char*)chaine);
                     }
                     if (oldformula && strcmp(schain,oldformula->GetTitle())) {
                        Int_t nprior = fNpar;
                        Analyze(oldformula->GetExpFormula(),err,fNpar);
                        // if the oldformula was using a normalized function (gausn or landaun) set also in this one
                        if (oldformula->IsNormalized()) SetBit(kNormalized);
                        if (err) return; // changes fNpar
                        fNpar = nprior;
                        find=1;
                        if (!err) {
                           Int_t npold = oldformula->GetNpar();
                           fNpar += npold;
                           for (Int_t ipar=0;ipar<npold;ipar++) {
                              fParams[ipar+fNpar-npold] = oldformula->GetParameter(ipar);
                           }
                        }
                     }
                  }
                  if (find == 0) {

   // Check if chaine is a defined variable.
   // Note that DefinedVariable can be overloaded

                     ctemp = chaine;
                     ctemp.ReplaceAll(escapedSlash, slash);
                     Int_t action;
                     k = DefinedVariable(ctemp,action);
                     if (k==-3) {
                        // Error message already issued
                        err = 1;
                     } else if (k==-2) {
                        err = 31;
                        chaine_error = ctemp;
                     } else if ( k >= 0 ) {
                        fExpr[fNoper] = ctemp;
                        actionCode = action;
                        actionParam = k;
                        SetAction(fNoper,actionCode,actionParam);
                        if (action==kDefinedString) fNstring++;
                        else if (k <kMAXFOUND && !fAlreadyFound.TestBitNumber(k)) {
                           fAlreadyFound.SetBitNumber(k);
                           fNval++;
                        }
                        fNoper++;
                     } else if (chaine(0,1) == "!") {
                        ctemp = chaine(1,lchain-1);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "!";
                        actionCode = kNot;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,1)=="\"" && chaine(chaine.Length()-1,1)=="\"") {
                        // It is a string !!!
                        fExpr[fNoper] = chaine(1,chaine.Length()-2);
                        actionCode = kStringConst;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                     } else if (chaine(0,4) == "cos(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "cos";
                        actionCode = kcos;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "sin(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "sin";
                        actionCode = ksin;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "tan(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "tan";
                        actionCode = ktan;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "acos(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "acos";
                        actionCode = kacos;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "asin(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "asin";
                        actionCode = kasin;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "atan(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "atan";
                        actionCode = katan;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "cosh(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "cosh";
                        actionCode = kcosh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "sinh(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "sinh";
                        actionCode = ksinh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "tanh(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "tanh";
                        actionCode = ktanh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,6) == "acosh(") {
                        ctemp = chaine(5,lchain-5);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "acosh";
                        actionCode = kacosh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,6) == "asinh(") {
                        ctemp = chaine(5,lchain-5);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "asinh";
                        actionCode = kasinh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,6) == "atanh(") {
                        ctemp = chaine(5,lchain-5);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "atanh";
                        actionCode = katanh;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,3) == "sq(") {
                        ctemp = chaine(2,lchain-2);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "sq";
                        actionCode = ksq;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "log(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "log";
                        actionCode = klog;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,6) == "log10(") {
                        ctemp = chaine(5,lchain-5);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "log10";
                        actionCode = klog10;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "exp(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "exp";
                        actionCode = kexp;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "abs(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "abs";
                        actionCode = kabs;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,5) == "sign(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "sign";
                        actionCode = ksign;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine(0,4) == "int(") {
                        ctemp = chaine(3,lchain-3);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "int";
                        actionCode = kint;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;
                     } else if (chaine == "rndm" || chaine(0,5) == "rndm(") {
                        fExpr[fNoper] = "rndm";
                        actionCode = krndm;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                     } else if (chaine(0,5) == "sqrt(") {
                        ctemp = chaine(4,lchain-4);
                        Analyze(ctemp.Data(),err,offset); if (err) return;
                        fExpr[fNoper] = "sqrt";
                        actionCode = ksqrt;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;;
                        if (!CheckOperands(fNoper-1,err)) return;

   // Look for an exponential

                     } else if ( chaine == "expo" || chaine(0,5)=="expo("
                            || (lchain==5 && chaine(1,4)=="expo")
                            || (lchain==6 && chaine(2,4)=="expo")
                            || chaine(1,5)=="expo(" || chaine(2,5)=="expo(" ) {
                        chaine1ST=chaine;
                        if (chaine(1,4) == "expo") {
                           ctemp=chaine(0,1);
                           if (ctemp=="x") {
                              inter2=0;
                              if (fNdim < 1) fNdim = 1; }
                           else if (ctemp=="y") {
                              inter2=1;
                              if (fNdim < 2) fNdim = 2; }
                           else if (ctemp=="z") {
                              inter2=2;
                              if (fNdim < 3) fNdim = 3; }
                           else if (ctemp=="t") {
                              inter2=3;
                              if (fNdim < 4) fNdim = 4; }
                           else {
                              err=26; // unknown name;
                              chaine_error=chaine1ST;
                           }
                           chaine=chaine(1,lchain-1);
                           lchain=chaine.Length();
                        } else inter2=0;
                        if (chaine(2,4) == "expo") {
                           if (chaine(0,2) != "xy") {
                              err=26; // unknown name
                              chaine_error=chaine1ST;
                           }
                           else {
                              inter2=5;
                              if (fNdim < 2) fNdim = 2;
                              chaine=chaine(2,lchain-2);
                              lchain=chaine.Length();
                           }
                        }
                        if (lchain == 4) {
                           if (fNpar>=gMAXPAR) err=7; // too many parameters
                           if (!err) {
                              fExpr[fNoper] = chaine1ST;
                              actionCode = kexpo + inter2;
                              actionParam = offset;
                              SetAction(fNoper,actionCode,actionParam);
                              if (inter2 == 5+offset && fNpar < 3+offset) fNpar = 3+offset;
                              if (fNpar < 2+offset) fNpar = 2+offset;
                              if (fNpar>=gMAXPAR) err=7; // too many parameters
                              if (!err) {
                                 fNoper++;
                                 if (fNdim < 1) fNdim = 1;
                                 if (fNpar == 2) SetNumber(200);
                              }
                           }
                        } else if (chaine(4,1) == "(") {
                           ctemp = chaine(5,lchain-6);
                           fExpr[fNoper] = chaine1ST;
                           for (j=0; j<ctemp.Length(); j++) {
                              t=ctemp[j];
                              if (strchr("0123456789",t)==0 && (ctemp(j,1)!="+" || j!=0)) {
                                 err=20;
                                 chaine_error=chaine1ST;
                              }
                           }
                           if (err==0) {
                              sscanf(ctemp.Data(),"%d",&inter);
                              if (inter>=0) {
                                 inter += offset;
                                 actionCode = kexpo + inter2;
                                 actionParam = inter;
                                 SetAction(fNoper,actionCode,actionParam);
                                 if (inter2 == 5) inter++;
                                 if (inter+2>fNpar) fNpar = inter+2;
                                 if (fNpar>=gMAXPAR) err=7; // too many parameters
                                 if (!err) fNoper++;
                                 if (fNpar == 2) SetNumber(200);
                              } else err=20;
                           } else err = 20; // non integer value for parameter number
                        } else {
                           err=26; // unknown name
                           chaine_error=chaine;
                        }

   // Look for gaus, xgaus,ygaus,xygaus

                     } else if (chaine=="gaus"
                            || (lchain==5 && chaine(1,4)=="gaus")
                            || (lchain==6 && chaine(2,4)=="gaus")
                            || chaine(0,5)=="gaus(" || chaine(1,5)=="gaus(" || chaine(2,5)=="gaus(") {
                        chaine1ST=chaine;
                        if (chaine(1,4) == "gaus") {
                           ctemp=chaine(0,1);
                           if (ctemp=="x") {
                              inter2=0;
                              if (fNdim < 1) fNdim = 1; }
                           else if (ctemp=="y") {
                              inter2=1;
                              if (fNdim < 2) fNdim = 2; }
                           else if (ctemp=="z") {
                              inter2=2;
                              if (fNdim < 3) fNdim = 3; }
                           else if (ctemp=="t") {
                              inter2=3;
                              if (fNdim < 4) fNdim = 4; }
                           else {
                              err=26; // unknown name
                              chaine_error=chaine1ST;
                           }
                           chaine=chaine(1,lchain-1);
                           lchain=chaine.Length();
                        } else inter2=0;
                        if (chaine(2,4) == "gaus") {
                           if (chaine(0,2) != "xy") {
                              err=26; // unknown name
                              chaine_error=chaine1ST;
                           }
                           else {
                              inter2=5;
                              if (fNdim < 2) fNdim = 2;
                              chaine=chaine(2,lchain-2);
                              lchain=chaine.Length();
                              SetNumber(110); // xygaus
                           }
                        }
                        if (lchain == 4 && err==0) {
                           if (fNpar>=gMAXPAR) err=7; // too many parameters
                           if (!err) {
                              fExpr[fNoper] = chaine1ST;
                              actionCode = kgaus + inter2;
                              actionParam = offset;
                              SetAction(fNoper,actionCode,actionParam);
                              if (inter2 == 5+offset && fNpar < 5+offset) fNpar = 5+offset;
                              if (3+offset>fNpar) fNpar = 3+offset;
                              if (fNpar>=gMAXPAR) err=7; // too many parameters
                              if (!err) {
                                 fNoper++;
                                 if (fNdim < 1) fNdim = 1;
                                 if (fNpar == 3) SetNumber(100);
                              }
                           }
                        } else if (chaine(4,1) == "(" && err==0) {
                           ctemp = chaine(5,lchain-6);
                           fExpr[fNoper] = chaine1ST;
                           for (j=0; j<ctemp.Length(); j++) {
                              t=ctemp[j];
                              if (strchr("0123456789",t)==0 && (ctemp(j,1)!="+" || j!=0)) {
                                 err=20;
                                 chaine_error=chaine1ST;
                              }
                           }
                           if (err==0) {
                              sscanf(ctemp.Data(),"%d",&inter);
                              if (inter >= 0) {
                                 inter += offset;
                                 actionCode = kgaus + inter2;
                                 actionParam = inter;
                                 SetAction(fNoper,actionCode,actionParam);
                                 if (inter2 == 5) inter += 2;
                                 if (inter+3>fNpar) fNpar = inter+3;
                                 if (fNpar>=gMAXPAR) err=7; // too many parameters
                                 if (!err) fNoper++;
                                 if(fNpar == 3) SetNumber(100);
                              } else err = 20; // non integer value for parameter number
                           }
                        } else if (err==0) {
                           err=26; // unknown name
                           chaine_error=chaine1ST;
                        }

   // Look for landau, xlandau,ylandau,xylandau

                     } else if (chaine=="landau" || (lchain==7 && chaine(1,6)=="landau")
                            || (lchain==8 && chaine(2,6)=="landau")
                            || chaine(0,7)=="landau(" || chaine(1,7)=="landau(" || chaine(2,7)=="landau(") {
                        chaine1ST=chaine;
                        if (chaine(1,6) == "landau") {
                           ctemp=chaine(0,1);
                           if (ctemp=="x") {
                              inter2=0;
                              if (fNdim < 1) fNdim = 1; }
                           else if (ctemp=="y") {
                              inter2=1;
                              if (fNdim < 2) fNdim = 2; }
                           else if (ctemp=="z") {
                              inter2=2;
                              if (fNdim < 3) fNdim = 3; }
                           else if (ctemp=="t") {
                              inter2=3;
                              if (fNdim < 4) fNdim = 4; }
                           else {
                              err=26; // unknown name
                              chaine_error=chaine1ST;
                           }
                           chaine=chaine(1,lchain-1);
                           lchain=chaine.Length();
                        } else inter2=0;
                        if (chaine(2,6) == "landau") {
                           if (chaine(0,2) != "xy") {
                              err=26; // unknown name
                              chaine_error=chaine1ST;
                           }
                           else {
                              inter2=5;
                              if (fNdim < 2) fNdim = 2;
                              chaine=chaine(2,lchain-2);
                              lchain=chaine.Length();
                              SetNumber(410);
                           }
                        }
                        if (lchain == 6 && err==0) {
                           if (fNpar>=gMAXPAR) err=7; // too many parameters
                           if (!err) {
                              fExpr[fNoper] = chaine1ST;
                              actionCode = klandau + inter2;
                              actionParam = offset;
                              SetAction(fNoper,actionCode,actionParam);
                              if (inter2 == 5+offset && fNpar < 5+offset) fNpar = 5+offset;
                              if (3+offset>fNpar) fNpar = 3+offset;
                              if (fNpar>=gMAXPAR) err=7; // too many parameters
                              if (!err) {
                                 fNoper++;
                                 if (fNdim < 1) fNdim = 1;
                                 if (fNpar == 3) SetNumber(400);
                              }
                           }
                        } else if (chaine(6,1) == "(" && err==0) {
                           ctemp = chaine(7,lchain-8);
                           fExpr[fNoper] = chaine1ST;
                           for (j=0; j<ctemp.Length(); j++) {
                              t=ctemp[j];
                              if (strchr("0123456789",t)==0 && (ctemp(j,1)!="+" || j!=0)) {
                                 err=20;
                                 chaine_error=chaine1ST;
                              }
                           }
                           if (err==0) {
                              sscanf(ctemp.Data(),"%d",&inter);
                              if (inter >= 0) {
                                 inter += offset;
                                 actionCode = klandau + inter2;
                                 actionParam = inter;
                                 SetAction(fNoper,actionCode,actionParam);
                                 if (inter2 == 5) inter += 2;
                                 if (inter+3>fNpar) fNpar = inter+3;
                                 if (fNpar>=gMAXPAR) err=7; // too many parameters
                                 if (!err) fNoper++;
                                 if (fNpar == 3) SetNumber(400);
                              } else err = 20; // non integer value for parameter number
                           }
                        } else if (err==0) {
                           err=26; // unknown name
                           chaine_error=chaine1ST;
                        }

   // Look for a polynomial

                     } else if (chaine(0,3) == "pol" || chaine(1,3) == "pol") {
                        chaine1ST=chaine;
                        if (chaine(1,3) == "pol") {
                           ctemp=chaine(0,1);
                           if (ctemp=="x") {
                              inter2=1;
                              if (fNdim < 1) fNdim = 1; }
                           else if (ctemp=="y") {
                              inter2=2;
                              if (fNdim < 2) fNdim = 2; }
                           else if (ctemp=="z") {
                              inter2=3;
                              if (fNdim < 3) fNdim = 3; }
                           else if (ctemp=="t") {
                              inter2=4;
                              if (fNdim < 4) fNdim = 4; }
                           else {
                              err=26; // unknown name;
                              chaine_error=chaine1ST;
                           }
                           chaine=chaine(1,lchain-1);
                           lchain=chaine.Length();
                        } else inter2=1;
                        if (chaine(lchain-1,1) == ")") {
                           nomb = 0;
                           for (j=3;j<lchain;j++) if (chaine(j,1)=="(" && nomb == 0) nomb = j;
                           if (nomb == 3) err = 23; // degree of polynomial not specified
                           if (nomb == 0) err = 40; // '(' is expected
                           ctemp = chaine(nomb+1,lchain-nomb-2);
                           for (j=0; j<ctemp.Length(); j++) {
                              t=ctemp[j];
                              if (strchr("0123456789",t)==0 && (ctemp(j,1)!="+" || j!=0)) {
                                 err=20;
                                 chaine_error=chaine1ST;
                              }
                           }
                           if (!err) {
                              sscanf(ctemp.Data(),"%d",&inter);
                              if (inter < 0) err = 20;
                           }
                        }
                        else {
                           nomb = lchain;
                           inter = 0;
                        }
                        if (!err) {
                           inter--;
                           ctemp = chaine(3,nomb-3);
                           if (sscanf(ctemp.Data(),"%d",&n) > 0) {
                              if (n < 0  ) err = 24; //Degree of polynomial must be positive
                              if (n >= 20) err = 25; //Degree of polynomial must be less than 20
                           } else err = 20;
                        }
                        if (!err) {
                           fExpr[fNoper] = chaine1ST;
                           actionCode = kpol+(inter2-1);
                           actionParam = n*100+inter+2;
                           SetAction(fNoper,actionCode,actionParam);
                           if (inter+n+1>=fNpar) fNpar = inter + n + 2;
                           if (fNpar>=gMAXPAR) err=7; // too many parameters
                           if (!err) {
                              fNoper++;
                              if (fNdim < 1) fNdim = 1;
                              SetNumber(300+n);
                           }
                        }

   // Look for pow,atan2,etc

                     } else if (chaine(0,4) == "pow(") {
                        compt = 4; nomb = 0; virgule = 0; nest=0;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "(") nest++;
                           else if (chaine(compt-1,1) == ")") nest--;
                           else if (chaine(compt-1,1) == "," && nest==0) {
                              nomb++;
                              if (nomb == 1 && virgule == 0) virgule = compt;
                           }
                        }
                        if (nomb != 1) err = 22; // There are plus or minus than 2 arguments for pow
                        else {
                           ctemp = chaine(4,virgule-5);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           UInt_t leftopr = fNoper-1;
                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "^";
                           actionCode = kpow;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;
                           if (!CheckOperands(leftopr,fNoper-1,err)) return;
                        }
                     } else if (chaine(0,7) == "strstr(") {
                        compt = 7; nomb = 0; virgule = 0; nest=0;
                        inString = false;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "\"") {
                              inString = !inString;
                           }  else if (!inString) {
                              if (chaine(compt-1,1) == "(") nest++;
                              else if (chaine(compt-1,1) == ")") nest--;
                              else if (chaine(compt-1,1) == "," && nest==0) {
                                 nomb++;
                                 if (nomb == 1 && virgule == 0) virgule = compt;
                              }
                           }
                        }
                        if (nomb != 1) err = 28; // There are plus or minus than 2 arguments for strstr
                        else {
                           ctemp = chaine(7,virgule-8);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           Int_t optloc = fNoper-1;

                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "strstr";
                           actionCode = kstrstr;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;

                           if ( !IsString(optloc) || !IsString(fNoper-2) ) {
                              err = 46;
                              chaine_error = "strstr";
                           }
                        }
                     } else if (chaine(0,4) == "min(") {
                        compt = 4; nomb = 0; virgule = 0; nest=0;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "(") nest++;
                           else if (chaine(compt-1,1) == ")") nest--;
                           else if (chaine(compt-1,1) == "," && nest==0) {
                              nomb++;
                              if (nomb == 1 && virgule == 0) virgule = compt;
                           }
                        }
                        if (nomb != 1) {
                           err = 44; // There are plus or minus than 2 arguments for min
                           err_hint = 3;
                        }
                        else {
                           ctemp = chaine(4,virgule-5);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           UInt_t leftopr = fNoper-1;
                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "min";
                           actionCode = kmin;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;
                           if (!CheckOperands(leftopr,fNoper-1,err)) return;
                        }
                     } else if (chaine(0,4) == "max(") {
                        compt = 4; nomb = 0; virgule = 0; nest=0;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "(") nest++;
                           else if (chaine(compt-1,1) == ")") nest--;
                           else if (chaine(compt-1,1) == "," && nest==0) {
                              nomb++;
                              if (nomb == 1 && virgule == 0) virgule = compt;
                           }
                        }
                        if (nomb != 1) {
                           err = 44; // There are plus or minus than 2 arguments for min
                           err_hint = 3;
                        }
                        else {
                           ctemp = chaine(4,virgule-5);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           UInt_t leftopr = fNoper-1;
                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "max";
                           actionCode = kmax;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;
                           if (!CheckOperands(leftopr,fNoper-1,err)) return;
                        }

                     } else if (chaine(0,6) == "atan2(") {
                        compt = 6; nomb = 0; virgule = 0; nest=0;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "(") nest++;
                           else if (chaine(compt-1,1) == ")") nest--;
                           else if (chaine(compt-1,1) == "," && nest==0) {
                              nomb++;
                              if (nomb == 1 && virgule == 0) virgule = compt;
                           }
                        }
                        if (nomb != 1) err = 21;  //{ There are plus or minus than 2 arguments for atan2
                        else {
                           ctemp = chaine(6,virgule-7);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           UInt_t leftopr = fNoper-1;
                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "atan2";
                           actionCode = katan2;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;
                           if (!CheckOperands(leftopr,fNoper-1,err)) return;
                        }
                     } else if (chaine(0,5) == "fmod(") {
                        compt = 5; nomb = 0; virgule = 0; nest=0;
                        while(compt != lchain) {
                           compt++;
                           if (chaine(compt-1,1) == "(") nest++;
                           else if (chaine(compt-1,1) == ")") nest--;
                           else if (chaine(compt-1,1) == "," && nest==0) {
                              nomb++;
                              if (nomb == 1 && virgule == 0) virgule = compt;
                           }
                        }
                        if (nomb != 1) {
                           err = 44; // There are plus or minus than 2 arguments for fmod
                           err_hint = 4;
                        }
                        else {
                           ctemp = chaine(5,virgule-6);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           UInt_t leftopr = fNoper-1;
                           ctemp = chaine(virgule,lchain-virgule-1);
                           Analyze(ctemp.Data(),err,offset); if (err) return;
                           fExpr[fNoper] = "fmod";
                           actionCode = kfmod;
                           SetAction(fNoper,actionCode,actionParam);
                           fNoper++;
                           if (!CheckOperands(leftopr,fNoper-1,err)) return;
                        }
                     } else if (AnalyzeFunction(chaine,err,offset) || err) { // The '||err' is to grab an error coming from AnalyzeFunction
                        if (err) {
                           chaine_error = chaine;
                        } else {
                           // We have a function call. Note that all the work was already,
                           // eventually done in AnalyzeFunction
                           //fprintf(stderr,"We found a foreign function in %s\n",chaine.Data());
                        }
                     } else if (chaine(0,1) == "[" && chaine(lchain-1,1) == "]") {
                        fExpr[fNoper] = chaine;
                        fNoper++;
                        ctemp = chaine(1,lchain-2);
                        for (j=0; j<ctemp.Length(); j++) {
                           t=ctemp[j];
                           if (strchr("0123456789",t)==0 && (ctemp(j,1)!="+" || j!=0)) {
                              err=20;
                              chaine_error=chaine1ST; // le numero ? de par[?] n'est pas un entier }
                           }
                        }
                        if (!err) {
                           sscanf(ctemp.Data(),"%d",&valeur);
                           actionCode = kParameter;
                           actionParam = offset + valeur;
                           SetAction(fNoper-1, actionCode, actionParam);
                           fExpr[fNoper-1] = "[";
                           fExpr[fNoper-1] = (fExpr[fNoper-1] + (long int)(valeur+offset)) + "]";
                        }
                     } else if (chaine == "pi") {
                        fExpr[fNoper] = "pi";
                        actionCode = kpi;
                        SetAction(fNoper,actionCode,actionParam);
                        fNoper++;
                     }
                     else {

   // None of the above.

                        err = 30;
                     }
                  }
               }
            }
         }
      }
   }

   // Overflows
   if (fNoper>=gMAXOP) err=6; // too many operators

   }

   // errors!
   if (err>1) {
      TString er = "";
      chaine_error = "\""+chaine_error+"\"";
      switch(err) {
         case  2 : er = " Invalid Floating Point Operation"; break;
         case  4 : er = " Empty String"; break;
         case  5 : er = " Invalid Syntax " + chaine_error; break;
         case  6 : er = " Too many operators !"; break;
         case  7 : er = " Too many parameters !"; break;
         case 10 : er = " z specified but not x and y"; break;
         case 11 : er = " z and y specified but not x"; break;
         case 12 : er = " y specified but not x"; break;
         case 13 : er = " z and x specified but not y"; break;
         case 20 : er = " Non integer value for parameter number : " + chaine_error; break;
         case 21 : er = " ATAN2 requires two arguments"; break;
         case 22 : er = " POW requires two arguments"; break;
         case 23 : er = " Degree of polynomial not specified"; break;
         case 24 : er = " Degree of polynomial must be positive"; break;
         case 25 : er = " Degree of polynomial must be less than 20"; break;
         case 26 : er = " Unknown name : " + chaine_error; break;
         case 27 : er = " Too many constants in expression"; break;
         case 28 : er = " strstr requires two arguments"; break;
         case 29 : er = " TFormula can only call interpreted and compiled functions that return a numerical type: " + chaine_error; break;
         case 30 : er = " Bad numerical expression : " + chaine_error; break;
         case 31 : er = " Part of the Variable " + chaine_error; er += " exists but some of it is not accessible or useable"; break;
         case 40 : er = " '(' is expected"; break;
         case 41 : er = " ')' is expected"; break;
         case 42 : er = " '[' is expected"; break;
         case 43 : er = " ']' is expected"; break;
         case 44 : er = " The function '" + chaine(0,err_hint) + "' requires two arguments."; break;
         case 45 : er = "The operator " + chaine_error + " requires a numerical operand."; break;
         case 46 : er = "Both operands of the operator " + chaine_error + " have to be either numbers or strings."; break;
         case 47 : er = chaine_error + " requires 2 string arguments"; break;
      }
      Error("Compile", "%s", er.Data());
      err=1;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Check whether the operand at 'oper-1' is compatible with the operation
/// at 'oper'.

Bool_t TFormula::CheckOperands(Int_t oper, Int_t &err)
{
   if ( IsString(oper-1) && !StringToNumber(oper-1) ) {
      Error("Compile","\"%s\" requires a numerical operand.",fExpr[oper].Data());
      err = 45;
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether the operands at 'leftoper' and 'oper-1' are compatible with
/// the operation at 'oper'.

Bool_t TFormula::CheckOperands(Int_t leftoper, Int_t oper, Int_t &err)
{
   if ( IsString(oper-1) || IsString(leftoper) ) {
      if (IsString(oper-1) && StringToNumber(oper-1)) {
          return kTRUE;
      }
      if (IsString(leftoper) && StringToNumber(leftoper)) {
         return kTRUE;
      }
      Error("Compile","\"%s\" requires two numerical operands.",fExpr[oper].Data());
      err = 46;
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Try to 'demote' a string into an array bytes.  If this is not possible,
/// return false.

Bool_t TFormula::StringToNumber(Int_t /* code */)
{
   // In TFormula proper, we can not handle array of bytes ...
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the objects.
///
/// Resets the object to its state before compilation.

void TFormula::Clear(Option_t * /*option*/ )
{
   ClearFormula();
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the objects.
///
/// Resets the object to its state before compilation.

void TFormula::ClearFormula(Option_t * /*option*/ )
{
   fNdim   = 0;
   fNpar   = 0;
   fNoper  = 0;
   fNconst = 0;
   fNumber = 0;
   fNstring= 0;
   fNval   = 0;

   if (fExpr)   { delete [] fExpr;   fExpr   = 0;}
   if (fNames)  { delete [] fNames;  fNames  = 0;}
   if (fOper)   { delete [] fOper;   fOper   = 0;}
   if (fConst)  { delete [] fConst;  fConst  = 0;}
   if (fParams) { delete [] fParams; fParams = 0;}
   fFunctions.Delete();
   fLinearParts.Delete();
   //
   //MI change
   if (fPredefined)    { delete [] fPredefined;    fPredefined    = 0;}
   if (fOperOffset)    { delete [] fOperOffset;    fOperOffset    = 0;}
   if (fExprOptimized) { delete [] fExprOptimized; fExprOptimized = 0;}
   if (fOperOptimized) { delete [] fOperOptimized; fOperOptimized = 0;}
   // should we also remove the object from the list?
   // gROOT->GetListOfFunctions()->Remove(this);
   // if we don't, what happens if it fails the new compilation?
}

namespace {
   template <class T>
   inline static void ResizeArrayIfAllocated(T*& oldArray, int newSize){

      // Don't do anything in this case.
      if (!oldArray || newSize <=0) return;

      T* newArray = new T[newSize];
      std::copy(oldArray, oldArray+newSize, newArray);
      delete [] oldArray;
      oldArray = newArray;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compile expression already stored in fTitle.
///
///   Loop on all subexpressions of formula stored in fTitle
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///     TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///     MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.
///
/// \image html TFormula_compile.png

Int_t TFormula::Compile(const char *expression)
{
   Int_t i,j,lc,valeur,err;
   TString ctemp;

   ClearFormula();

   // If expression is not empty, take it, otherwise take the title
   if (strlen(expression)) SetTitle(expression);

   TString chaine = GetTitle();

   if (chaine.Contains(";")) {
      char *sctemp = new char[chaine.Length()+1];
      strlcpy(sctemp,chaine.Data(),chaine.Length()+1);
      char *semicol = (char*)strstr(sctemp,";");
      if (semicol) *semicol = 0;
      chaine = sctemp;
      delete [] sctemp;
   }

   // if the function is linear, process it and fill the array of linear parts
   if (TestBit(kLinear)){
      ProcessLinear(chaine);
   }

     // see static function SetMaxima to change the gMAX.. default values
   fExpr   = new TString[gMAXOP];
   fConst  = new Double_t[gMAXCONST];
   fParams = new Double_t[gMAXPAR];
   fNames  = new TString[gMAXPAR];
   fOper   = new Int_t[gMAXOP];
   for (i=0; i<gMAXPAR; i++) {
      fParams[i] = 0;
      fNames[i] = "";
   }
   for (i=0; i<gMAXOP; i++) {
      fExpr[i] = "";
      fOper[i] = 0;
   }
   for (i=0; i<gMAXCONST; i++)
      fConst[i] = 0;

   // Substitution of some operators to C++ style
   Bool_t inString = false;
   for (i=1; i<=chaine.Length(); i++) {
      lc =chaine.Length();
      if (chaine(i-1,1) == "\"") inString = !inString;
      if (inString) continue;
      if (chaine(i-1,2) == "**") {
         chaine = chaine(0,i-1) + "^" + chaine(i+1,lc-i-1);
         i=0;
      } else if (chaine(i-1,2) == "++") {
         chaine = chaine(0,i) + chaine(i+1,lc-i-1);
         i=0;
      } else if (chaine(i-1,2) == "+-" || chaine(i-1,2) == "-+") {
         chaine = chaine(0,i-1) + "-" + chaine(i+1,lc-i-1);
         i=0;
      } else if (chaine(i-1,2) == "--") {
         chaine = chaine(0,i-1) + "+" + chaine(i+1,lc-i-1);
         i=0;
      } else if (chaine(i-1,2) == "->") {
         chaine = chaine(0,i-1) + "." + chaine(i+1,lc-i-1);
         i=0;
      } else if (chaine(i-1,1) == "[") {
         for (j=1;j<=chaine.Length()-i;j++) {
            if (chaine(j+i-1,1) == "]" || j+i > chaine.Length()) break;
         }
         ctemp = chaine(i,j-1);
         valeur=0;
         sscanf(ctemp.Data(),"%d",&valeur);
         if (valeur >= fNpar) fNpar = valeur+1;
      } else if (chaine(i-1,1) == " ") {
         chaine = chaine(0,i-1)+chaine(i,lc-i);
         i=0;
      }
   }
   err = 0;
   Analyze((const char*)chaine,err);

   // if no parameters delete arrays fParams and fNames
   if (!fNpar) {
      delete [] fParams; fParams = 0;
      delete [] fNames;  fNames = 0;
   }

   // if no errors, copy local parameters to formula objects
   if (!err) {
      if (fNdim <= 0) fNdim = 1;
      if (chaine.Length() > 4)
      {
         if ( GetNumber() != 400 &&
              GetNumber() != 410 &&
              GetNumber() != 110 )
            SetNumber(0);
         else if ( GetNumber() == 110 && chaine.Length() > 6 )
            SetNumber(0);
         else if ( GetNumber() == 410 && chaine.Length() > 8 )
            SetNumber(0);
      }
      // if formula is a gaussian, set parameter names
      if (GetNumber() == 100) {
         SetParName(0,"Constant");
         SetParName(1,"Mean");
         SetParName(2,"Sigma");
      }
      // if formula is a 2D gaussian, set parameter names
      if (GetNumber() == 110){
         SetParName(0,"Constant");
         SetParName(1,"MeanX");
         SetParName(2,"SigmaX");
         SetParName(3,"MeanY");
         SetParName(4,"SigmaY");
      }
      // if formula is an exponential, set parameter names
      if (GetNumber() == 200) {
         SetParName(0,"Constant");
         SetParName(1,"Slope");
      }
      // if formula is a polynom, set parameter names
      if (GetNumber() == 300+fNpar) {
         for (i=0;i<fNpar;i++) SetParName(i,Form("p%d",i));
      }
      // if formula is a landau, set parameter names
      if (GetNumber() == 400) {
         SetParName(0,"Constant");
         SetParName(1,"MPV");
         SetParName(2,"Sigma");
      }
      // if formula is a 2D landau, set parameter names
      if (GetNumber() == 410) {
         SetParName(0,"Constant");
         SetParName(1,"MPVX");
         SetParName(2,"SigmaX");
         SetParName(3,"MPVY");
         SetParName(4,"SigmaY");
      }
   }

   // Here we shrink the arrays allocated like this:
   //    fExpr   = new TString[gMAXOP];
   //    fConst  = new Double_t[gMAXCONST];
   //    fParams = new Double_t[gMAXPAR];
   //    fNames  = new TString[gMAXPAR];
   //    fOper   = new Int_t[gMAXOP];
   // fParams and fNames may be already 0, so we have to check.
   if (!err){
      ResizeArrayIfAllocated(fExpr, fNoper);
      ResizeArrayIfAllocated(fConst, fNconst);
      ResizeArrayIfAllocated(fParams, fNpar);
      ResizeArrayIfAllocated(fNames, fNpar);
      ResizeArrayIfAllocated(fOper, fNoper);
      }


   if (err) { fNdim = 0; return 1; }
   //   Convert(5);
   //
   //MI change
   if (!IsA()->GetBaseClass("TTreeFormula")) {
      Optimize();
   }
   //
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this formula.

void TFormula::Copy(TObject &obj) const
{
   Int_t i;
   ((TFormula&)obj).ClearFormula();
   TNamed::Copy(obj);
   ((TFormula&)obj).fNdim   = fNdim;
   ((TFormula&)obj).fNpar   = fNpar;
   ((TFormula&)obj).fNoper  = fNoper;
   ((TFormula&)obj).fNconst = fNconst;
   ((TFormula&)obj).fNumber = fNumber;
   ((TFormula&)obj).fNval   = fNval;
   ((TFormula&)obj).fExpr   = 0;
   ((TFormula&)obj).fConst  = 0;
   ((TFormula&)obj).fParams = 0;
   ((TFormula&)obj).fNames  = 0;
   if (fExpr && fNoper) {
      ((TFormula&)obj).fExpr = new TString[fNoper];
      for (i=0;i<fNoper;i++)  ((TFormula&)obj).fExpr[i]   = fExpr[i];
   }
   if (fOper && fNoper) {
     ((TFormula&)obj).fOper = new Int_t[fNoper];
      for (i=0;i<fNoper;i++)  ((TFormula&)obj).fOper[i]   = fOper[i];
   }
   if (fConst && fNconst) {
      ((TFormula&)obj).fConst = new Double_t[fNconst];
      for (i=0;i<fNconst;i++) ((TFormula&)obj).fConst[i]  = fConst[i];
   }
   if (fParams && fNpar) {
      ((TFormula&)obj).fParams = new Double_t[fNpar];
      for (i=0;i<fNpar;i++)   ((TFormula&)obj).fParams[i] = fParams[i];
   }
   if (fNames && fNpar) {
      ((TFormula&)obj).fNames = new TString[fNpar];
      for (i=0;i<fNpar;i++)   ((TFormula&)obj).fNames[i]  = fNames[i];
   }

   TIter next(&fFunctions);
   TObject *fobj;
   while ( (fobj = next()) ) {
      ((TFormula&)obj).fFunctions.Add( fobj->Clone() );
   }
   //
   // MI change
   //
   //
   if (fNoper) {
      if(fExprOptimized) {
         ((TFormula&)obj).fExprOptimized   = new TString[fNoper];
         for (i=0;i<fNoper;i++)  ((TFormula&)obj).fExprOptimized[i]   = fExprOptimized[i];
      }
      if (fOperOptimized) {
         ((TFormula&)obj).fOperOptimized   = new Int_t[fNoper];
         for (i=0;i<fNoper;i++)  ((TFormula&)obj).fOperOptimized[i]   = fOperOptimized[i];
      }
      if (fPredefined) {
         ((TFormula&)obj).fPredefined      = new ROOT::v5::TFormulaPrimitive*[fNoper];
         for (i=0;i<fNoper;i++) {((TFormula&)obj).fPredefined[i] = fPredefined[i];}
      }
      if (fOperOffset) {
         ((TFormula&)obj).fOperOffset         = new TOperOffset[fNoper];
         for (i=0;i<fNoper;i++) {((TFormula&)obj).fOperOffset[i] = fOperOffset[i];}
      }
   }
   ((TFormula&)obj).fNOperOptimized = fNOperOptimized;
   ((TFormula&)obj).fOptimal = fOptimal;

}

////////////////////////////////////////////////////////////////////////////////
/// Return address of string corresponding to special code.
///
///   This member function is inactive in the TFormula class.
///   It may be redefined in derived classes.
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///     TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///     MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.

char *TFormula::DefinedString(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return value corresponding to special code.
///
///   This member function is inactive in the TFormula class.
///   It may be redefined in derived classes.
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///     TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///     MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.

Double_t TFormula::DefinedValue(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if expression is in the list of defined variables.
///
///   This member function can be overloaded in derived classes
///
///   If you overload this member function, you also HAVE TO
///   never call the constructor:
///
/// ~~~ {.cpp}
///     TFormula::TFormula(const char *name,const char *expression)
/// ~~~
///
///   and write your own constructor
///
/// ~~~ {.cpp}
///     MyClass::MyClass(const char *name,const char *expression) : TFormula()
/// ~~~
///
///   which has to call the TFormula default constructor and whose implementation
///   should be similar to the implementation of the normal TFormula constructor
///
///   This is necessary because the normal TFormula constructor call indirectly
///   the virtual member functions Analyze, DefaultString, DefaultValue
///   and DefaultVariable.
///
///   The expected returns values are
///   -   -2 :  the name has been recognized but won't be usable
///   -   -1 :  the name has not been recognized
///   -  >=0 :  the name has been recognized, return the action parameter.

Int_t TFormula::DefinedVariable(TString &chaine,Int_t &action)
{
   action = kVariable;
   if (chaine == "x") {
      if (fNdim < 1) fNdim = 1;
      return 0;
   } else if (chaine == "y") {
      if (fNdim < 2) fNdim = 2;
      return 1;
   } else if (chaine == "z") {
      if (fNdim < 3) fNdim = 3;
      return 2;
   } else if (chaine == "t") {
      if (fNdim < 4) fNdim = 4;
      return 3;
   }
   // MI change
   // extended defined variable (MI)
   //
   if (chaine.Data()[0]=='x'){
      if (chaine.Data()[1]=='[' && chaine.Data()[3]==']'){
         const char ch0 = '0';
         Int_t dim = chaine.Data()[2]-ch0;
         if (dim<0) return -1;
         if (dim>9) return -1;
         if (fNdim<=dim) fNdim = dim+1;
         return dim;
      }
      if (chaine.Data()[1]=='[' && chaine.Data()[4]==']'){
         const char ch0 = '0';
         Int_t dim = (chaine.Data()[2]-ch0)*10+(chaine.Data()[3]-ch0);
         if (dim<0) return -1;
         if (dim>99) return -1;
         if (fNdim<=dim) fNdim = dim+1;
         return dim;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate this formula.
///
///   The current value of variables x,y,z,t is passed through x, y, z and t.
///   The parameters used will be the ones in the array params if params is given
///    otherwise parameters will be taken from the stored data members fParams

Double_t TFormula::Eval(Double_t x, Double_t y, Double_t z, Double_t t) const
{
   Double_t xx[4];
   xx[0] = x;
   xx[1] = y;
   xx[2] = z;
   xx[3] = t;
   return ((TFormula*)this)->EvalPar(xx);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate this formula.
///
///   The current value of variables x,y,z,t is passed through the pointer x.
///   The parameters used will be the ones in the array params if params is given
///    otherwise parameters will be taken from the stored data members fParams
///
/// \image html TFormula_eval.png

Double_t TFormula::EvalParOld(const Double_t *x, const Double_t *uparams)
{
   Int_t i,j;
   // coverity[uninit] the tab value of tab is guaranteed to be set properly by the control flow.
   Double_t tab[kMAXFOUND];
   const char *stringStack[gMAXSTRINGFOUND] = {0};
   Double_t param_calc[kMAXFOUND];
   char *string_calc[gMAXSTRINGFOUND] = {0};
   Int_t precalculated = 0;
   Int_t precalculated_str = 0;
   Double_t *params;

   if (uparams) {
      params = const_cast<Double_t*>(uparams);
   } else {
      params = fParams;
   }
   UInt_t pos    = 0;
   UInt_t strpos = 0;

   for (i=0; i<fNoper; ++i) {

      const int oper = fOper[i];
      const int opcode = oper >> kTFOperShift;

      switch(opcode) {

         case kParameter  : { pos++; tab[pos-1] = params[ oper & kTFOperMask ]; continue; }
         case kConstant   : { pos++; tab[pos-1] = fConst[ oper & kTFOperMask ]; continue; }
         case kVariable   : { pos++; tab[pos-1] = x[ oper & kTFOperMask ]; continue; }
         case kStringConst: { strpos++; stringStack[strpos-1] = (char*)fExpr[i].Data(); pos++; tab[pos-1] = 0; continue; }

         case kAdd        : pos--; tab[pos-1] += tab[pos]; continue;
         case kSubstract  : pos--; tab[pos-1] -= tab[pos]; continue;
         case kMultiply   : pos--; tab[pos-1] *= tab[pos]; continue;
         case kDivide     : pos--; if (tab[pos] == 0) tab[pos-1] = 0; //  division by 0
                           else               tab[pos-1] /= tab[pos];
                           continue;
         case kModulo     : {pos--;
                              Long64_t int1((Long64_t)tab[pos-1]);
                              Long64_t int2((Long64_t)tab[pos]);
                              tab[pos-1] = Double_t(int1%int2);
                              continue;}

         case kcos  : tab[pos-1] = TMath::Cos(tab[pos-1]); continue;
         case ksin  : tab[pos-1] = TMath::Sin(tab[pos-1]); continue;
         case ktan  : if (TMath::Cos(tab[pos-1]) == 0) {tab[pos-1] = 0;} // { tangente indeterminee }
                        else tab[pos-1] = TMath::Tan(tab[pos-1]);
                        continue;
         case kacos : if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} //  indetermination
                        else tab[pos-1] = TMath::ACos(tab[pos-1]);
                        continue;
         case kasin : if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} //  indetermination
                        else tab[pos-1] = TMath::ASin(tab[pos-1]);
                        continue;
         case katan : tab[pos-1] = TMath::ATan(tab[pos-1]); continue;
         case kcosh : tab[pos-1] = TMath::CosH(tab[pos-1]); continue;
         case ksinh : tab[pos-1] = TMath::SinH(tab[pos-1]); continue;
         case ktanh : if (TMath::CosH(tab[pos-1]) == 0) {tab[pos-1] = 0;} // { tangente indeterminee }
                      else tab[pos-1] = TMath::TanH(tab[pos-1]);
                      continue;
         case kacosh: if (tab[pos-1] < 1) {tab[pos-1] = 0;} //  indetermination
                      else tab[pos-1] = TMath::ACosH(tab[pos-1]);
                      continue;
         case kasinh: tab[pos-1] = TMath::ASinH(tab[pos-1]); continue;
         case katanh: if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} // indetermination
                      else tab[pos-1] = TMath::ATanH(tab[pos-1]);
                      continue;
         case katan2: pos--; tab[pos-1] = TMath::ATan2(tab[pos-1],tab[pos]); continue;

         case kfmod : pos--; tab[pos-1] = fmod(tab[pos-1],tab[pos]); continue;
         case kpow  : pos--; tab[pos-1] = TMath::Power(tab[pos-1],tab[pos]); continue;
         case ksq   : tab[pos-1] = tab[pos-1]*tab[pos-1]; continue;
         case ksqrt : tab[pos-1] = TMath::Sqrt(TMath::Abs(tab[pos-1])); continue;

         case kstrstr : strpos -= 2; pos-=2; pos++;
                        if (strstr(stringStack[strpos],stringStack[strpos+1])) tab[pos-1]=1;
                        else tab[pos-1]=0;
                        continue;

         case kmin : pos--; tab[pos-1] = TMath::Min(tab[pos-1],tab[pos]); continue;
         case kmax : pos--; tab[pos-1] = TMath::Max(tab[pos-1],tab[pos]); continue;

         case klog  : if (tab[pos-1] > 0) tab[pos-1] = TMath::Log(tab[pos-1]);
                      else {tab[pos-1] = 0;} //{indetermination }
                      continue;
         case kexp  : { Double_t dexp = tab[pos-1];
                        if (dexp < -700) {tab[pos-1] = 0; continue;}
                        if (dexp >  700) {tab[pos-1] = TMath::Exp(700); continue;}
                        tab[pos-1] = TMath::Exp(dexp); continue;  }
         case klog10: if (tab[pos-1] > 0) tab[pos-1] = TMath::Log10(tab[pos-1]);
                      else {tab[pos-1] = 0;} //{indetermination }
                      continue;

         case kpi   : pos++; tab[pos-1] = TMath::ACos(-1); continue;

         case kabs  : tab[pos-1] = TMath::Abs(tab[pos-1]); continue;
         case ksign : if (tab[pos-1] < 0) tab[pos-1] = -1; else tab[pos-1] = 1; continue;
         case kint  : tab[pos-1] = Double_t(Int_t(tab[pos-1])); continue;

         case kSignInv: tab[pos-1] = -1 * tab[pos-1]; continue;

         case krndm : pos++; tab[pos-1] = gRandom->Rndm(); continue;

         case kAnd  : pos--; if (tab[pos-1]!=0 && tab[pos]!=0) tab[pos-1]=1;
                             else tab[pos-1]=0;
                      continue;
         case kOr   : pos--; if (tab[pos-1]!=0 || tab[pos]!=0) tab[pos-1]=1;
                             else tab[pos-1]=0;
                      continue;
         case kEqual: pos--; if (tab[pos-1] == tab[pos]) tab[pos-1]=1;
                             else tab[pos-1]=0;
                      continue;
         case kNotEqual : pos--; if (tab[pos-1] != tab[pos]) tab[pos-1]=1;
                                 else tab[pos-1]=0;
                          continue;
         case kLess     : pos--; if (tab[pos-1] < tab[pos]) tab[pos-1]=1;
                          else tab[pos-1]=0;
                          continue;
         case kGreater  : pos--; if (tab[pos-1] > tab[pos]) tab[pos-1]=1;
                          else tab[pos-1]=0;
                          continue;

         case kLessThan: pos--; if (tab[pos-1]<=tab[pos]) tab[pos-1]=1;
                         else tab[pos-1]=0;
                         continue;
         case kGreaterThan: pos--; if (tab[pos-1]>=tab[pos]) tab[pos-1]=1;
                            else tab[pos-1]=0;
                            continue;
         case kNot : if (tab[pos-1]!=0) tab[pos-1] = 0; else tab[pos-1] = 1;
                     continue;

         case kStringEqual : strpos -= 2; pos -=2 ; pos++;
                        if (!strcmp(stringStack[strpos+1],stringStack[strpos])) tab[pos-1]=1;
                        else tab[pos-1]=0;
                        continue;
         case kStringNotEqual: strpos -= 2; pos -= 2; pos++;
                        if (strcmp(stringStack[strpos+1],stringStack[strpos])) tab[pos-1]=1;
                        else tab[pos-1]=0;
                        continue;

         case kBitAnd : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) & ((Int_t) tab[pos]); continue;
         case kBitOr  : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) | ((Int_t) tab[pos]); continue;
         case kLeftShift : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) <<((Int_t) tab[pos]); continue;
         case kRightShift: pos--; tab[pos-1]= ((Int_t) tab[pos-1]) >>((Int_t) tab[pos]); continue;

         case kJump   : i = (oper & kTFOperMask); continue;
         case kJumpIf : pos--; if (!tab[pos]) i = (oper & kTFOperMask); continue;

         case kBoolOptimize: {
            // boolean operation optimizer

            int param = (oper & kTFOperMask);
            Bool_t skip = kFALSE;
            int op = param % 10; // 1 is && , 2 is ||

            if (op == 1 && (!tab[pos-1]) ) {
               // &&: skip the right part if the left part is already false

               skip = kTRUE;

               // Preserve the existing behavior (i.e. the result of a&&b is
               // either 0 or 1)
               tab[pos-1] = 0;

            } else if (op == 2 && tab[pos-1] ) {
               // ||: skip the right part if the left part is already true

               skip = kTRUE;

               // Preserve the existing behavior (i.e. the result of a||b is
               // either 0 or 1)
               tab[pos-1] = 1;

            }

            if (skip) {
               int toskip = param / 10;
               i += toskip;
            }
            continue;
         }

      }

      switch(opcode) {

         #define R__EXPO(var)                                                 \
         {                                                                    \
            pos++; int param = (oper & kTFOperMask);                          \
            tab[pos-1] = TMath::Exp(params[param]+params[param+1]*x[var]);  \
            continue;                                                         \
         }
         // case kexpo:
         case kxexpo: R__EXPO(0);
         case kyexpo: R__EXPO(1);
         case kzexpo: R__EXPO(2);
         case kxyexpo:{ pos++; int param = (oper & kTFOperMask);
                        tab[pos-1] = TMath::Exp(params[param]+params[param+1]*x[0]+params[param+2]*x[1]);
                        continue;  }

         #define R__GAUS(var)                                                    \
         {                                                                       \
            pos++; int param = (oper & kTFOperMask);                             \
            tab[pos-1] = params[param]*TMath::Gaus(x[var],params[param+1],params[param+2],IsNormalized()); \
            continue;                                                            \
         }

         // case kgaus:
         case kxgaus: R__GAUS(0);
         case kygaus: R__GAUS(1);
         case kzgaus: R__GAUS(2);
         case kxygaus: {pos++; int param = (oper & kTFOperMask);
                        Double_t intermede1;
                        if (params[param+2] == 0) {
                           intermede1=1e10;
                        } else {
                           intermede1=Double_t((x[0]-params[param+1])/params[param+2]);
                        }
                        Double_t intermede2;
                        if (params[param+4] == 0) {
                           intermede2=1e10;
                        } else {
                           intermede2=Double_t((x[1]-params[param+3])/params[param+4]);
                        }
                        tab[pos-1] = params[param]*TMath::Exp(-0.5*(intermede1*intermede1+intermede2*intermede2));
                        continue; }

         #define R__LANDAU(var)                                                                  \
         {                                                                                       \
            pos++; const int param = (oper & kTFOperMask);                                       \
            tab[pos-1] = params[param]*TMath::Landau(x[var],params[param+1],params[param+2],IsNormalized()); \
            continue;                                                                            \
         }
         // case klandau:
         case kxlandau: R__LANDAU(0);
         case kylandau: R__LANDAU(1);
         case kzlandau: R__LANDAU(2);
         case kxylandau: { pos++; int param = oper&0x7fffff /* ActionParams[i] */ ;
                           Double_t intermede1=TMath::Landau(x[0], params[param+1], params[param+2],IsNormalized());
                           Double_t intermede2=TMath::Landau(x[1], params[param+3], params[param+4],IsNormalized());
                           tab[pos-1] = params[param]*intermede1*intermede2;
                           continue;
         }

         #define R__POLY(var)                                                  \
         {                                                                     \
            pos++; int param = (oper & kTFOperMask);                           \
            tab[pos-1] = 0; Double_t intermede = 1;                            \
            Int_t inter = param/100; /* arrondit */                            \
            Int_t int1= param-inter*100-1; /* aucune simplification ! (sic) */ \
            for (j=0 ;j<inter+1;j++) {                                         \
               tab[pos-1] += intermede*params[j+int1];                        \
               intermede *= x[var];                                            \
            }                                                                  \
            continue;                                                          \
         }
         // case kpol:
         case kxpol: R__POLY(0);
         case kypol: R__POLY(1);
         case kzpol: R__POLY(2);

         case kDefinedVariable : {
            if (!precalculated) {
               precalculated = 1;
               for(j=0;j<fNval;j++) param_calc[j]=DefinedValue(j);
            }
            pos++; tab[pos-1] = param_calc[(oper & kTFOperMask)];
            continue;
         }

         case kDefinedString : {
            int param = (oper & kTFOperMask);
            if (!precalculated_str) {
               precalculated_str=1;
               for (j=0;j<fNstring;j++) string_calc[j]=DefinedString(j);
            }
            strpos++; stringStack[strpos-1] = string_calc[param];
            pos++; tab[pos-1] = 0;
            continue;
         }

         case kFunctionCall: {
            // an external function call

            int param = (oper & kTFOperMask);
            int fno   = param / 1000;
            int nargs = param % 1000;

            // Retrieve the function
            TMethodCall *method = (TMethodCall*)fFunctions.At(fno);

            // Set the arguments
            method->ResetParam();
            if (nargs) {
               UInt_t argloc = pos-nargs;
               for(j=0;j<nargs;j++,argloc++,pos--) {
                  method->SetParam(tab[argloc]);
               }
            }
            pos++;
            Double_t ret;
            method->Execute(ret);
            tab[pos-1] = ret; // check for the correct conversion!

            continue;
         };
      }
      if (!TestBit(kOptimizationError)) {
         SetBit(kOptimizationError);
         Warning("EvalParOld","Found an unsupported opcode (%d)",oper >> kTFOperShift);
      }
   }
   Double_t result0 = tab[0];
   return result0;

}

////////////////////////////////////////////////////////////////////////////////
/// Reconstruct the formula expression from the internal TFormula member variables
///
///   This function uses the internal member variables of TFormula to
///   construct the mathematical expression associated with the TFormula
///   instance. This function can be used to get an expanded version of the
///   expression originally assigned to the TFormula instance, i.e. that
///   the string returned by GetExpFormula() doesn't depend on other
///   TFormula object names.
///
///  if option contains "p" the returned string will contain the formula
///  expression with symbolic parameters, eg [0] replaced by the actual value
///  of the parameter. Example:
///  if expression in formula is: "[0]*(x>-[1])+[2]*exp(-[3]*x)"
///  and parameters are 3.25,-4.01,4.44,-0.04, GetExpFormula("p") will return:
///   "(3.25*(x>+4.01))+(4.44*exp(+0.04*x))"

TString TFormula::GetExpFormula(Option_t *option) const
{
   if (fNoper>0) {
      TString* tab=new TString[fNoper];
      Bool_t* ismulti=new Bool_t[fNoper];
      Int_t spos=0;

      ismulti[0]=kFALSE;
      Int_t optype;
      Int_t j;
      Int_t ternaryend = -1;
      for(Int_t i=0;i<fNoper;i++){
         optype = GetAction(i);

         if (ternaryend==i) {
            // The ? and : have been added to tab[spos-2]
            if(ismulti[spos-1]){
               tab[spos-2]=tab[spos-2]+"("+tab[spos-1]+")";
            } else {
               tab[spos-2]=tab[spos-2]+tab[spos-1];
            }
            spos--;
            // Do not call continue since we need to
            // do the rest of the loop.
         }

         // Boolean optimization breakpoint
         if (optype==kBoolOptimize) { // -3) {
            continue;
         }

         //Sign inversion
         if (optype==kSignInv) { // -1) {
            tab[spos-1]="-("+tab[spos-1]+")";
            // i++;
            continue;
         }

         //Simple name (parameter,pol0,landau, etc)
         if (kexpo<=optype && optype<=kzpol) { // >=0) {
            tab[spos]=fExpr[i];
            ismulti[spos]=kFALSE;
            spos++;
            continue;
         }

         //constants, variables x,y,z,t, pi
         if ((optype<=151 && optype>=140 && optype!=145) || (optype == 40)) {
            tab[spos]=fExpr[i];
            ismulti[spos]=kFALSE;
            spos++;
            continue;
         }

         //Basic operators (+,-,*,/,==,^,etc)
         if(((optype>0 && optype<6) || optype==20 ||
             (((optype>59 && optype<69) || (optype >75 && optype<82)) && spos>=2))) {
             // if(optype==-20 && spos>=2){
            if(ismulti[spos-2]){
               tab[spos-2]="("+tab[spos-2]+")";
            }
            if(ismulti[spos-1]){
               tab[spos-2]+=fExpr[i]+("("+tab[spos-1]+")");
            }else{
               tab[spos-2]+=fExpr[i]+tab[spos-1];
            }
            ismulti[spos-2]=kTRUE;
            spos--;
            continue;
         }
         //Ternary condition
         if (optype==kJumpIf) {
            if(ismulti[spos-1]){
               tab[spos-1]="("+tab[spos-1]+")?";
            } else {
               tab[spos-1]=tab[spos-1]+"?";
            }
            continue;
         }
         if (optype==kJump) {
            if(ismulti[spos-1]){
               tab[spos-2]=tab[spos-2]+"("+tab[spos-1]+"):";
            } else {
               tab[spos-2]=tab[spos-2]+tab[spos-1]+":";
            }
            ternaryend = GetActionParam(i) + 1;
            spos--;
            continue;
         }

         //Functions
         int offset = 0;
         TString funcname = fExpr[i];
         if((optype>9  && optype<16) ||
            (optype>20 && optype<23) ||
            (optype>29 && optype<34) ||
            (optype>40 && optype<44) ||
            (optype>69 && optype<76) ||
            (optype==145)) {
            //Functions with the format func(x)
            offset = -1;
         }

         if((optype>15 && optype<20) ||
            (optype>22 && optype<26)) {
            //Functions with the format func(x,y)
            offset = -2;
         }
         if(optype==145) {
            int param = (fOper[i] & kTFOperMask);
            //int fno   = param / 1000;
            int nargs = param % 1000;
            offset = -nargs;
            // The function name contains return type and parameters types we need
            // to trim them.
            int depth;
            for(j=0, depth=0;j<funcname.Length();++j) {
               switch (funcname[j]) {
                  case '<':
                     ++depth; break;
                  case '>':
                     --depth; break;
                  case ' ':
                     if (depth==0) {
                        funcname.Remove(0,j+1);
                        j = funcname.Length();
                        break;
                     }
               }
            }
            Ssiz_t ind = funcname.First('(');
            funcname.Remove(ind);
         }
         if (offset > 0) {
            Error("GetExpFormula","Internal error, number of argument found is %d",-offset);
         } else if (offset == 0) {
            tab[spos]=funcname+"()";
            ismulti[spos]=kFALSE;
            spos += 1;
            continue;
         } else if (offset<=0 && (spos+offset>=0)) {
            tab[spos+offset]=funcname+("("+tab[spos+offset]);
            for (j=offset+1; j<0; j++){
               tab[spos+offset]+=","+tab[spos+j];
            }
            tab[spos+offset]+=")";
            ismulti[spos+offset]=kFALSE;
            spos += offset+1;
            continue;
         }
      }
      if (ternaryend==fNoper) {
         // The ? and : have been added to tab[spos-2]
         if(ismulti[spos-1]){
            tab[spos-2]=tab[spos-2]+"("+tab[spos-1]+")";
         } else {
            tab[spos-2]=tab[spos-2]+tab[spos-1];
         }
         spos--;
      }

      TString ret = "";
      if (spos > 0) ret = tab[spos-1];
      delete[] tab;
      delete[] ismulti;

      //if option "p" is specified, return the real values of parameters instead of [0]
      TString opt = option;
      opt.ToLower();
      if (opt.Contains("p")) {
         char pb[13];
         char pbv[100];
         for (j=0;j<fNpar;j++) {
            snprintf(pb,sizeof(pb),"[%d]",j);
            snprintf(pbv,100,"%g",fParams[j]);
            ret.ReplaceAll(pb,pbv);
         }
         ret.ReplaceAll("--","+");
         ret.ReplaceAll("+-","-");
      }
      return ret;
   } else{
      TString ret="";
      return ret;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return linear part.

const TObject* TFormula::GetLinearPart(Int_t i)
{
   if (!fLinearParts.IsEmpty())
      return fLinearParts.UncheckedAt(i);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of parameter number ipar.

Double_t TFormula::GetParameter(Int_t ipar) const
{
   if (ipar <0 || ipar >= fNpar) return 0;
   return fParams[ipar];
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of parameter named parName.

Double_t TFormula::GetParameter(const char *parName) const
{
   const Double_t kNaN = 1e-300;
   Int_t index = GetParNumber(parName);
   if (index==-1) {
      Error("TFormula", "Parameter %s not found", parName);
      return kNaN;
   }
   return GetParameter(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Return name of one parameter.

const char *TFormula::GetParName(Int_t ipar) const
{
   if (ipar <0 || ipar >= fNpar) return "";
   if (fNames[ipar].Length() > 0) return (const char*)fNames[ipar];
   return Form("p%d",ipar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return parameter number by name.

Int_t TFormula::GetParNumber(const char *parName) const
{
   if (!parName)
      return -1;

   for (Int_t i=0; i<fNpar; i++) {
      if (!strcmp(GetParName(i),parName)) return i;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the expression at the index 'oper' has to be treated as a string

Bool_t TFormula::IsString(Int_t oper) const
{
   return GetAction(oper) == kStringConst;
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this formula with its attributes.

void TFormula::Print(Option_t *) const
{
   Int_t i;
   Printf(" %20s : %s Ndim= %d, Npar= %d, Noper= %d",GetName(),GetTitle(), fNdim,fNpar,fNoper);
   for (i=0;i<fNoper;i++) {
      Printf(" fExpr[%d] = %s  action = %d action param = %d ",
             i,(const char*)fExpr[i],GetAction(i),GetActionParam(i));
   }
   //MI change
   //
   if (fNOperOptimized>0){
      Printf("Optimized expression");
      for (i=0;i<fNOperOptimized;i++) {
         Printf(" fExpr[%d] = %s\t\t  action = %d action param = %d ",
            i,(const char*)fExprOptimized[i],GetActionOptimized(i),GetActionParamOptimized(i));
      }
   }

   if (!fNames) return;
   if (!fParams) return;
   for (i=0;i<fNpar;i++) {
      Printf(" Par%3d  %20s = %g",i,GetParName(i),fParams[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If the formula is for linear fitting, change the title to
/// normal and fill the LinearParts array

void TFormula::ProcessLinear(TString &formula)
{
   TString formula2(formula);
   char repl[20];
   char *pch;
   Int_t nf, offset, replsize;
   //replace "++" with "+[i]*"
   pch= (char*)strstr(formula.Data(), "++");
   if (pch)
      formula.Insert(0, "[0]*(");
   pch= (char*)strstr(formula.Data(), "++");
   if (pch){
      //if there are "++", replaces them with +[i]*
      nf = 1;
      while (pch){
         snprintf(repl,20, ")+[%d]*(", nf);
         offset = pch-formula.Data();
         if (nf<10) replsize = 7;
         else if (nf<100) replsize = 8;
         else replsize = 9;
         formula.Replace(pch-formula.Data(), 2, repl, replsize);
         pch = (char*)strstr(formula.Data()+offset, "++");
         nf++;
      }
      formula.Append(')', 1);
   } else {
      //if there are no ++, create a new string with ++ instead of +[i]*
      formula2=formula2(4, formula2.Length()-4);
      pch= (char*)strchr(formula2.Data(), '[');
      snprintf(repl,20, "++");
      nf = 1;
      while (pch){
         offset = pch-formula2.Data()-1;
         if (nf<10) replsize = 5;
         else replsize = 6;
         formula2.Replace(pch-formula2.Data()-1, replsize, repl, 2);
         pch = (char*)strchr(formula2.Data()+offset, '[');
         nf++;
      }
   }

   fLinearParts.Expand(nf);
   //break up the formula and fill the array of linear parts
   TString replaceformula;
   formula2 = formula2.ReplaceAll("++", 2, "|", 1);
   TObjArray *oa = formula2.Tokenize("|");
   TString replaceformula_name;
   for (Int_t i=0; i<nf; i++) {
      replaceformula = ((TObjString *)oa->UncheckedAt(i))->GetString();
      replaceformula_name = "f_linear_";
      replaceformula_name.Append(replaceformula);
      TFormula *f = new TFormula(replaceformula_name.Data(), replaceformula.Data());
      if (!f) {
         Error("TFormula", "f_linear not allocated");
         return;
      }
      {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfFunctions()->Remove(f);
      }
      f->SetBit(kNotGlobal, 1);
      fLinearParts.Add(f);
   }
   oa->Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize parameter number ipar.

void TFormula::SetParameter(const char *name, Double_t value)
{
   Int_t ipar = GetParNumber(name);
   if (ipar <0 || ipar >= fNpar) return;
   fParams[ipar] = value;
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize parameter number ipar.

void TFormula::SetParameter(Int_t ipar, Double_t value)
{
   if (ipar <0 || ipar >= fNpar) return;
   fParams[ipar] = value;
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize array of all parameters.
/// See also the next function with the same name.

void TFormula::SetParameters(const Double_t *params)
{
   for (Int_t i=0; i<fNpar;i++) {
      fParams[i] = params[i];
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize up to 11 parameters
/// All arguments except THE FIRST TWO are optional
/// In case of a function with only one parameter, call this function with p1=0.
/// Minimum two arguments are required to differentiate this function
/// from the SetParameters(cont Double_t *params)

void TFormula::SetParameters(Double_t p0,Double_t p1,Double_t p2,Double_t p3,Double_t p4
                       ,Double_t p5,Double_t p6,Double_t p7,Double_t p8,Double_t p9,Double_t p10)
{
   if (fNpar > 0) fParams[0] = p0;
   if (fNpar > 1) fParams[1] = p1;
   if (fNpar > 2) fParams[2] = p2;
   if (fNpar > 3) fParams[3] = p3;
   if (fNpar > 4) fParams[4] = p4;
   if (fNpar > 5) fParams[5] = p5;
   if (fNpar > 6) fParams[6] = p6;
   if (fNpar > 7) fParams[7] = p7;
   if (fNpar > 8) fParams[8] = p8;
   if (fNpar > 9) fParams[9] = p9;
   if (fNpar >10) fParams[10]= p10;
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Set name of parameter number ipar

void TFormula::SetParName(Int_t ipar, const char *name)
{
   if (ipar <0 || ipar >= fNpar) return;
   fNames[ipar] = name;
}

////////////////////////////////////////////////////////////////////////////////
/// Set up to 11 parameter names.

void TFormula::SetParNames(const char*name0,const char*name1,const char*name2,const char*name3,const char*name4,
                     const char*name5,const char*name6,const char*name7,const char*name8,const char*name9,const char*name10)
{
   if (fNpar > 0) fNames[0] = name0;
   if (fNpar > 1) fNames[1] = name1;
   if (fNpar > 2) fNames[2] = name2;
   if (fNpar > 3) fNames[3] = name3;
   if (fNpar > 4) fNames[4] = name4;
   if (fNpar > 5) fNames[5] = name5;
   if (fNpar > 6) fNames[6] = name6;
   if (fNpar > 7) fNames[7] = name7;
   if (fNpar > 8) fNames[8] = name8;
   if (fNpar > 9) fNames[9] = name9;
   if (fNpar >10) fNames[10]= name10;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TFormula::Streamer(TBuffer &b, const TClass *onfile_class)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v==6) {
         Error("Streamer","version 6 is not supported");
         return;
      }
      Streamer(b, v, R__s, R__c, onfile_class);

   } else {
      b.WriteClassBuffer(TFormula::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TFormula::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v==6) {
         Error("Streamer","version 6 is not supported");
         return;
      }
      Streamer(b, v, R__s, R__c, nullptr);

   } else {
      b.WriteClassBuffer(TFormula::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// specialized streamer function being able to read old TF1 versions as TF1Old in memory

void TFormula::Streamer(TBuffer &b, Int_t v, UInt_t R__s, UInt_t R__c, const TClass *onfile_class)
{

   //printf("Reading TFormula - version %d \n",v);
   if (v > 3 ) {
      b.ReadClassBuffer(TFormula::Class(), this, v, R__s, R__c, onfile_class);
      if (!TestBit(kNotGlobal)) {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfFunctions()->Add(this);
      }

      // We need to reinstate (if possible) the TMethodCall.
      if (fFunctions.GetLast()>=0) {
         // Compiles will reset the parameter values so we need
         // to temporarily keep them
         Double_t *param = fParams;
         TString *names = fNames;
         Int_t npar = fNpar;
         fParams = 0;
         fNames = 0;
         if (Compile()) {
            Error("Streamer","error compiling formula");
            return;
         }
         for (Int_t i = 0; i<npar && i<fNpar; ++i) fParams[i] = param[i];
         delete [] param;
         delete [] fNames;
         fNames = names;
      } else if (v<6) {
         Convert(v);
      }
      Optimize();
      return;
   }
   // version smaller or equal to 3
   // process old versions before automatic schema evolution
   TNamed::Streamer(b);
   b >> fNdim;
   b >> fNumber;
   if (v > 1) b >> fNval;
   if (v > 2) b >> fNstring;
   fNpar   = b.ReadArray(fParams);
   fOper = new Int_t[gMAXOP];
   fNoper  = b.ReadArray(fOper);
   fNconst = b.ReadArray(fConst);
   if (fNoper) {
      fExpr   = new TString[fNoper];
   }
   if (fNpar) {
      fNames  = new TString[fNpar];
   }
   Int_t i;
   for (i=0;i<fNoper;i++)  fExpr[i].Streamer(b);
   for (i=0;i<fNpar;i++)   fNames[i].Streamer(b);
   {
      R__LOCKGUARD(gROOTMutex);
      if (gROOT->GetListOfFunctions()->FindObject(GetName())) return;
      gROOT->GetListOfFunctions()->Add(this);
   }
   b.CheckByteCount(R__s, R__c, TFormula::IsA());

   Convert(v);
   // end of old versions

}

void TFormula::Convert(UInt_t /* fromVersion */)
{
   // Convert the fOper of a TFormula version fromVersion to the current in memory version

   enum {
      kOldexpo         =  1000,
      kOldgaus         =  2000,
      kOldlandau       =  4000,
      kOldxylandau     =  4500,
      kOldConstants    =  50000,
      kOldStrings      =  80000,
      kOldVariable     = 100000,
      kOldTreeString   = 105000,
      kOldFormulaVar   = 110000,
      kOldBoolOptimize = 120000,
      kOldFunctionCall = 200000
   };
   int i,j;

   for (i=0,j=0; i<fNoper; ++i,++j) {
      Int_t action = fOper[i];
      Int_t newActionCode = 0;
      Int_t newActionParam = 0;

      if ( action == 0) {
         // Sign Inversion

         newActionCode = kSignInv;

         Float_t aresult = 99.99;
         sscanf((const char*)fExpr[i],"%g",&aresult);
         R__ASSERT((aresult+1)<0.001);

         ++i; // skip the implied multiplication.

         // For consistency and for Optimize to work correctly
         // we need to remove the "-1" string in fExpr
         for (int z=i; z<fNoper; ++z) {
            fExpr[z-1] = fExpr[z];
         }

      } else  if ( action < 100 ) {
         // basic operators and mathematical library

         newActionCode = action;

      } else if (action >= kOldFunctionCall) {
         // Function call

         newActionCode = kFunctionCall;
         newActionParam = action-kOldFunctionCall;

      } else if (action >= kOldBoolOptimize) {
         // boolean operation optimizer

         newActionCode = kBoolOptimize;
         newActionParam = action-kOldBoolOptimize;

      } else if (action >= kOldFormulaVar) {
         // a variable

         newActionCode = kVariable;
         newActionParam = action-kOldFormulaVar;

      } else if (action >= kOldTreeString) {
         // a tree string

         newActionCode = kDefinedString;
         newActionParam = action-kOldTreeString;

      } else if (action >= kOldVariable) {
         // a tree variable

         newActionCode = kDefinedVariable;
         newActionParam = action-kOldVariable;

      } else if (action == kOldStrings) {
         // String

         newActionCode = kStringConst;

      } else if (action >= kOldConstants) {
         // numerical value

         newActionCode = kConstant;
         newActionParam = action-kOldConstants;

      } else if (action > 10000 && action < kOldConstants) {
         // Polynomial

         int var = action/10000; //arrondit
         newActionCode = kpol + (var-1);
         newActionParam = action - var*10000;

      } else if (action >= 4600) {

         Error("Convert","Unsupported value %d",action);

      } else if (action > kOldxylandau) {
         // xylandau

         newActionCode = kxylandau;
         newActionParam = action - (kOldxylandau+1);

      } else if (action > kOldlandau) {
         // landau, xlandau, ylandau or zlandau

         newActionCode = klandau;
         int var = action/100-40;
         if (var) newActionCode += var;
         newActionParam = action - var*100 - (kOldlandau+1);

      } else if (action > 2500 && action < 2600) {
         // xygaus

         newActionCode = kxygaus;
         newActionParam = action-2501;

      }  else if (action > 2000 && action < 2500) {
         //  gaus, xgaus, ygaus or zgaus

         newActionCode = kgaus;
         int var = action/100-20;
         if (var) newActionCode += var;
         newActionParam = action - var*100 - (kOldgaus+1);

      } else if (action > 1500 && action < 1600) {
         // xyexpo

         newActionCode = kxyexpo;
         newActionParam = action-1501;

      } else if (action > 1000 && action < 1500) {
         // expo or xexpo or yexpo or zexpo

         newActionCode = kexpo;
         int var = action/100-10;
         if (var) newActionCode += var;
         newActionParam = action - var*100 - (kOldexpo+1);

      } else if (action > 100 && action < 200) {
         // Parameter substitution

         newActionCode = kParameter;
         newActionParam = action - 101;
      }

      SetAction( j, newActionCode, newActionParam );

   }
   if (i!=j) {
      fNoper -= (i-j);
   }

}

////////////////////////////////////////////////////////////////////////////////
///  TOper offset  - helper class for TFormula*
///                     specify type of operand
///                     fTypeX   = kVariable
///                              = kParameter
///                              = kConstant
///                     fOffestX = offset in corresponding array

TOperOffset::TOperOffset()
{
   fType0=0;
   fType1=0;
   fType2=0;
   fType3=0;
   fOffset0=0;
   fOffset1=0;
   fOffset2=0;
   fOffset3=0;
   fOldAction=0;
   fToJump=0;
}

////////////////////////////////////////////////////////////////////////////////
///  MakePrimitive
///  find TFormulaPrimitive replacement for some operands

void  TFormula::MakePrimitive(const char *expr, Int_t pos)
{
   TString cbase(expr);
   cbase.ReplaceAll("Double_t ","");
   int paran = cbase.First("(");
   // int nargs = 0;
   if (paran>0) {
      //nargs = 1;
      cbase[paran]=0;
   }

   if (cbase=="<") cbase="XlY";
   if (cbase=="<=") cbase="XleY";
   if (cbase==">") cbase="XgY";
   if (cbase==">=") cbase="XgeY";
   if (cbase=="==" && GetActionOptimized(pos)!=kStringEqual) cbase="XeY";
   if (cbase=="!=" && GetActionOptimized(pos)!=kStringNotEqual) cbase="XneY";

   ROOT::v5::TFormulaPrimitive *prim = ROOT::v5::TFormulaPrimitive::FindFormula(cbase ,paran>0 ? cbase.Data() + paran + 1 : (const char*)0);
   if (prim) {
      fPredefined[pos] = prim;
      if (prim->fType==10) {
         SetActionOptimized(pos, kFD1);
      }
      if (prim->fType==110) {
         SetActionOptimized(pos, kFD2);
      }
      if (prim->fType==1110) {
         SetActionOptimized(pos, kFD3);
      }
      if (prim->fType==-1) {
         SetActionOptimized(pos, kFDM);
      }
      if (prim->fType==0){
         SetActionOptimized(pos,kConstant,fNconst);
         fConst[fNconst] = prim->Eval(0);
         fNconst++;
      }
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// MI include
///
/// Optimize formula
///  -  Minimize the number of operands
///     1.  several operands are glued together
///     2.  some primitive functions glued together - exemp. (x+y) => PlusXY(x,y)
///     3.  maximize number of standard calls minimizing number of jumps in Eval cases
///     4.  variables, parameters and constants are mapped - using fOperOfssets0
///          Eval procedure use direct acces to data (only one corresponding case statement in eval procedure)
/// ~~~ {.cpp}
///          pdata[operand={Var,Par,Const}][offset]
///          pdata[fOperOffsets0[i]][fOperOffset1[i+1]]
/// ~~~
///  - The fastest evaluation function is chosen at the end
///     1. fOptimal := pointer to the fastest function for given evaluation string
/// ~~~ {.cpp}
///             switch(GetActionOptimized(0)){
///               case kData : {fOptimal= (TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive0; break;}
///               case kUnary : {fOptimal= (TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive1; break;}
///               case kBinary : {fOptimal= (TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive2; break;}
///               case kThree : {fOptimal= (TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive3; break;}
///               case kFDM : {fOptimal= (TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive4; break;}
///             }
/// ~~~
///     2. ex.
///           -  fOptimal = ::EvalPrimitive0 - if it return only variable, constant or parameter
///           -           = ::EvalParameter1 - if only one unary operation
///           -           = ::EvalPrimitive2 - if only one binary operation

void TFormula::Optimize()
{
   //
   // Initialize data members
   //

   Int_t i;

   if (fPredefined)    { delete [] fPredefined;    fPredefined    = 0;}
   if (fOperOffset)    { delete [] fOperOffset;    fOperOffset    = 0;}
   if (fExprOptimized) { delete [] fExprOptimized; fExprOptimized = 0;}
   if (fOperOptimized) { delete [] fOperOptimized; fOperOptimized = 0;}

   fExprOptimized   = new TString[fNoper];
   fOperOptimized   = new Int_t[fNoper];
   fPredefined      = new ROOT::v5::TFormulaPrimitive*[fNoper];
   fOperOffset      = new TOperOffset[fNoper];
   for (i=0; i<fNoper; i++) {
      fExprOptimized[i]   = fExpr[i] ;
      fOperOptimized[i]   = fOper[i];
      fPredefined[i]= 0;
   }

   //
   //Make primitives
   //
   for (i=0;i<fNoper;i++){
      if (fExprOptimized[i].Data()) {
         MakePrimitive(fExprOptimized[i].Data(), i);
      }
   }
   //
   Int_t maxfound   = fNoper+1;
   Int_t *offset    = new Int_t[maxfound*16];
   Int_t *optimized = new Int_t[maxfound];
   //
   //
   ROOT::v5::TFormulaPrimitive*  primitive[10];
   primitive[0] =  ROOT::v5::TFormulaPrimitive::FindFormula("PlusXY");
   primitive[1] =  ROOT::v5::TFormulaPrimitive::FindFormula("MinusXY");
   primitive[2] =  ROOT::v5::TFormulaPrimitive::FindFormula("MultXY");
   primitive[3] =  ROOT::v5::TFormulaPrimitive::FindFormula("DivXY");
   primitive[4] =  ROOT::v5::TFormulaPrimitive::FindFormula("XpYpZ");
   primitive[5] =  ROOT::v5::TFormulaPrimitive::FindFormula("XxYxZ");
   primitive[6] =  ROOT::v5::TFormulaPrimitive::FindFormula("XxYpZ");
   primitive[7] =  ROOT::v5::TFormulaPrimitive::FindFormula("XpYxZ");
   primitive[8] =  ROOT::v5::TFormulaPrimitive::FindFormula("Pow2");
   primitive[9] =  ROOT::v5::TFormulaPrimitive::FindFormula("Pow3");
   //
   // set data pointers
   //
   for (i=0;i<fNoper;i++) optimized[i]=0;
   //
   for (i=0;i<fNoper;i++){
      Int_t actionparam = GetActionParamOptimized(i);
      Int_t action = GetActionOptimized(i);

      if (action==kBoolOptimize){
         //
         // optimize booleans
         //
         fOperOffset[i].fType1     = actionparam/10;           //  operands to skip
         fOperOffset[i].fOffset0   = actionparam%10;           //  1 is && , 2 is ||   - operand
         fOperOffset[i].fToJump    = i+fOperOffset[i].fType1;  //  where we should  jump
         continue;
      }
      if (action==kJump || action==kJumpIf) {
         // Ternary conditional operator
         fOperOffset[i].fType1     = action;
         fOperOffset[i].fToJump    = actionparam;
      }
      //
      if (action==kConstant&&i<fNoper-2){
         //
         // get offsets for kFDM operands
         //
         if (GetActionOptimized(i+1)==kConstant && GetActionOptimized(i+2)==kFDM){
            optimized[i]=1;
            optimized[i+1]=1;
            i+=2;
            fOperOffset[i].fType0=actionparam;
            fOperOffset[i].fOffset0=GetActionParamOptimized(i-1);
            Int_t offset2 =  int(fConst[fOperOffset[i].fOffset0]+0.4);
            fOperOffset[i].fOffset0=offset2;
            Int_t nparmax = offset2+fPredefined[i]->fNParameters;
            if (nparmax>fNpar){ // increase expected number of parameters
               fNpar=nparmax;
            }
            continue;
         }
      }
      switch(action){
      case kVariable : {action=kData; fOperOffset[i].fType0=0; break;}
      case kParameter: {action=kData; fOperOffset[i].fType0=1; break;}
      case kConstant : {action=kData; fOperOffset[i].fType0=2; break;}
      }
      //
      fOperOffset[i].fOffset0 = GetActionParamOptimized(i);
      SetActionOptimized(i,action, actionparam);    //set common data option
   }
   //
   //
   fNOperOptimized = fNoper;
   //
   for (i=0; i<fNoper; ++i)
   {
      //
      if (!(GetActionOptimized(i)== kData)) continue;
      offset[0] = fOperOffset[i].fType0;       //
      offset[1] = fOperOptimized[i] & kTFOperMask;   // offset

      if ((i+1) >= fNoper) continue;

      if (GetActionOptimized(i+1)==kFD1){
         optimized[i] = 1; // to be optimized
         i++;
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         SetActionOptimized(i  ,kUnary);
         continue;
      }
      if (GetActionOptimized(i+1)==kAdd){
         optimized[i] = 1; // to be optimized
         i++;
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         SetActionOptimized(i  ,kPlusD);
         continue;
      }
      if (GetActionOptimized(i+1)==kMultiply){
         optimized[i] = 1; // to be optimized
         i++;
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         SetActionOptimized(i,kMultD);
         continue;
      }

      if ((i+2) >= fNoper) continue;

      //
      //Binary operators
      if (!(GetActionOptimized(i+1)== kData))  continue;
      offset[2] = fOperOffset[i+1].fType0;
      offset[3] = fOperOptimized[i+1] & kTFOperMask;   // offset
      //
      if (GetActionOptimized(i+2)==kFD2 || GetActionOptimized(i+2)==kAdd ||GetActionOptimized(i+2)==kSubstract||
         GetActionOptimized(i+2)==kMultiply || GetActionOptimized(i+2)==kDivide){

         optimized[i] = 1; // to be optimized
         optimized[i+1] = 1; // to be optimized
         i+=2;
         //
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         fOperOffset[i].fType1   = offset[2];
         fOperOffset[i].fOffset1 = offset[3];
         fOperOffset[i].fType2   = GetActionOptimized(i);  //remember old action
         if (GetActionOptimized(i)==kAdd)       {fPredefined[i] = primitive[0];}
         if (GetActionOptimized(i)==kSubstract) {fPredefined[i] = primitive[1];}
         if (GetActionOptimized(i)==kMultiply) {
            fPredefined[i]=primitive[2];
            if (offset[0]==offset[2]&&offset[1]==offset[3]) {
               fPredefined[i] = primitive[8];
               SetActionOptimized(i,kUnary);
               continue;
            }
         }
         if (GetActionOptimized(i)==kDivide) {
            fPredefined[i] = primitive[3];
         }
         SetActionOptimized(i,kBinary);
         continue;
      }

      if ((i+3) >= fNoper) continue;

      //
      //operator 3
      //
      if (!(GetActionOptimized(i+2)== kData))  continue;
      offset[4] = fOperOffset[i+2].fType0;
      offset[5] = fOperOptimized[i+2] & kTFOperMask;   // offset
      //
      if (GetActionOptimized(i+3)==kFD3|| (  (GetActionOptimized(i+3)==kAdd||GetActionOptimized(i+3)==kMultiply) &&
         (GetActionOptimized(i+4)==kAdd||GetActionOptimized(i+4)==kMultiply) ) ){
         optimized[i+0]   = 1; // to be optimized
         optimized[i+1] = 1; // to be optimized
         optimized[i+2] = 1; // to be optimized
         i+=3;
         //
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         fOperOffset[i].fType1   = offset[2];
         fOperOffset[i].fOffset1 = offset[3];
         fOperOffset[i].fType2   = offset[4];
         fOperOffset[i].fOffset2 = offset[5];
         //
         fOperOffset[i].fOldAction = GetActionOptimized(i);  //remember old action
         if (GetActionOptimized(i)==kFD3) {
            SetActionOptimized(i,kThree);
            continue;
         }
         Int_t action=0;
         Int_t action2=kThree;
         if (GetActionOptimized(i)==kAdd&&GetActionOptimized(i+1)==kAdd)      action=4;
         if (GetActionOptimized(i)==kMultiply&&GetActionOptimized(i+1)==kMultiply) {
            action=5;
            if (offset[0]==offset[2]&&offset[1]==offset[3]&&offset[0]==offset[4]&&offset[1]==offset[5]){
               fPredefined[i]=primitive[9];
               action2=kUnary;
               action =9;
            }
         }
         if (GetActionOptimized(i)==kAdd&&GetActionOptimized(i+1)==kMultiply) action=6;
         if (GetActionOptimized(i)==kMultiply&&GetActionOptimized(i+1)==kAdd) action=7;
         //
         optimized[i]=1;
         i++;
         fOperOffset[i].fType0   = offset[0];
         fOperOffset[i].fOffset0 = offset[1];
         fOperOffset[i].fType1 = offset[2];
         fOperOffset[i].fOffset1 = offset[3];
         fOperOffset[i].fType2 = offset[4];
         fOperOffset[i].fOffset2 = offset[5];
         fPredefined[i]=primitive[action];
         SetActionOptimized(i,action2);
         continue;
      }
   }
   //
   //
   Int_t operO=0;
   TString expr="";
   Int_t *map0 = new Int_t[maxfound];   //remapping of the operands
   Int_t *map1 = new Int_t[maxfound];   //remapping of the operands
   for (i=0;i<fNoper;i++){
      map0[i]     =  operO;
      map1[operO] =  i;
      fOperOptimized[operO] = fOperOptimized[i];
      fPredefined[operO]    = fPredefined[i];
      fOperOffset[operO]    = fOperOffset[i];
      expr += fExprOptimized[i];
      if (optimized[i]==0){
         fExprOptimized[operO] = expr;
         expr = "";
         operO++;
      }else{
         expr += ",";
      }
   }
   //
   // Recalculate long jump for Boolean optimize
   //
   for (i=0; i<fNOperOptimized; i++){
      Int_t optaction = GetActionOptimized(i);
      if (optaction==kBoolOptimize){
         Int_t oldpos = fOperOffset[i].fToJump;
         Int_t newpos = oldpos==fNoper ? fNOperOptimized : map0[oldpos];
         fOperOffset[i].fToJump = newpos;   // new position to jump
         Int_t actionop = GetActionParamOptimized(i) % 10;
         switch (actionop) {
            case 1: SetActionOptimized(i,kBoolOptimizeAnd,newpos);  break;
            case 2: SetActionOptimized(i,kBoolOptimizeOr,newpos); break;
         }
      } else if (optaction==kJump || optaction==kJumpIf) {
         Int_t oldpos = fOperOffset[i].fToJump;
         Int_t newpos = oldpos==fNoper ? fNOperOptimized : map0[oldpos];
         fOperOffset[i].fToJump = newpos;   // new position to jump
         SetActionOptimized(i,optaction,newpos);
      }
   }


   fNOperOptimized = operO;
   //
   fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalParFast;
   if (fNOperOptimized==1) {
      switch(GetActionOptimized(0)){
         case kData   : {fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive0; break;}
         case kUnary  : {fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive1; break;}
         case kBinary : {fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive2; break;}
         case kThree  : {fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive3; break;}
         case kFDM    : {fOptimal= (ROOT::v5::TFormulaPrimitive::TFuncG)&TFormula::EvalPrimitive4; break;}
      }
   }

   delete [] map1;
   delete [] map0;
   delete [] offset;
   delete [] optimized;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive(const Double_t *x, const Double_t *params)
{
   const Double_t  *pdata[3] = {x,(params!=0)?params:fParams, fConst};
   Double_t result = pdata[fOperOffset->fType0][fOperOffset->fOffset0];
   switch((fOperOptimized[0] >> kTFOperShift)) {
      case kData          : return result;
      case kUnary         : return (fPredefined[0]->fFunc10)(pdata[fOperOffset->fType0][fOperOffset->fOffset0]);
      case kBinary         :return (fPredefined[0]->fFunc110)(result,
                               pdata[fOperOffset->fType1][fOperOffset->fOffset1]);

      case kThree         :return (fPredefined[0]->fFunc1110)(result, pdata[fOperOffset->fType1][fOperOffset->fOffset1],
                              pdata[fOperOffset->fType2][fOperOffset->fOffset2]);
      case kFDM         : return (fPredefined[0]->fFuncG)((Double_t*)&x[fOperOffset->fType0],
                             (Double_t*)&params[fOperOffset->fOffset0]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive0(const Double_t *x, const Double_t *params)
{
   const Double_t  *pdata[3] = {x,(params!=0)?params:fParams, fConst};
   return  pdata[fOperOffset->fType0][fOperOffset->fOffset0];
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive1(const Double_t *x, const Double_t *params)
{
   const Double_t  *pdata[3] = {x,(params!=0)?params:fParams, fConst};
   return (fPredefined[0]->fFunc10)(pdata[fOperOffset->fType0][fOperOffset->fOffset0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive2(const Double_t *x, const Double_t *params)
{
   const Double_t  *pdata[3] = {x,(params!=0)?params:fParams, fConst};
   return (fPredefined[0]->fFunc110)(pdata[fOperOffset->fType0][fOperOffset->fOffset0],
      pdata[fOperOffset->fType1][fOperOffset->fOffset1]);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive3(const Double_t *x, const Double_t *params)
{
   const Double_t  *pdata[3] = {x,(params!=0)?params:fParams, fConst};
   return (fPredefined[0]->fFunc1110)(pdata[fOperOffset->fType0][fOperOffset->fOffset0], pdata[fOperOffset->fType1][fOperOffset->fOffset1],
      pdata[fOperOffset->fType2][fOperOffset->fOffset2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate primitive formula

Double_t TFormula::EvalPrimitive4(const Double_t *x, const Double_t *params)
{
   const Double_t *par = (params!=0)?params:fParams;
   return (fPredefined[0]->fFuncG)((Double_t*)&x[fOperOffset->fType0],
      (Double_t*)&par[fOperOffset->fOffset0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate this formula.
///
/// The current value of variables x,y,z,t is passed through the pointer x.
/// The parameters used will be the ones in the array params if params is given
/// otherwise parameters will be taken from the stored data members fParams
///
/// \image html TFormula_eval.png

Double_t TFormula::EvalParFast(const Double_t *x, const Double_t *uparams)
{
   const Double_t  *pdata[3] = {x,(uparams!=0)?uparams:fParams, fConst};
   //
   Int_t i,j;
   Double_t tab[kMAXFOUND] = {0};
   const char *stringStack[gMAXSTRINGFOUND] = {0};
   Double_t param_calc[kMAXFOUND];
   char *string_calc[gMAXSTRINGFOUND] = {0};
   Int_t precalculated = 0;
   Int_t precalculated_str = 0;

   Double_t *params;

   if (uparams) {
      //for (j=0;j<fNpar;j++) fParams[j] = params[j];
      params = const_cast<Double_t*>(uparams);
   } else {
      params = fParams;
   }

   //if (params) {
   // for (j=0;j<fNpar;j++) fParams[j] = params[j];
   //}
   UInt_t pos    = 0;
   UInt_t strpos = 0;
   //   for (i=0; i<fNoper; ++i) {
   for (i=0; i<fNOperOptimized; ++i) {
      //
      const int oper = fOperOptimized[i];
      const int opcode = oper >> kTFOperShift;

      switch(opcode) {  // FREQUENTLY USED OPERATION
         case kData      : tab[pos] = pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0]; pos++;continue;
         case kPlusD     : tab[pos-1]+= pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0]; continue;
         case kMultD     : tab[pos-1]*= pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0]; continue;
         case kAdd       : pos--; tab[pos-1] += tab[pos]; continue;
         case kSubstract : pos--; tab[pos-1] -= tab[pos]; continue;
         case kMultiply  : pos--; tab[pos-1] *= tab[pos]; continue;
         case kDivide    : pos--; if (tab[pos] == 0) tab[pos-1] = 0; //  division by 0
                              else     tab[pos-1] /= tab[pos];
                              continue;
         case kUnary     : tab[pos] = (fPredefined[i]->fFunc10)(pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0]); pos++;continue;
         case kBinary    : tab[pos] = (fPredefined[i]->fFunc110)(pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0],
                              pdata[fOperOffset[i].fType1][fOperOffset[i].fOffset1]);pos++;continue;

         case kThree     : tab[pos]   = (fPredefined[i]->fFunc1110)(pdata[fOperOffset[i].fType0][fOperOffset[i].fOffset0],
                              pdata[fOperOffset[i].fType1][fOperOffset[i].fOffset1],
                              pdata[fOperOffset[i].fType2][fOperOffset[i].fOffset2]); pos++; continue;

         case kFDM       : tab[pos] = (fPredefined[i]->fFuncG)(&x[fOperOffset[i].fType0],&params[fOperOffset[i].fOffset0]); pos++;continue;
         case kFD1       : tab[pos-1]   =(fPredefined[i]->fFunc10)(tab[pos-1]); continue;
         case kFD2       :    pos--; tab[pos-1]   = (fPredefined[i]->fFunc110)(tab[pos-1],tab[pos]); continue;
         case kFD3       :    pos-=2; tab[pos-1]   = (fPredefined[i]->fFunc1110)(tab[pos-1],tab[pos],tab[pos+1]); continue;
      }
      //
      switch(opcode) {
         case kBoolOptimizeAnd: {
            if (!tab[pos-1]) i=fOperOffset[i].fToJump;
            continue;
         }
         case kBoolOptimizeOr: {
            if (tab[pos-1])  i=fOperOffset[i].fToJump;
            continue;
         }
         case kAnd  : pos--; tab[pos-1] = (bool)tab[pos];  continue;  // use the fact that other were check before - see bool optimize
         case kOr   : pos--; tab[pos-1] = (bool)tab[pos];  continue;
      }
      switch(opcode) {
         //    case kabs  : tab[pos-1] = TMath::Abs(tab[pos-1]); continue;
         case kabs  : if (tab[pos-1]<0) tab[pos-1]=-tab[pos-1]; continue;
         case ksign : if (tab[pos-1] < 0) tab[pos-1] = -1; else tab[pos-1] = 1; continue;
         case kint  : tab[pos-1] = Double_t(Int_t(tab[pos-1])); continue;
         case kpow  : pos--; tab[pos-1] = TMath::Power(tab[pos-1],tab[pos]); continue;

         case kModulo     : {pos--;
            Long64_t int1((Long64_t)tab[pos-1]);
            Long64_t int2((Long64_t)tab[pos]);
            tab[pos-1] = Double_t(int1%int2);
            continue;}


         case kStringConst: { strpos++; stringStack[strpos-1] = (char*)fExprOptimized[i].Data(); pos++; tab[pos-1] = 0; continue; }
         case kfmod : pos--; tab[pos-1] = fmod(tab[pos-1],tab[pos]); continue;

         case kstrstr : strpos -= 2; pos-=2; pos++;
            if (strstr(stringStack[strpos],stringStack[strpos+1])) tab[pos-1]=1;
            else tab[pos-1]=0;
            continue;
         case kpi   : pos++; tab[pos-1] = TMath::ACos(-1); continue;


         case kSignInv: tab[pos-1] = -1 * tab[pos-1]; continue;

         case krndm : pos++; tab[pos-1] = gRandom->Rndm(); continue;


         case kEqual: pos--; if (tab[pos-1] == tab[pos]) tab[pos-1]=1;
                      else tab[pos-1]=0;
                      continue;
         case kNotEqual : pos--; if (tab[pos-1] != tab[pos]) tab[pos-1]=1;
                          else tab[pos-1]=0;
                          continue;
         case kNot : if (tab[pos-1]!=0) tab[pos-1] = 0; else tab[pos-1] = 1;
                     continue;

         case kStringEqual : strpos -= 2; pos -=2 ; pos++;
            if (!strcmp(stringStack[strpos+1],stringStack[strpos])) tab[pos-1]=1;
            else tab[pos-1]=0;
            continue;
         case kStringNotEqual: strpos -= 2; pos -= 2; pos++;
            if (strcmp(stringStack[strpos+1],stringStack[strpos])) tab[pos-1]=1;
            else tab[pos-1]=0;
            continue;

         case kBitAnd : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) & ((Int_t) tab[pos]); continue;
         case kBitOr  : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) | ((Int_t) tab[pos]); continue;
         case kLeftShift : pos--; tab[pos-1]= ((Int_t) tab[pos-1]) <<((Int_t) tab[pos]); continue;
         case kRightShift: pos--; tab[pos-1]= ((Int_t) tab[pos-1]) >>((Int_t) tab[pos]); continue;

         case kJump   : i = (oper & kTFOperMask); continue;
         case kJumpIf : pos--; if (!tab[pos]) i = (oper & kTFOperMask); continue;

         case kBoolOptimize: {
            // boolean operation optimizer

            int param = (oper & kTFOperMask);
            int op = param % 10; // 1 is && , 2 is ||

            if (op == 1 && (!tab[pos-1]) ) {
               // &&: skip the right part if the left part is already false

               i +=  param / 10;

               // Preserve the existing behavior (i.e. the result of a&&b is
               // either 0 or 1)
               tab[pos-1] = 0;

            } else if (op == 2 && tab[pos-1] ) {
               // ||: skip the right part if the left part is already true

               i +=  param / 10;

               // Preserve the existing behavior (i.e. the result of a||b is
               // either 0 or 1)
               tab[pos-1] = 1;

            }

            continue;
         }

      }
      switch(opcode) {

#define R__EXPO(var)                                                         \
         {                                                                   \
           pos++; int param = (oper & kTFOperMask);                          \
           tab[pos-1] = TMath::Exp(params[param]+params[param+1]*x[var]);    \
           continue;                                                         \
         }
         // case kexpo:
         case kxexpo: R__EXPO(0);
         case kyexpo: R__EXPO(1);
         case kzexpo: R__EXPO(2);
         case kxyexpo:{  pos++; int param = (oper & kTFOperMask);
            tab[pos-1] = TMath::Exp(params[param]+params[param+1]*x[0]+params[param+2]*x[1]);
            continue;  }
#ifdef R__GAUS
#undef R__GAUS
#endif
#define R__GAUS(var)                                                                           \
                     {                                                                         \
                     pos++; int param = (oper & kTFOperMask);                                  \
                     tab[pos-1] = params[param]*TMath::Gaus(x[var],params[param+1],          \
                                                             params[param+2],IsNormalized()); \
                     continue;                                                                 \
                     }

                     // case kgaus:
         case kxgaus: R__GAUS(0);
         case kygaus: R__GAUS(1);
         case kzgaus: R__GAUS(2);
         case kxygaus: { pos++; int param = (oper & kTFOperMask);
            Double_t intermede1;
            if (params[param+2] == 0) {
               intermede1=1e10;
            } else {
               intermede1=Double_t((x[0]-params[param+1])/params[param+2]);
            }
            Double_t intermede2;
            if (params[param+4] == 0) {
               intermede2=1e10;
            } else {
               intermede2=Double_t((x[1]-params[param+3])/params[param+4]);
            }
            tab[pos-1] = params[param]*TMath::Exp(-0.5*(intermede1*intermede1+intermede2*intermede2));
            continue; }

#define R__LANDAU(var)                                                                  \
                      {                                                                                       \
                      pos++; const int param = (oper & kTFOperMask);                                       \
                      tab[pos-1] = params[param]*TMath::Landau(x[var],params[param+1],params[param+2],IsNormalized()); \
                      continue;                                                                            \
                      }
                      // case klandau:
         case kxlandau: R__LANDAU(0);
         case kylandau: R__LANDAU(1);
         case kzlandau: R__LANDAU(2);
         case kxylandau: { pos++; int param = oper&0x7fffff /* ActionParams[i] */ ;
            Double_t intermede1=TMath::Landau(x[0], params[param+1], params[param+2],IsNormalized());
            Double_t intermede2=TMath::Landau(x[1], params[param+3], params[param+4],IsNormalized());
            tab[pos-1] = params[param]*intermede1*intermede2;
            continue;
                        }

#define R__POLY(var)                                                                       \
                        {                                                                  \
                        pos++; int param = (oper & kTFOperMask);                           \
                        tab[pos-1] = 0; Double_t intermede = 1;                            \
                        Int_t inter = param/100; /* arrondit */                            \
                        Int_t int1= param-inter*100-1; /* aucune simplification ! (sic) */ \
                        for (j=0 ;j<inter+1;j++) {                                         \
                        tab[pos-1] += intermede*params[j+int1];                           \
                        intermede *= x[var];                                               \
                        }                                                                  \
                        continue;                                                          \
                        }
                        // case kpol:
         case kxpol: R__POLY(0);
         case kypol: R__POLY(1);
         case kzpol: R__POLY(2);

         case kDefinedVariable : {
            if (!precalculated) {
               precalculated = 1;
               for(j=0;j<fNval;j++) param_calc[j]=DefinedValue(j);
            }
            pos++; tab[pos-1] = param_calc[(oper & kTFOperMask)];
            continue;
            }

         case kDefinedString : {
            int param = (oper & kTFOperMask);
            if (!precalculated_str) {
               precalculated_str=1;
               for (j=0;j<fNstring;j++) string_calc[j]=DefinedString(j);
            }
            strpos++; stringStack[strpos-1] = string_calc[param];
            pos++; tab[pos-1] = 0;
            continue;
            }

         case kFunctionCall: {
            // an external function call

            int param = (oper & kTFOperMask);
            int fno   = param / 1000;
            int nargs = param % 1000;

            // Retrieve the function
            TMethodCall *method = (TMethodCall*)fFunctions.At(fno);

            // Set the arguments
            method->ResetParam();
            if (nargs) {
               UInt_t argloc = pos-nargs;
               for(j=0;j<nargs;j++,argloc++,pos--) {
                  method->SetParam(tab[argloc]);
               }
            }
            pos++;
            Double_t ret;
            method->Execute(ret);
            tab[pos-1] = ret; // check for the correct conversion!

            continue;
         };
      }
      if (!TestBit(kOptimizationError)) {
         SetBit(kOptimizationError);
         Warning("EvalParFast","Found an unsupported optmized opcode (%d)",oper >> kTFOperShift);
      }
   }
   Double_t result0 = tab[0];
   return result0;

}

////////////////////////////////////////////////////////////////////////////////
/// Pre compile function

Int_t TFormula::PreCompile()
{
   TString str = fTitle;
   if (str.Length()<3) return 1;
   if (str[str.Length()-1]!='+'&&str[str.Length()-2]!='+') return 1;
   str[str.Length()-2]=0;
   TString funName("preformula_");
   funName += fName;
   if (ROOT::v5::TFormulaPrimitive::FindFormula(funName)) return 0;
   TString fileName;
   fileName.Form("/tmp/%s.C",funName.Data());

   FILE *hf;
   hf = fopen(fileName.Data(),"w");
   if (hf == 0) {
      Error("PreCompile","Unable to open the file %s for writing.",fileName.Data());
      return 1;
   }
   fprintf(hf,   "/////////////////////////////////////////////////////////////////////////\n");
   fprintf(hf,   "//   This code has been automatically generated \n");
   //
   fprintf(hf,   "Double_t %s(Double_t *x, Double_t *p){",funName.Data());
   fprintf(hf,   "return (%s);\n}",str.Data());

   //   fprintf("ROOT::v5::TFormulaPrimitive::AddFormula(new ROOT::v5::TFormulaPrimitive(\"%s::%s\",\"%s::%s\",(ROOT::v5::TFormulaPrimitive::GenFunc0)%s::%s));\n",
   //             clname,method->GetName(),clname,method->GetName(),clname,method->GetName());
   fclose(hf);

   return 0;


}

////////////////////////////////////////////////////////////////////////////////
/// static function to set the maximum value of 3 parameters
///
///  - maxop    : maximum number of operations
///  - maxpar   : maximum number of parameters
///  - maxconst : maximum number of constants
///
/// None of these parameters cannot be less than 10 (default is 1000)
/// call this function to increase one or all maxima when processing
/// very complex formula, eg TFormula::SetMaxima(100000,1000,1000000);
/// If you process many functions with a small number of operations/parameters
/// you may gain some memory and performance by decreasing these values.

void TFormula::SetMaxima(Int_t maxop, Int_t maxpar, Int_t maxconst)
{
   gMAXOP    = TMath::Max(10,maxop);
   gMAXPAR   = TMath::Max(10,maxpar);
   gMAXCONST = TMath::Max(10,maxconst);
}

////////////////////////////////////////////////////////////////////////////////
/// static function to get the maximum value of 3 parameters
///  -maxop    : maximum number of operations
///  -maxpar   : maximum number of parameters
///  -maxconst : maximum number of constants

void TFormula::GetMaxima(Int_t& maxop, Int_t& maxpar, Int_t& maxconst)
{
   maxop = gMAXOP;
   maxpar = gMAXPAR;
   maxconst = gMAXCONST;
}

   } // end namespace v5

} // end namespace ROOT
