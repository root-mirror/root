// @(#)root/fumili:$Id$
// Author: Stanislav Nesterov  07/05/2003


/** \class TFumili

### FUMILI minimization package

FUMILI is based on ideas, proposed by I.N. Silin [See NIM A440, 2000 (p431)].
It was converted from FORTRAN to C  by Sergey Yaschenko <s.yaschenko@fz-juelich.de>


FUMILI is used to minimize Chi-square function or to search maximum of
likelihood function.

Experimentally measured values \f$F_i\f$ are fitted with theoretical
functions \f$f_i({\vec x}_i,\vec\theta\,\,)\f$, where \f${\vec x}_i\f$ are
coordinates, and \f$\vec\theta\f$ -- vector of parameters.

For better convergence Chi-square function has to be the following form

\f[
{\chi^2\over2}={1\over2}\sum^n_{i=1}\left(f_i(\vec
x_i,\vec\theta\,\,)-F_i\over\sigma_i\right)^2 \tag{1}
\f]

where \f$\sigma_i\f$ are errors of measured function.

The minimum condition is

\f[
{\partial\chi^2\over\partial\theta_i}=\sum^n_{j=1}{1\over\sigma^2_j}\cdot
{\partial f_j\over\partial\theta_i}\left[f_j(\vec
x_j,\vec\theta\,\,)-F_j\right]=0,\qquad i=1\ldots m\tag{2}
\f]

where m is the quantity of parameters.

Expanding left part of (2) over parameter increments and
retaining only linear terms one gets

\f[
\left(\partial\chi^2\over\theta_i\right)_{\vec\theta={\vec\theta}^0}
+\sum_k\left(\partial^2\chi^2\over\partial\theta_i\partial\theta_k\right)_{
\vec\theta={\vec\theta}^0}\cdot(\theta_k-\theta_k^0)
= 0\tag{3}
\f]

Here \f${\vec\theta}_0\f$ is some initial value of parameters. In general case:

\f[
{\partial^2\chi^2\over\partial\theta_i\partial\theta_k}=
\sum^n_{j=1}{1\over\sigma^2_j}{\partial f_j\over\theta_i}
{\partial f_j\over\theta_k} +
\sum^n_{j=1}{(f_j - F_j)\over\sigma^2_j}\cdot
{\partial^2f_j\over\partial\theta_i\partial\theta_k}\tag{4}
\f]

In FUMILI algorithm for second derivatives of Chi-square approximate
expression is used when last term in (4) is discarded. It is often
done, not always wittingly, and sometimes causes troubles, for example,
if user wants to limit parameters with positive values by writing down
\f$\theta_i^2\f$ instead of \f$\theta_i\f$. FUMILI will fail if one tries
minimize \f$\chi^2 = g^2(\vec\theta)\f$ where g is arbitrary function.

Approximate value is:
\f[{\partial^2\chi^2\over\partial\theta_i\partial\theta_k}\approx
Z_{ik}=
\sum^n_{j=1}{1\over\sigma^2_j}{\partial f_j\over\theta_i}
{\partial f_j\over\theta_k}\tag{5}
\f]

Then the equations for parameter increments are
\f[\left(\partial\chi^2\over\partial\theta_i\right)_{\vec\theta={\vec\theta}^0}
+\sum_k Z_{ik}\cdot(\theta_k-\theta^0_k) = 0,
\qquad i=1\ldots m\tag{6}
\f]

Remarkable feature of algorithm is the technique for step
restriction. For an initial value of parameter \f${\vec\theta}^0\f$ a
parallelepiped \f$P_0\f$ is built with the center at \f${\vec\theta}^0\f$ and
axes parallel to coordinate axes \f$\theta_i\f$. The lengths of
parallelepiped sides along i-th axis is \f$2b_i\f$, where \f$b_i\f$ is such a
value that the functions \f$f_j(\vec\theta)\f$ are quasi-linear all over
the parallelepiped.

FUMILI takes into account simple linear inequalities in the form:
\f[
\theta_i^{\rm min}\le\theta_i\le\theta^{\rm max}_i\tag{7}
\f]

They form parallelepiped \f$P\f$ (\f$P_0\f$ may be deformed by \f$P\f$).
Very similar step formulae are used in FUMILI for negative logarithm
of the likelihood function with the same idea - linearization of
function argument.

*/


#include "TFumili.h"

#include <iostream>
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TList.h"
#include "TVirtualFitter.h"


extern void H1FitChisquareFumili(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);
extern void H1FitLikelihoodFumili(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);
extern void GraphFitChisquareFumili(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);


ClassImp(TFumili);

TFumili *gFumili=0;
// Machine dependent values  fiXME!!
// But don't set min=max=0 if param is unlimited
static const Double_t gMAXDOUBLE=1e300;
static const Double_t gMINDOUBLE=-1e300;

////////////////////////////////////////////////////////////////////////////////

TFumili::TFumili(Int_t maxpar)
{//----------- FUMILI constructor ---------
   // maxpar is the maximum number of parameters used with TFumili object
   //
   fMaxParam = TMath::Max(maxpar,25);
   if (fMaxParam>200) fMaxParam=200;
   BuildArrays();

   fNumericDerivatives = true;
   fLogLike = false;
   fNpar    = fMaxParam;
   fGRAD    = false;
   fWARN    = true;
   fDEBUG   = false;
   fNlog    = 0;
   fSumLog  = 0;
   fNED1    = 0;
   fNED2    = 0;
   fNED12   = fNED1+fNED2;
   fEXDA    = 0;
   fFCN     = 0;
   fNfcn    = 0;
   fRP      = 1.e-15; //precision
   fS       = 1e10;
   fEPS     =0.01;
   fENDFLG  = 0;
   fNlimMul = 2;
   fNmaxIter= 150;
   fNstepDec= 3;
   fLastFixed = -1;

   fAKAPPA = 0.;
   fGT = 0.;
   for (int i = 0; i<5; ++i) fINDFLG[i] = 0;

   SetName("Fumili");
   gFumili = this;
   gROOT->GetListOfSpecials()->Add(gFumili);
}

////////////////////////////////////////////////////////////////////////////////
///
///   Allocates memory for internal arrays. Called by TFumili::TFumili
///

void TFumili::BuildArrays(){
   fCmPar      = new Double_t[fMaxParam];
   fA          = new Double_t[fMaxParam];
   fPL0        = new Double_t[fMaxParam];
   fPL         = new Double_t[fMaxParam];
   fParamError = new Double_t[fMaxParam];
   fDA         = new Double_t[fMaxParam];
   fAMX        = new Double_t[fMaxParam];
   fAMN        = new Double_t[fMaxParam];
   fR          = new Double_t[fMaxParam];
   fDF         = new Double_t[fMaxParam];
   fGr         = new Double_t[fMaxParam];
   fANames     = new TString[fMaxParam];

   //   fX = new Double_t[10];

   Int_t zSize = fMaxParam*(fMaxParam+1)/2;
   fZ0 = new Double_t[zSize];
   fZ  = new Double_t[zSize];

   for (Int_t i=0;i<fMaxParam;i++) {
      fA[i] =0.;
      fDF[i]=0.;
      fAMN[i]=gMINDOUBLE;
      fAMX[i]=gMAXDOUBLE;
      fPL0[i]=.1;
      fPL[i] =.1;
      fParamError[i]=0.;
      fANames[i]=Form("%d",i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TFumili destructor

TFumili::~TFumili() {
   DeleteArrays();
   if (gROOT && !gROOT->TestBit(TObject::kInvalidObject))
      gROOT->GetListOfSpecials()->Remove(this);
   if (gFumili == this) gFumili = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// return a chisquare equivalent

Double_t TFumili::Chisquare(Int_t npar, Double_t *params) const
{
   Double_t amin = 0;
   H1FitChisquareFumili(npar,params,amin,params,1);
   return 2*amin;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Resets all parameter names, values and errors to zero
///
/// Argument opt is ignored
///
/// NB: this procedure doesn't reset parameter limits

void TFumili::Clear(Option_t *)
{
   fNpar = fMaxParam;
   fNfcn = 0;
   for (Int_t i=0;i<fNpar;i++) {
      fA[i]   =0.;
      fDF[i]  =0.;
      fPL0[i] =.1;
      fPL[i]  =.1;
      fAMN[i] = gMINDOUBLE;
      fAMX[i] = gMAXDOUBLE;
      fParamError[i]=0.;
      fANames[i]=Form("%d",i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deallocates memory. Called from destructor TFumili::~TFumili

void TFumili::DeleteArrays(){
   delete[] fCmPar;
   delete[] fANames;
   delete[] fDF;
   // delete[] fX;
   delete[] fZ0;
   delete[] fZ;
   delete[] fGr;
   delete[] fA;
   delete[] fPL0;
   delete[] fPL;
   delete[] fDA;
   delete[] fAMN;
   delete[] fAMX;
   delete[] fParamError;
   delete[] fR;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Calculates partial derivatives of theoretical function
///
/// Input:
///   - fX  - vector of data point
///
/// Output:
///   - DF - array of derivatives
///
/// ARITHM.F: Converted from CERNLIB

void TFumili::Derivatives(Double_t *df,Double_t *fX){
   Double_t ff,ai,hi,y,pi;
   y = EvalTFN(df,fX);
   for (Int_t i=0;i<fNpar;i++) {
      df[i]=0;
      if(fPL0[i]>0.) {
         ai = fA[i]; // save current parameter value
         hi = 0.01*fPL0[i]; // diff step
         pi = fRP*TMath::Abs(ai);
         if (hi<pi) hi = pi; // if diff step is less than precision
         fA[i] = ai+hi;

         if (fA[i]>fAMX[i]) { // if param is out of limits
            fA[i] = ai-hi;
            hi = -hi;
            if (fA[i]<fAMN[i]) { // again out of bounds
               fA[i] = fAMX[i];   // set param to high limit
               hi = fAMX[i]-ai;
               if (fAMN[i]-ai+hi<0) { // if hi < (ai-fAMN)
                  fA[i]=fAMN[i];
                  hi=fAMN[i]-ai;
               }
            }
         }
         ff = EvalTFN(df,fX);
         df[i] = (ff-y)/hi;
         fA[i] = ai;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate the minimisation function
///
///  Input parameters:
///   - npar:    number of currently variable parameters
///   - par:     array of (constant and variable) parameters
///   - flag:    Indicates what is to be calculated
///   - grad:    array of gradients
///
///  Output parameters:
///   - fval:    The calculated function value.
///   - grad:    The vector of first derivatives.
///
/// The meaning of the parameters par is of course defined by the user,
/// who uses the values of those parameters to calculate their function value.
/// The starting values must be specified by the user.
///
/// Inside FCN user has to define Z-matrix by means TFumili::GetZ
///  and TFumili::Derivatives,
/// set theoretical function by means of TFumili::SetUserFunc,
/// but first - pass number of parameters by TFumili::SetParNumber
///
/// Later values are determined by Fumili as it searches for the minimum
/// or performs whatever analysis is requested by the user.
///
/// The default function calls the function specified in SetFCN

Int_t TFumili::Eval(Int_t& npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
   if (fFCN) (*fFCN)(npar,grad,fval,par,flag);
   return npar;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate theoretical function
/// - df: array of partial derivatives
/// - X:  vector of theoretical function argument

Double_t TFumili::EvalTFN(Double_t * /*df*/, Double_t *X)
{
   // for the time being disable possibility to compute derivatives
   //if(fTFN)
   //  return (*fTFN)(df,X,fA);
   //else if(fTFNF1) {

   TF1 *f1 = (TF1*)fUserFunc;
   return f1->EvalPar(X,fA);
   //}
   //return 0.;
}

////////////////////////////////////////////////////////////////////////////////
///
///  Execute MINUIT commands. MINImize, SIMplex, MIGrad and FUMili all
///  will call TFumili::Minimize method.
///
///  For full command list see
///  MINUIT. Reference Manual. CERN Program Library Long Writeup D506.
///
///  Improvement and errors calculation are not yet implemented as well
///  as Monte-Carlo seeking and minimization.
///  Contour commands are also unsupported.
///
/// - command   : command string
/// - args      : array of arguments
/// - nargs     : number of arguments

Int_t TFumili::ExecuteCommand(const char *command, Double_t *args, Int_t nargs){
   TString comand = command;
   static TString clower = "abcdefghijklmnopqrstuvwxyz";
   static TString cupper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
   const Int_t nntot = 40;
   const char *cname[nntot] = {
      "MINImize  ",    //  0    checked
      "SEEk      ",    //  1    none
      "SIMplex   ",    //  2    checked same as 0
      "MIGrad    ",    //  3    checked  same as 0
      "MINOs     ",    //  4    none
      "SET xxx   ",    //  5 lot of stuff
      "SHOw xxx  ",    //  6 -----------
      "TOP of pag",    //  7 .
      "fiX       ",   //   8 .
      "REStore   ",   //   9 .
      "RELease   ",   //   10 .
      "SCAn      ",   //   11  not yet implemented
      "CONtour   ",   //   12  not yet implemented
      "HESse     ",   //   13  not yet implemented
      "SAVe      ",   //   14  obsolete
      "IMProve   ",   //   15  not yet implemented
      "CALl fcn  ",   //   16 .
      "STAndard  ",   //   17 .
      "END       ",   //   18 .
      "EXIt      ",   //   19 .
      "RETurn    ",   //   20 .
      "CLEar     ",   //   21 .
      "HELP      ",   //   22 not yet implemented
      "MNContour ",   //   23 not yet implemented
      "STOp      ",   //   24 .
      "JUMp      ",   //   25 not yet implemented
      "          ",   //
      "          ",   //
      "FUMili    ",   //    28 checked same as 0
      "          ",   //
      "          ",  //
      "          ",  //
      "          ",  //
      "COVARIANCE",  // 33
      "PRINTOUT  ",  // 34
      "GRADIENT  ",  // 35
      "MATOUT    ",  // 36
      "ERROR DEF ",  // 37
      "LIMITS    ",  // 38
      "PUNCH     "}; // 39


   fCword = comand;
   fCword.ToUpper();
   if (nargs<=0) fCmPar[0] = 0;
   Int_t i;
   for(i=0;i<fMaxParam;i++) {
      if(i<nargs) fCmPar[i] = args[i];
   }
   /*
   fNmaxIter = int(fCmPar[0]);
   if (fNmaxIter <= 0) {
       fNmaxIter = fNpar*10 + 20 + fNpar*M*5;
   }
   fEPS = fCmPar[1];
   */
   //*-*-               look for command in list CNAME . . . . . . . . . .
   TString ctemp = fCword(0,3);
   Int_t ind;
   for (ind = 0; ind < nntot; ++ind) {
      if (strncmp(ctemp.Data(),cname[ind],3) == 0) break;
   }
   if (ind==nntot) return -3; // Unknown command - input ignored
   if (fCword(0,4) == "MINO") ind=3;
   switch (ind) {
      case 0:  case 3: case 2: case 28:
         // MINImize [maxcalls] [tolerance]
         // also SIMplex, MIGrad  and  FUMili
         if(nargs>=1)
         fNmaxIter=TMath::Max(Int_t(fCmPar[0]),fNmaxIter); // fiXME!!
         if(nargs==2)
            fEPS=fCmPar[1];
         return Minimize();
      case 1:
         // SEEk not implemented in this package
         return -10;

      case 4: // MINos errors analysis not implemented
         return -10;

      case 5: case 6: // SET xxx & SHOW xxx
         return ExecuteSetCommand(nargs);

      case 7: // Obsolete command
         Printf("1");
         return 0;
      case 8: // fiX <parno> ....
         if (nargs<1) return -1; // No parameters specified
         for (i=0;i<nargs;i++) {
            Int_t parnum = Int_t(fCmPar[i])-1;
            FixParameter(parnum);
         }
         return 0;
      case 9: // REStore <code>
         if (nargs<1) return 0;
         if(fCmPar[0]==0.)
            for (i=0;i<fNpar;i++)
               ReleaseParameter(i);
         else
            if(fCmPar[0]==1.) {
               ReleaseParameter(fLastFixed);
               std::cout <<fLastFixed<<std::endl;
            }
         return 0;
      case 10: // RELease <parno> ...
         if (nargs<1) return -1; // No parameters specified
         for (i=0;i<nargs;i++) {
            Int_t parnum = Int_t(fCmPar[i])-1;
            ReleaseParameter(parnum);
         }
         return 0;
      case 11: // SCAn not implemented
         return -10;
      case 12: // CONt not implemented
         return -10;

      case 13: // HESSe not implemented
         return -10;
      case 14: // SAVe
         Printf("SAVe command is obsolete");
         return -10;
      case 15: // IMProve not implemented
         return -10;
      case 16: // CALl fcn <iflag>
         {if(nargs<1) return -1;
         Int_t flag = Int_t(fCmPar[0]);
         Double_t fval;
         Eval(fNpar,fGr,fval,fA,flag);
         return 0;}
      case 17: // STAndard must call function STAND
         return 0;
      case 18:   case 19:
      case 20:  case 24: {
         Double_t fval;
         Int_t flag = 3;
         Eval(fNpar,fGr,fval,fA,flag);
         return 0;
      }
      case 21:
         Clear();
         return 0;
      case 22: //HELp not implemented
      case 23: //MNContour not implemented
      case 25: // JUMp not implemented
         return -10;
      case 26:   case 27:  case 29:  case 30:  case 31:  case 32:
         return 0; // blank commands
      case 33:   case 34:   case 35:  case 36:   case 37:  case 38:
      case 39:
         Printf("Obsolete command. Use corresponding SET command instead");
         return -10;
      default:
         break;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Called from TFumili::ExecuteCommand in case
/// of "SET xxx" and "SHOW xxx".

Int_t TFumili::ExecuteSetCommand(Int_t nargs){
   static Int_t nntot = 30;
   static const char *cname[30] = {
      "FCN value ", // 0 .
      "PARameters", // 1 .
      "LIMits    ", // 2 .
      "COVariance", // 3 .
      "CORrelatio", // 4 .
      "PRInt levl", // 5 not implemented yet
      "NOGradient", // 6 .
      "GRAdient  ", // 7 .
      "ERRor def ", // 8 not sure how to implement - by time being ignored
      "INPut file", // 9 not implemented
      "WIDth page", // 10 not implemented yet
      "LINes page", // 11 not implemented yet
      "NOWarnings", // 12 .
      "WARnings  ", // 13 .
      "RANdom gen", // 14 not implemented
      "TITle     ", // 15 ignored
      "STRategy  ", // 16 ignored
      "EIGenvalue", // 17 not implemented yet
      "PAGe throw", // 18 ignored
      "MINos errs", // 19 not implemented yet
      "EPSmachine", // 20 .
      "OUTputfile", // 21 not implemented
      "BATch     ", // 22 ignored
      "INTeractiv", // 23 ignored
      "VERsion   ", // 24 .
      "reserve   ", // 25 .
      "NODebug   ", // 26 .
      "DEBug     ", // 27 .
      "SHOw      ", // 28 err
      "SET       "};// 29 err

   TString  cfname, cmode, ckind,  cwarn, copt, ctemp, ctemp2;
   Int_t i, ind;
   Bool_t setCommand=kFALSE;
   for (ind = 0; ind < nntot; ++ind) {
      ctemp  = cname[ind];
      ckind  = ctemp(0,3);
      ctemp2 = fCword(4,6);
      if (strstr(ctemp2.Data(),ckind.Data())) break;
   }
   ctemp2 = fCword(0,3);
   if(ctemp2.Contains("SET")) setCommand=true;
   if(ctemp2.Contains("HEL") || ctemp2.Contains("SHO")) setCommand=false;

   if (ind>=nntot) return -3;

   switch(ind) {
      case 0: // SET FCN value illegal // SHOw only
         if(!setCommand) Printf("FCN=%f",fS);
         return 0;
      case 1: // PARameter <parno> <value>
         {
         if (nargs<2 && setCommand) return -1;
         Int_t parnum;
         Double_t val;
         if(setCommand) {
            parnum = Int_t(fCmPar[0])-1;
            val= fCmPar[1];
            if(parnum<0 || parnum>=fNpar) return -2; //no such parameter
            fA[parnum] = val;
         } else {
            if (nargs>0) {
               parnum = Int_t(fCmPar[0])-1;
               if(parnum<0 || parnum>=fNpar) return -2; //no such parameter
               Printf("Parameter %s = %E",fANames[parnum].Data(),fA[parnum]);
            } else
               for (i=0;i<fNpar;i++)
               Printf("Parameter %s = %E",fANames[i].Data(),fA[i]);

         }
         return 0;
         }
      case 2: // LIMits [parno] [ <lolim> <uplim> ]
         {
         Int_t parnum;
         Double_t lolim,uplim;
         if (nargs<1) {
            for(i=0;i<fNpar;i++)
               if(setCommand) {
                  fAMN[i] = gMINDOUBLE;
                  fAMX[i] = gMAXDOUBLE;
               } else
                  Printf("Limits for param %s: Low=%E, High=%E",
                         fANames[i].Data(),fAMN[i],fAMX[i]);
         } else {
            parnum = Int_t(fCmPar[0])-1;
            if(parnum<0 || parnum>=fNpar)return -1;
            if(setCommand) {
               if(nargs>2) {
                  lolim = fCmPar[1];
                  uplim = fCmPar[2];
                  if(uplim==lolim) return -1;
                  if(lolim>uplim) {
                     Double_t tmp = lolim;
                     lolim = uplim;
                     uplim = tmp;
                  }
               } else {
                  lolim = gMINDOUBLE;
                  uplim = gMAXDOUBLE;
               }
               fAMN[parnum] = lolim;
               fAMX[parnum] = uplim;
            } else
               Printf("Limits for param %s Low=%E, High=%E",
                    fANames[parnum].Data(),fAMN[parnum],fAMX[parnum]);
            }
            return 0;
         }
      case 3:
         {
         if(setCommand) return 0;
         Printf("\nCovariant matrix ");
         Int_t l = 0,nn=0,nnn=0;
         for (i=0;i<fNpar;i++) if(fPL0[i]>0.) nn++;
         for (i=0;i<nn;i++) {
            for(;fPL0[nnn]<=0.;nnn++) { }
            printf("%5s: ",fANames[nnn++].Data());
            for (Int_t j=0;j<=i;j++)
               printf("%11.2E",fZ[l++]);
            std::cout<<std::endl;
         }
         std::cout<<std::endl;
         return 0;
         }
      case 4:
         if(setCommand) return 0;
         Printf("\nGlobal correlation factors (maximum correlation of the parameter\n  with arbitrary linear combination of other parameters)");
         for(i=0;i<fNpar;i++) {
            printf("%5s: ",fANames[i].Data());
            printf("%11.3E\n",TMath::Sqrt(1-1/((fR[i]!=0.)?fR[i]:1.)) );
         }
         std::cout<<std::endl;
         return 0;
      case 5:   // PRIntout not implemented
         return -10;
      case 6: // NOGradient
         if(!setCommand) return 0;
         fGRAD = false;
         return 0;
      case 7: // GRAdient
         if(!setCommand) return 0;
         fGRAD = true;
         return 0;
      case 8: // ERRordef - now ignored
         return 0;
      case 9: // INPut - not implemented
         return -10;
      case 10: // WIDthpage - not implemented
         return -10;
      case 11: // LINesperpage - not implemented
         return -10;
      case 12: //NOWarnings
         if(!setCommand) return 0;
         fWARN = false;
         return 0;
      case 13: // WARnings
         if(!setCommand) return 0;
         fWARN = true;
         return 0;
      case 14: // RANdomgenerator - not implemented
         return -10;
      case 15: // TITle - ignored
         return 0;
      case 16: // STRategy - ignored
         return 0;
      case 17: // EIGenvalues - not implemented
         return -10;
      case 18: // PAGethrow - ignored
         return 0;
      case 19: // MINos errors - not implemented
         return -10;
      case 20: //EPSmachine
         if(!setCommand) {
            Printf("Relative floating point precision RP=%E",fRP);
         } else
            if (nargs>0) {
               Double_t pres=fCmPar[0];
               if (pres<1e-5 && pres>1e-34) fRP=pres;
            }
         return 0;
      case 21: // OUTputfile - not implemented
         return -10;
      case 22: // BATch - ignored
         return 0;
      case 23: // INTerative - ignored
         return 0;
      case 24: // VERsion
         if(setCommand) return 0;
         Printf("FUMILI-ROOT version 0.1");
         return 0;
      case 25: // reserved
         return 0;
      case 26: // NODebug
         if(!setCommand) return 0;
         fDEBUG = false;
         return 0;
      case 27: // DEBug
         if(!setCommand) return 0;
         fDEBUG = true;
         return 0;
      case 28:
      case 29:
         return -3;
      default:
         break;
   }
   return -3;
}

////////////////////////////////////////////////////////////////////////////////
/// Fixes parameter number ipar

void TFumili::FixParameter(Int_t ipar) {
   if(ipar>=0 && ipar<fNpar && fPL0[ipar]>0.) {
      fPL0[ipar] = -fPL0[ipar];
      fLastFixed = ipar;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the covariance matrix

Double_t *TFumili::GetCovarianceMatrix() const
{
   return fZ;

}

////////////////////////////////////////////////////////////////////////////////
/// Return element i,j from the covariance matrix

Double_t TFumili::GetCovarianceMatrixElement(Int_t i, Int_t j) const
{
   if (!fZ) return 0;
   if (i < 0 || i >= fNpar || j < 0 || j >= fNpar) {
      Error("GetCovarianceMatrixElement","Illegal arguments i=%d, j=%d",i,j);
      return 0;
   }
   return fZ[j+fNpar*i];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the total number of parameters (free + fixed)

Int_t TFumili::GetNumberTotalParameters() const
{
   return fNpar;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of free parameters

Int_t TFumili::GetNumberFreeParameters() const
{
   Int_t nfree = fNpar;
   for (Int_t i=0;i<fNpar;i++) {
      if (IsFixed(i)) nfree--;
   }
   return nfree;
}

////////////////////////////////////////////////////////////////////////////////
/// Return error of parameter ipar

Double_t TFumili::GetParError(Int_t ipar) const
{
   if (ipar<0 || ipar>=fNpar) return 0;
   else return fParamError[ipar];
}

////////////////////////////////////////////////////////////////////////////////
/// Return current value of parameter ipar

Double_t TFumili::GetParameter(Int_t ipar) const
{
   if (ipar<0 || ipar>=fNpar) return 0;
   else return fA[ipar];
}

////////////////////////////////////////////////////////////////////////////////
/// Get various ipar parameter attributes:
///
/// - cname:    parameter name
/// - value:    parameter value
/// - verr:     parameter error
/// - vlow:     lower limit
/// - vhigh:    upper limit
///
/// WARNING! parname must be suitably dimensioned in the calling function.

Int_t TFumili::GetParameter(Int_t ipar,char *cname,Double_t &value,Double_t &verr,Double_t &vlow, Double_t &vhigh) const
{
   if (ipar<0 || ipar>=fNpar) {
      value = 0;
      verr  = 0;
      vlow  = 0;
      vhigh = 0;
      return -1;
   }
   strcpy(cname,fANames[ipar].Data());
   value = fA[ipar];
   verr  = fParamError[ipar];
   vlow  = fAMN[ipar];
   vhigh = fAMX[ipar];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return name of parameter ipar

const char *TFumili::GetParName(Int_t ipar) const
{
   if (ipar < 0 || ipar > fNpar) return "";
   return fANames[ipar].Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Return errors after MINOs
/// not implemented

Int_t TFumili::GetErrors(Int_t ipar,Double_t &eplus, Double_t &eminus, Double_t &eparab, Double_t &globcc) const
{
   eparab = 0;
   globcc = 0;
   if (ipar<0 || ipar>=fNpar) {
      eplus  = 0;
      eminus = 0;
      return -1;
   }
   eplus=fParamError[ipar];
   eminus=-eplus;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return global fit parameters
///  - amin     : chisquare
///  - edm      : estimated distance to minimum
///  - errdef
///  - nvpar    : number of variable parameters
///  - nparx    : total number of parameters

Int_t TFumili::GetStats(Double_t &amin, Double_t &edm, Double_t &errdef, Int_t &nvpar, Int_t &nparx) const
{
   amin   = 2*fS;
   edm    = fGT; //
   errdef = 0; // ??
   nparx  = fNpar;
   nvpar  = 0;
   for(Int_t ii=0; ii<fNpar; ii++) {
      if(fPL0[ii]>0.) nvpar++;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Sum(log(i) i=0,n
/// used by log-likelihood fits

Double_t TFumili::GetSumLog(Int_t n)
{
   if (n < 0) return 0;
   if (n > fNlog) {
      if (fSumLog) delete [] fSumLog;
      fNlog = 2*n+1000;
      fSumLog = new Double_t[fNlog+1];
      Double_t fobs = 0;
      for (Int_t j=0;j<=fNlog;j++) {
         if (j > 1) fobs += TMath::Log(j);
         fSumLog[j] = fobs;
      }
   }
   if (fSumLog) return fSumLog[n];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Inverts packed diagonal matrix Z by square-root method.
///  Matrix elements corresponding to
/// fix parameters are removed.
///
/// - n: number of variable parameters

void TFumili::InvertZ(Int_t n)
{
   static Double_t am = 3.4e138;
   static Double_t rp = 5.0e-14;
   Double_t  ap, aps, c, d;
   Double_t *r_1=fR;
   Double_t *pl_1=fPL;
   Double_t *z_1=fZ;
   Int_t i, k, l, ii, ki, li, kk, ni, ll, nk, nl, ir, lk;
   if (n < 1) {
      return;
   }
   --pl_1;
   --r_1;
   --z_1;
   aps = am / n;
   aps = sqrt(aps);
   ap = 1.0e0 / (aps * aps);
   ir = 0;
   for (i = 1; i <= n; ++i) {
      L1:
         ++ir;
         if (pl_1[ir] <= 0.0e0) goto L1;
         else                   goto L2;
      L2:
         ni = i * (i - 1) / 2;
         ii = ni + i;
         k = n + 1;
         if (z_1[ii] <= rp * TMath::Abs(r_1[ir]) || z_1[ii] <= ap) {
            goto L19;
         }
         z_1[ii] = 1.0e0 / sqrt(z_1[ii]);
         nl = ii - 1;
      L3:
         if (nl - ni <= 0) goto L5;
         else              goto L4;
         L4:
         z_1[nl] *= z_1[ii];
         if (TMath::Abs(z_1[nl]) >= aps) {
            goto L16;
         }
         --nl;
         goto L3;
      L5:
         if (i - n >= 0) goto L12;
         else            goto L6;
      L6:
         --k;
         nk = k * (k - 1) / 2;
         nl = nk;
         kk = nk + i;
         d = z_1[kk] * z_1[ii];
         c = d * z_1[ii];
         l = k;
      L7:
         ll = nk + l;
         li = nl + i;
         z_1[ll] -= z_1[li] * c;
         --l;
         nl -= l;
         if (l - i <= 0) goto L9;
         else            goto L7;
      L8:
         ll = nk + l;
         li = ni + l;
         z_1[ll] -= z_1[li] * d;
      L9:
         --l;
         if (l <= 0) goto L10;
         else        goto L8;
      L10:
         z_1[kk] = -c;
         if (k - i - 1 <= 0) goto L11;
         else                goto L6;
      L11:
         ;
   }
   L12:
      for (i = 1; i <= n; ++i) {
         for (k = i; k <= n; ++k) {
            nl = k * (k - 1) / 2;
            ki = nl + i;
            d = 0.0e0;
            for (l = k; l <= n; ++l) {
               li = nl + i;
               lk = nl + k;
               d += z_1[li] * z_1[lk];
               nl += l;
            }
            ki = k * (k - 1) / 2 + i;
            z_1[ki] = d;
         }
      }
   L15:
      return;
   L16:
      k = i + nl - ii;
      ir = 0;
      for (i = 1; i <= k; ++i) {
         L17:
            ++ir;
            if (pl_1[ir] <= 0.0e0) {
               goto L17;
            }
      }
   L19:
      pl_1[ir] = -2.0e0;
      r_1[ir] = 0.0e0;
      fINDFLG[0] = ir - 1;
      goto L15;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if parameter ipar is fixed, kFALSE otherwise)

Bool_t TFumili::IsFixed(Int_t ipar) const
{
   if(ipar < 0 || ipar >= fNpar) {
      Warning("IsFixed","Illegal parameter number :%d",ipar);
      return kFALSE;
   }
   if (fPL0[ipar] < 0) return kTRUE;
   else                return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Main minimization procedure
///
/// This function is called after setting theoretical function
/// by means of TFumili::SetUserFunc and initializing parameters.
/// Optionally one can set FCN function (see TFumili::SetFCN and TFumili::Eval)
/// If FCN is undefined then user has to provide data arrays by calling
///  TFumili::SetData procedure.
///
/// TFumili::Minimize return following values:
///  -  0  - fit is converged
///  - -2  - function is not decreasing (or bad derivatives)
///  - -3  - error estimations are infinite
///  - -4  - maximum number of iterations is exceeded

Int_t TFumili::Minimize()
{
   Int_t i;
   // Flag3 - is fit is chi2 or likelihood? 0 - chi2, 1 - likelihood
   fINDFLG[2]=0;
   //
   // Are the parameters outside of the boundaries ?
   //
   Int_t parn;

   if(fFCN) {
      Eval(parn,fGr,fS,fA,9); fNfcn++;
   }
   for( i = 0; i < fNpar; i++) {
      if(fA[i] > fAMX[i]) fA[i] = fAMX[i];
      if(fA[i] < fAMN[i]) fA[i] = fAMN[i];
   }

   Int_t nn2, n, fixFLG,  ifix1, fi, nn3, nn1, n0;
   Double_t t1;
   Double_t sp, t, olds=0;
   Double_t bi, aiMAX=0, amb;
   Double_t afix, sigi, akap;
   Double_t alambd, al, bm, abi, abm;
   Int_t l1, k, ifix;

   nn2=0;

   // Number of parameters;
   n=fNpar;
   fixFLG=0;

   // Exit flag
   fENDFLG=0;

   // Flag2
   fINDFLG[1] = 0;
   ifix1=-1;
   fi=0;
   nn3=0;

   // Initialize param.step limits
   for( i=0; i < n; i++) {
      fR[i]=0.;
      if ( fEPS > 0.) fParamError[i] = 0.;
      fPL[i] = fPL0[i];
   }

L3: // Start Iteration

   nn1 = 1;
   t1 = 1.;

L4: // New iteration

   // fS - objective function value - zero first
   fS = 0.;
   // n0 - number of variable parameters in fit
   n0 = 0;
   for( i = 0; i < n; i++) {
      fGr[i]=0.; // zero gradients
      if (fPL0[i] > .0) {
         n0=n0+1;
         // new iteration - new parallelepiped
         if (fPL[i] > .0) fPL0[i]=fPL[i];
      }
   }
   Int_t nn0;
   // Calculate number of fZ-matrix elements as nn0=1+2+..+n0
   nn0 = n0*(n0+1)/2;
   // if (nn0 >= 1) ????
   // fZ-matrix is initialized
   for( i=0; i < nn0; i++) fZ[i]=0.;

   // Flag1
   fINDFLG[0] = 0;
   Int_t ijkl=1;

   // Calculate fS - objective function, fGr - gradients, fZ - fZ-matrix
   if(fFCN) {
      Eval(parn,fGr,fS,fA,2);
      fNfcn++;
   } else
      ijkl = SGZ();
   if(!ijkl) return 10;
   if (ijkl == -1) fINDFLG[0]=1;

   // sp - scaled on fS machine precision
   sp=fRP*TMath::Abs(fS);

   // save fZ-matrix
   for( i=0; i < nn0; i++) fZ0[i] = fZ[i];
   if (nn3 > 0) {
      if (nn1 <= fNstepDec) {
         t=2.*(fS-olds-fGT);
         if (fINDFLG[0] == 0) {
            if (TMath::Abs(fS-olds) <= sp && -fGT <= sp) goto L19;
            if(        0.59*t < -fGT) goto L19;
            t = -fGT/t;
            if (t < 0.25 ) t = 0.25;
         }
         else   t = 0.25;
         fGT = fGT*t;
         t1 = t1*t;
         nn2=0;
         for( i = 0; i < n; i++) {
            if (fPL[i] > 0.) {
               fA[i]=fA[i]-fDA[i];
               fPL[i]=fPL[i]*t;
               fDA[i]=fDA[i]*t;
               fA[i]=fA[i]+fDA[i];
            }
         }
         nn1=nn1+1;
         goto L4;
      }
   }

L19:

   if(fINDFLG[0] != 0) {
      fENDFLG=-4;
      printf("trying to execute an illegal jump at L85\n");
      //goto L85;
   }


   Int_t k1, k2, i1, j, l;
   k1 = 1;
   k2 = 1;
   i1 = 1;
   // In this cycle we removed from fZ contributions from fixed parameters
   // We'll get fixed parameters after boundary check
   for( i = 0; i < n; i++) {
      if (fPL0[i] > .0) {
         // if parameter was fixed - release it
         if (fPL[i] == 0.) fPL[i]=fPL0[i];
         if (fPL[i] > .0) { // ??? it is already non-zero
            // if derivative is negative and we above maximum
            // or vice versa then fix parameter again and increment k1 by i1
            if ((fA[i] >= fAMX[i] && fGr[i] < 0.) ||
                   (fA[i] <= fAMN[i] && fGr[i] > 0.)) {
               fPL[i] = 0.;
               k1 = k1 + i1; // i1 stands for fZ-matrix row-number multiplier
               ///  - skip this row
               //  in case we are fixing parameter number i
            } else {
               for( j=0; j <= i; j++) {// cycle on columns of fZ-matrix
                  if (fPL0[j] > .0) {
                     // if parameter is not fixed then fZ = fZ0
                     // Now matrix fZ of other dimension
                     if (fPL[j] > .0) {
                        fZ[k2 -1] = fZ0[k1 -1];
                        k2=k2+1;
                     }
                     k1=k1+1;
                  }
               }
            }
         }
         else k1 = k1 + i1; // In case of negative fPL[i] - after mconvd
         i1=i1+1;  // Next row of fZ0
      }
   }

   // INVERT fZ-matrix (mconvd() procedure)
   i1 = 1;
   l  = 1;
   for( i = 0; i < n; i++) {// extract diagonal elements to fR-vector
      if (fPL[i] > .0) {
         fR[i] = fZ[l - 1];
         i1 = i1+1;
         l = l + i1;
      }
   }

   n0 = i1 - 1;
   InvertZ(n0);

   // fZ matrix now is inverted
   if (fINDFLG[0] != 0) { // problems
      // some PLs now have negative values, try to reduce fZ-matrix again
      fINDFLG[0] = 0;
      fINDFLG[1] = 1; // errors can be infinite
      fixFLG = fixFLG + 1;
      fi = 0;
      goto L19;
   }

   // ... CALCULATE THEORETICAL STEP TO MINIMUM
   i1 = 1;
   for( i = 0; i < n; i++) {
      fDA[i]=0.; // initial step is zero
      if (fPL[i] > .0) {   // for non-fixed parameters
         l1=1;
         for( l = 0; l < n; l++) {
            if (fPL[l] > .0) {
               // Calculate offset of Z^-1(i1,l1) element in packed matrix
               // because we skip fixed param numbers we need also i,l
               if (i1 <= l1 ) k=l1*(l1-1)/2+i1;
               else k=i1*(i1-1)/2+l1;
               // dA_i = \sum (-Z^{-1}_{il}*grad(fS)_l)
               fDA[i]=fDA[i]-fGr[l]*fZ[k - 1];
               l1=l1+1;
            }
         }
         i1=i1+1;
      }
   }
   //          ... CHECK FOR PARAMETERS ON BOUNDARY

   afix=0.;
   ifix = -1;
   i1 = 1;
   l = i1;
   for( i = 0; i < n; i++)
      if (fPL[i] > .0) {
         sigi = TMath::Sqrt(TMath::Abs(fZ[l - 1])); // calculate \sqrt{Z^{-1}_{ii}}
         fR[i] = fR[i]*fZ[l - 1];      // Z_ii * Z^-1_ii
         if (fEPS > .0) fParamError[i]=sigi;
         if ((fA[i] >= fAMX[i] && fDA[i] > 0.) || (fA[i] <= fAMN[i]
                                               && fDA[i] < .0)) {
            // if parameter out of bounds and if step is making things worse

            akap = TMath::Abs(fDA[i]/sigi);
            // let's found maximum of dA/sigi - the worst of parameter steps
            if (akap > afix) {
               afix=akap;
               ifix=i;
               ifix1=i;
            }
         }
         i1=i1+1;
         l=l+i1;
      }
   if (ifix != -1) {
      // so the worst parameter is found - fix it and exclude,
      //  reduce fZ-matrix again
      fPL[ifix] = -1.;
      fixFLG = fixFLG + 1;
      fi = 0;
      //.. REPEAT CALCULATION OF THEORETICAL STEP AFTER fiXING EACH PARAMETER
      goto L19;
   }

   //... CALCULATE STEP CORRECTION FACTOR

   alambd = 1.;
   fAKAPPA = 0.;
   Int_t imax;
   imax = -1;


   for( i = 0; i < n; i++) {
      if (fPL[i] > .0) {
         bm = fAMX[i] - fA[i];
         abi = fA[i] + fPL[i]; // upper  parameter limit
         abm = fAMX[i];
         if (fDA[i] <= .0) {
            bm = fA[i] - fAMN[i];
            abi = fA[i] - fPL[i]; // lower parameter limit
            abm = fAMN[i];
         }
         bi = fPL[i];
         // if parallelepiped boundary is crossing limits
         // then reduce it (deforming)
         if ( bi > bm) {
            bi = bm;
            abi = abm;
         }
         // if calculated step is out of bounds
         if ( TMath::Abs(fDA[i]) > bi) {
            // decrease step splitter alambdA if needed
            al = TMath::Abs(bi/fDA[i]);
            if (alambd > al) {
               imax=i;
               aiMAX=abi;
               alambd=al;
            }
         }
         // fAKAPPA - parameter will be <fEPS if fit is converged
         akap = TMath::Abs(fDA[i]/fParamError[i]);
         if (akap > fAKAPPA) fAKAPPA=akap;
      }
   }
   //... CALCULATE NEW CORRECTED STEP
   fGT = 0.;
   amb = 1.e18;
   // alambd - multiplier to split theoretical step dA
   if (alambd > .0) amb = 0.25/alambd;
   for( i = 0; i < n; i++) {
      if (fPL[i] > .0) {
         if (nn2 > fNlimMul ) {
            if (TMath::Abs(fDA[i]/fPL[i]) > amb ) {
               fPL[i] = 4.*fPL[i]; // increase parallelepiped
               t1=4.; // flag - that fPL was increased
            }
         }
         // cut step
         fDA[i] = fDA[i]*alambd;
         // expected function value change in next iteration
         fGT = fGT + fDA[i]*fGr[i];
      }
   }

   //.. CHECK IF MINIMUM ATTAINED AND SET EXIT MODE
   // if expected fGT smaller than precision
   // and other stuff
   if (-fGT <= sp && t1 < 1. && alambd < 1.)fENDFLG = -1; // function is not decreasing

   if (fENDFLG >= 0) {
      if (fAKAPPA < TMath::Abs(fEPS)) { // fit is converging
         if (fixFLG == 0)
            fENDFLG=1; // successful fit
         else {// we have fixed parameters
            if (fENDFLG == 0) {
               //... CHECK IF fiXING ON BOUND IS CORRECT
               fENDFLG = 1;
               fixFLG = 0;
               ifix1=-1;
               // release fixed parameters
               for( i = 0; i < fNpar; i++) fPL[i] = fPL0[i];
               fINDFLG[1] = 0;
               // and repeat iteration
               goto L19;
            } else {
               if( ifix1 >= 0) {
                  fi = fi + 1;
                  fENDFLG = 0;
               }
            }
         }
      } else { // fit does not converge
         if( fixFLG != 0) {
            if( fi > fixFLG ) {
               //... CHECK IF fiXING ON BOUND IS CORRECT
               fENDFLG = 1;
               fixFLG = 0;
               ifix1=-1;
               for( i = 0; i < fNpar; i++) fPL[i] = fPL0[i];
               fINDFLG[1] = 0;
               goto L19;
            } else {
               fi = fi + 1;
               fENDFLG = 0;
            }
         } else {
            fi = fi + 1;
            fENDFLG = 0;
         }
      }
   }

// L85:
   // iteration number limit is exceeded
   if(fENDFLG == 0 && nn3 >= fNmaxIter) fENDFLG=-3;

   // fit errors are infinite;
   if(fENDFLG > 0 && fINDFLG[1] > 0) fENDFLG=-2;

   //MONITO (fS,fNpar,nn3,IT,fEPS,fGT,fAKAPPA,alambd);
   if (fENDFLG == 0) { // make step
      for ( i = 0; i < n; i++) fA[i] = fA[i] + fDA[i];
      if (imax >= 0) fA[imax] = aiMAX;
      olds=fS;
      nn2=nn2+1;
      nn3=nn3+1;
   } else {
      // fill covariant matrix VL
      // fill parameter error matrix up
      Int_t il;
      il = 0;
      for( Int_t ip = 0; ip < fNpar; ip++) {
         if( fPL0[ip] > .0) {
            for( Int_t jp = 0; jp <= ip; jp++) {
               if(fPL0[jp] > .0) {
                  //         VL[ind(ip,jp)] = fZ[il];
                  il = il + 1;
               }
            }
         }
      }
      return fENDFLG - 1;
   }
   goto L3;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints fit results.
///
/// ikode is the type of printing parameters
/// p is function value
///
/// - ikode = 1   - print values, errors and limits
/// - ikode = 2   - print values, errors and steps
/// - ikode = 3   - print values, errors, steps and derivatives
/// - ikode = 4   - print only values and errors

void TFumili::PrintResults(Int_t ikode,Double_t p) const
{
   TString exitStatus="";
   TString xsexpl="";
   TString colhdu[3],colhdl[3],cx2,cx3;
   switch (fENDFLG) {
   case 1:
      exitStatus="CONVERGED";
      break;
   case -1:
      exitStatus="CONST FCN";
      xsexpl="****\n* FUNCTION IS NOT DECREASING OR BAD DERIVATIVES\n****";
      break;
   case -2:
      exitStatus="ERRORS INF";
      xsexpl="****\n* ESTIMATED ERRORS ARE INfiNITE\n****";
      break;
   case -3:
      exitStatus="MAX ITER.";
      xsexpl="****\n* MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED\n****";
      break;
   case -4:
      exitStatus="ZERO PROBAB";
      xsexpl="****\n* PROBABILITY OF LIKLIHOOD FUNCTION IS NEGATIVE OR ZERO\n****";
      break;
   default:
      exitStatus="UNDEfiNED";
      xsexpl="****\n* fiT IS IN PROGRESS\n****";
      break;
   }
   if (ikode == 1) {
      colhdu[0] = "              ";
      colhdl[0] = "      ERROR   ";
      colhdu[1] = "      PHYSICAL";
      colhdu[2] = " LIMITS       ";
      colhdl[1] = "    NEGATIVE  ";
      colhdl[2] = "    POSITIVE  ";
   }
   if (ikode == 2) {
      colhdu[0] = "              ";
      colhdl[0] = "      ERROR   ";
      colhdu[1] = "    INTERNAL  ";
      colhdl[1] = "    STEP SIZE ";
      colhdu[2] = "    INTERNAL  ";
      colhdl[2] = "      VALUE   ";
   }
   if (ikode == 3) {
      colhdu[0] = "              ";
      colhdl[0] = "      ERROR   ";
      colhdu[1] = "       STEP   ";
      colhdl[1] = "       SIZE   ";
      colhdu[2] = "       fiRST  ";
      colhdl[2] = "    DERIVATIVE";
   }
   if (ikode == 4) {
      colhdu[0] = "    PARABOLIC ";
      colhdl[0] = "      ERROR   ";
      colhdu[1] = "        MINOS ";
      colhdu[2] = "ERRORS        ";
      colhdl[1] = "   NEGATIVE   ";
      colhdl[2] = "   POSITIVE   ";
   }
   if(fENDFLG<1)Printf("%s", xsexpl.Data());
   Printf(" FCN=%g FROM FUMILI  STATUS=%-10s %9d CALLS OF FCN",
         p,exitStatus.Data(),fNfcn);
   Printf(" EDM=%g ",-fGT);
   Printf("  EXT PARAMETER              %-14s%-14s%-14s",
         (const char*)colhdu[0].Data()
         ,(const char*)colhdu[1].Data()
         ,(const char*)colhdu[2].Data());
   Printf("  NO.   NAME          VALUE  %-14s%-14s%-14s",
         (const char*)colhdl[0].Data()
         ,(const char*)colhdl[1].Data()
         ,(const char*)colhdl[2].Data());

   for (Int_t i=0;i<fNpar;i++) {
      if (ikode==3) {
         cx2 = Form("%14.5e",fDA[i]);
         cx3 = Form("%14.5e",fGr[i]);

      }
      if (ikode==1) {
         cx2 = Form("%14.5e",fAMN[i]);
         cx3 = Form("%14.5e",fAMX[i]);
      }
      if (ikode==2) {
         cx2 = Form("%14.5e",fDA[i]);
         cx3 = Form("%14.5e",fA[i]);
      }
      if(ikode==4){
         cx2 = " *undefined*  ";
         cx3 = " *undefined*  ";
      }
      if(fPL0[i]<=0.) { cx2="    *fixed*   ";cx3=""; }
      Printf("%4d %-11s%14.5e%14.5e%-14s%-14s",i+1
           ,fANames[i].Data(),fA[i],fParamError[i]
           ,cx2.Data(),cx3.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Releases parameter number ipar

void TFumili::ReleaseParameter(Int_t ipar) {
   if(ipar>=0 && ipar<fNpar && fPL0[ipar]<=0.) {
      fPL0[ipar] = -fPL0[ipar];
      if (fPL0[ipar] == 0. || fPL0[ipar]>=1.) fPL0[ipar]=.1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets pointer to data array provided by user.
/// Necessary if SetFCN is not called.
///
/// - numpoints:    number of experimental points
/// - vecsize:      size of data point vector + 2
///                 (for N-dimensional fit vecsize=N+2)
/// - exdata:       data array with following format
///
///   - exdata[0] = ExpValue_0     - experimental data value number 0
///   - exdata[1] = ExpSigma_0     - error of value number 0
///   - exdata[2] = X_0[0]
///   - exdata[3] = X_0[1]
///
///   - exdata[vecsize-1] = X_0[vecsize-3]
///   - exdata[vecsize]   = ExpValue_1
///   - exdata[vecsize+1] = ExpSigma_1
///   - exdata[vecsize+2] = X_1[0]
///
///   - exdata[vecsize*(numpoints-1)] = ExpValue_(numpoints-1)
///
///   - exdata[vecsize*numpoints-1] = X_(numpoints-1)[vecsize-3]

void TFumili::SetData(Double_t *exdata,Int_t numpoints,Int_t vecsize){
   if(exdata){
      fNED1 = numpoints;
      fNED2 = vecsize;
      fEXDA = exdata;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// ret fit method (chisquare or log-likelihood)

void TFumili::SetFitMethod(const char *name)
{
   if (!strcmp(name,"H1FitChisquare"))    SetFCN(H1FitChisquareFumili);
   if (!strcmp(name,"H1FitLikelihood"))   SetFCN(H1FitLikelihoodFumili);
   if (!strcmp(name,"GraphFitChisquare")) SetFCN(GraphFitChisquareFumili);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets for parameter number ipar initial parameter value,
/// name parname, initial error verr and limits vlow and vhigh
/// - If vlow = vhigh but not equil to zero, parameter will be fixed.
/// - If vlow = vhigh = 0, parameter is released and its limits are discarded

Int_t TFumili::SetParameter(Int_t ipar,const char *parname,Double_t value,Double_t verr,Double_t vlow, Double_t vhigh) {
   if (ipar<0 || ipar>=fNpar) return -1;
   fANames[ipar] = parname;
   fA[ipar] = value;
   fParamError[ipar] = verr;
   if(vlow<vhigh) {
      fAMN[ipar] = vlow;
      fAMX[ipar] = vhigh;
   } else {
      if(vhigh<vlow) {
         fAMN[ipar] = vhigh;
         fAMX[ipar] = vlow;
      }
      if(vhigh==vlow) {
         if(vhigh==0.) {
            ReleaseParameter(ipar);
            fAMN[ipar] = gMINDOUBLE;
            fAMX[ipar] = gMAXDOUBLE;
         }
         if(vlow!=0) FixParameter(ipar);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///  Evaluates objective function ( chi-square ), gradients and
///  Z-matrix using data provided by user via TFumili::SetData

Int_t TFumili::SGZ()
{
   fS = 0.;
   Int_t i,j,l,k2=1,k1,ki=0;
   Double_t *x  = new Double_t[fNED2];
   Double_t *df = new Double_t[fNpar];
   Int_t nx = fNED2-2;
   for (l=0;l<fNED1;l++) { // cycle on all exp. points
      k1 = k2;
      if (fLogLike) {
         fNumericDerivatives = kTRUE;
         nx  = fNED2;
         k1 -= 2;
      }

      for (i=0;i<nx;i++){
         ki  += 1+i;
         x[i] = fEXDA[ki];
      }
      //  Double_t y = ARITHM(df,x);
      Double_t y = EvalTFN(df,x);
      if(fNumericDerivatives) Derivatives(df,x);
      Double_t sig=1.;
      if(fLogLike) { // Likelihood method
         if(y>0.) {
            fS = fS - log(y);
            y  = -y;
            sig= y;
         } else { //
            delete [] x;
            delete [] df;
            fS = 1e10;
            return -1; // indflg[0] = 1;
         }
      } else { // Chi2 method
         sig = fEXDA[k2]; // sigma of experimental point
         y = y - fEXDA[k1-1]; // f(x_i) - F_i
         fS = fS + (y*y/(sig*sig))*.5; // simple chi2/2
      }
      Int_t n = 0;
      for (i=0;i<fNpar;i++) {
         if (fPL0[i]>0){
            df[n]   = df[i]/sig; // left only non-fixed param derivatives div by Sig
            fGr[i] += df[n]*(y/sig);
            n++;
         }
      }
      l = 0;
      for (i=0;i<n;i++)
         for (j=0;j<=i;j++)
            fZ[l++] += df[i]*df[j];
      k2 += fNED2;
   }

   delete[] df;
   delete[] x;
   return 1;
}


////////////////////////////////////////////////////////////////////////////////
///  Minimization function for H1s using a Chisquare method.
///  Default method (function evaluated at center of bin)
///  for each point the cache contains the following info
///    - 1D : bc,e,xc  (bin content, error, x of center of bin)
///    - 2D : bc,e,xc,yc
///    - 3D : bc,e,xc,yc,zc

void TFumili::FitChisquare(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
   Foption_t fitOption = GetFitOption();
   if (fitOption.Integral) {
      FitChisquareI(npar,gin,f,u,flag);
      return;
   }
   Double_t cu,eu,fu,fsum;
   Double_t x[3];
   Double_t *zik=0;
   Double_t *pl0=0;

   TH1 *hfit = (TH1*)GetObjectFit();
   TF1 *f1   = (TF1*)GetUserFunc();
   Int_t nd  = hfit->GetDimension();
   Int_t j;

   npar = f1->GetNpar();
   SetParNumber(npar);
   if(flag == 9) return;
   zik = GetZ();
   pl0 = GetPL0();

   Double_t *df = new Double_t[npar];
   f1->InitArgs(x,u);
   f = 0;

   Int_t npfit = 0;
   Double_t *cache = fCache;
   for (Int_t i=0;i<fNpoints;i++) {
      if (nd > 2) x[2]     = cache[4];
      if (nd > 1) x[1]     = cache[3];
      x[0]     = cache[2];
      cu  = cache[0];
      TF1::RejectPoint(kFALSE);
      fu = f1->EvalPar(x,u);
      if (TF1::RejectedPoint()) {cache += fPointSize; continue;}
      eu = cache[1];
      Derivatives(df,x);
      Int_t n = 0;
      fsum = (fu-cu)/eu;
      if (flag!=1) {
         for (j=0;j<npar;j++) {
            if (pl0[j]>0) {
               df[n] = df[j]/eu;
               // left only non-fixed param derivatives / by Sigma
               gin[j] += df[n]*fsum;
               n++;
            }
         }
         Int_t l = 0;
         for (j=0;j<n;j++)
            for (Int_t k=0;k<=j;k++)
               zik[l++] += df[j]*df[k];
      }
      f += .5*fsum*fsum;
      npfit++;
      cache += fPointSize;
   }
   f1->SetNumberFitPoints(npfit);
   delete [] df;
}

////////////////////////////////////////////////////////////////////////////////
///  Minimization function for H1s using a Chisquare method.
///  The "I"ntegral method is used
///  for each point the cache contains the following info
///    - 1D : bc,e,xc,xw  (bin content, error, x of center of bin, x bin width of bin)
///    - 2D : bc,e,xc,xw,yc,yw
///    - 3D : bc,e,xc,xw,yc,yw,zc,zw

void TFumili::FitChisquareI(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
   Double_t cu,eu,fu,fsum;
   Double_t x[3];
   Double_t *zik=0;
   Double_t *pl0=0;

   TH1 *hfit = (TH1*)GetObjectFit();
   TF1 *f1   = (TF1*)GetUserFunc();
   Int_t nd  = hfit->GetDimension();
   Int_t j;

   f1->InitArgs(x,u);
   npar = f1->GetNpar();
   SetParNumber(npar);
   if(flag == 9) return;
   zik = GetZ();
   pl0 = GetPL0();

   Double_t *df=new Double_t[npar];
   f = 0;

   Int_t npfit = 0;
   Double_t *cache = fCache;
   for (Int_t i=0;i<fNpoints;i++) {
      cu  = cache[0];
      TF1::RejectPoint(kFALSE);
      f1->SetParameters(u);
      if (nd < 2) {
         fu = f1->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3])/cache[3];
      } else if (nd < 3) {
         fu = ((TF2*)f1)->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3],cache[4] - 0.5*cache[5],cache[4] + 0.5*cache[5])/(cache[3]*cache[5]);
      } else {
         fu = ((TF3*)f1)->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3],cache[4] - 0.5*cache[5],cache[4] + 0.5*cache[5],cache[6] - 0.5*cache[7],cache[6] + 0.5*cache[7])/(cache[3]*cache[5]*cache[7]);
      }
      if (TF1::RejectedPoint()) {cache += fPointSize; continue;}
      eu = cache[1];
      Derivatives(df,x);
      Int_t n = 0;
      fsum = (fu-cu)/eu;
      if (flag!=1) {
         for (j=0;j<npar;j++) {
            if (pl0[j]>0){
               df[n] = df[j]/eu;
               // left only non-fixed param derivatives / by Sigma
               gin[j] += df[n]*fsum;
               n++;
            }
         }
         Int_t l = 0;
         for (j=0;j<n;j++)
            for (Int_t k=0;k<=j;k++)
               zik[l++] += df[j]*df[k];
      }
      f += .5*fsum*fsum;
      npfit++;
      cache += fPointSize;
   }
   f1->SetNumberFitPoints(npfit);
   delete[] df;
}

////////////////////////////////////////////////////////////////////////////////
///  Minimization function for H1s using a Likelihood method.
///     Basically, it forms the likelihood by determining the Poisson
///     probability that given a number of entries in a particular bin,
///     the fit would predict it's value.  This is then done for each bin,
///     and the sum of the logs is taken as the likelihood.
///
///  Default method (function evaluated at center of bin)
///  for each point the cache contains the following info
///    - 1D : bc,e,xc  (bin content, error, x of center of bin)
///    - 2D : bc,e,xc,yc
///    - 3D : bc,e,xc,yc,zc

void TFumili::FitLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
   Foption_t fitOption = GetFitOption();
   if (fitOption.Integral) {
      FitLikelihoodI(npar,gin,f,u,flag);
      return;
   }
   Double_t cu,fu,fobs,fsub;
   Double_t dersum[100];
   Double_t x[3];
   Int_t icu;

   TH1 *hfit = (TH1*)GetObjectFit();
   TF1 *f1   = (TF1*)GetUserFunc();
   Int_t nd = hfit->GetDimension();
   Int_t j;
   Double_t *zik = GetZ();
   Double_t *pl0 = GetPL0();

   npar = f1->GetNpar();
   SetParNumber(npar);
   if(flag == 9) return;
   Double_t *df=new Double_t[npar];
   if (flag == 2) for (j=0;j<npar;j++) dersum[j] = gin[j] = 0;
   f1->InitArgs(x,u);
   f = 0;

   Int_t npfit = 0;
   Double_t *cache = fCache;
   for (Int_t i=0;i<fNpoints;i++) {
      if (nd > 2) x[2] = cache[4];
      if (nd > 1) x[1] = cache[3];
      x[0]     = cache[2];
      cu  = cache[0];
      TF1::RejectPoint(kFALSE);
      fu = f1->EvalPar(x,u);
      if (TF1::RejectedPoint()) {cache += fPointSize; continue;}
      if (flag == 2) {
         for (j=0;j<npar;j++) {
            dersum[j] += 1; //should be the derivative
            //grad[j] += dersum[j]*(fu-cu)/eu; dersum[j] = 0;
         }
      }
      if (fu < 1.e-9) fu = 1.e-9;
      icu  = Int_t(cu);
      fsub = -fu +icu*TMath::Log(fu);
      fobs = GetSumLog(icu);
      fsub -= fobs;
      Derivatives(df,x);
      int n=0;
      // Here we need gradients of Log likelihood function
      //
      for (j=0;j<npar;j++) {
         if (pl0[j]>0){
            df[n]   = df[j]*(icu/fu-1);
            gin[j] -= df[n];
            n++;
         }
      }
      Int_t l = 0;
      // Z-matrix here - production of first derivatives
      //  of log-likelihood function
      for (j=0;j<n;j++)
         for (Int_t k=0;k<=j;k++)
            zik[l++] += df[j]*df[k];

      f -= fsub;
      npfit++;
      cache += fPointSize;
   }
   f *= 2;
   f1->SetNumberFitPoints(npfit);
   delete[] df;
}

////////////////////////////////////////////////////////////////////////////////
///  Minimization function for H1s using a Likelihood method.
///     Basically, it forms the likelihood by determining the Poisson
///     probability that given a number of entries in a particular bin,
///     the fit would predict it's value.  This is then done for each bin,
///     and the sum of the logs is taken as the likelihood.
///
///  The "I"ntegral method is used
///  for each point the cache contains the following info
///    - 1D : bc,e,xc,xw  (bin content, error, x of center of bin, x bin width of bin)
///    - 2D : bc,e,xc,xw,yc,yw
///    - 3D : bc,e,xc,xw,yc,yw,zc,zw

void TFumili::FitLikelihoodI(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
   Double_t cu,fu,fobs,fsub;
   Double_t dersum[100];
   Double_t x[3];
   Int_t icu;

   TH1 *hfit = (TH1*)GetObjectFit();
   TF1 *f1   = (TF1*)GetUserFunc();
   Int_t nd = hfit->GetDimension();
   Int_t j;
   Double_t *zik = GetZ();
   Double_t *pl0 = GetPL0();

   Double_t *df=new Double_t[npar];

   npar = f1->GetNpar();
   SetParNumber(npar);
   if(flag == 9) {delete [] df; return;}
   if (flag == 2) for (j=0;j<npar;j++) dersum[j] = gin[j] = 0;
   f1->InitArgs(x,u);
   f = 0;

   Int_t npfit = 0;
   Double_t *cache = fCache;
   for (Int_t i=0;i<fNpoints;i++) {
      if (nd > 2) x[2] = cache[4];
      if (nd > 1) x[1] = cache[3];
      x[0]     = cache[2];
      cu  = cache[0];
      TF1::RejectPoint(kFALSE);
      if (nd < 2) {
         fu = f1->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3])/cache[3];
      } else if (nd < 3) {
         fu = ((TF2*)f1)->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3],cache[4] - 0.5*cache[5],cache[4] + 0.5*cache[5])/(cache[3]*cache[5]);
      } else {
         fu = ((TF3*)f1)->Integral(cache[2] - 0.5*cache[3],cache[2] + 0.5*cache[3],cache[4] - 0.5*cache[5],cache[4] + 0.5*cache[5],cache[6] - 0.5*cache[7],cache[6] + 0.5*cache[7])/(cache[3]*cache[5]*cache[7]);
      }
      if (TF1::RejectedPoint()) {cache += fPointSize; continue;}
      if (flag == 2) {
         for (j=0;j<npar;j++) {
            dersum[j] += 1; //should be the derivative
            //grad[j] += dersum[j]*(fu-cu)/eu; dersum[j] = 0;
         }
      }
      if (fu < 1.e-9) fu = 1.e-9;
      icu  = Int_t(cu);
      fsub = -fu +icu*TMath::Log(fu);
      fobs = GetSumLog(icu);
      fsub -= fobs;
      Derivatives(df,x);
      int n=0;
      // Here we need gradients of Log likelihood function
      //
      for (j=0;j<npar;j++) {
         if (pl0[j]>0){
            df[n]   = df[j]*(icu/fu-1);
            gin[j] -= df[n];
            n++;
         }
      }
      Int_t l = 0;
      // Z-matrix here - production of first derivatives
      //  of log-likelihood function
      for (j=0;j<n;j++)
         for (Int_t k=0;k<=j;k++)
            zik[l++] += df[j]*df[k];

      f -= fsub;
      npfit++;
      cache += fPointSize;
   }
   f *= 2;
   f1->SetNumberFitPoints(npfit);
   delete[] df;
}


//______________________________________________________________________________
//
//  STATIC functions
//______________________________________________________________________________

////////////////////////////////////////////////////////////////////////////////
/// Minimization function for H1s using a Chisquare method.

void H1FitChisquareFumili(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
   TFumili *hFitter = (TFumili*)TVirtualFitter::GetFitter();
   hFitter->FitChisquare(npar, gin, f, u, flag);
}

////////////////////////////////////////////////////////////////////////////////
/// Minimization function for H1s using a Likelihood method.
///     Basically, it forms the likelihood by determining the Poisson
///     probability that given a number of entries in a particular bin,
///     the fit would predict it's value.  This is then done for each bin,
///     and the sum of the logs is taken as the likelihood.
///     PDF:  P=exp(-f(x_i))/[F_i]!*(f(x_i))^[F_i]
///    where F_i - experimental value, f(x_i) - expected theoretical value
///    [F_i] - integer part of F_i.
///    drawback is that if F_i>Int_t - GetSumLog will fail
///    for big F_i is faster to use Euler's Gamma-function

void H1FitLikelihoodFumili(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{

   TFumili *hFitter = (TFumili*)TVirtualFitter::GetFitter();
   hFitter->FitLikelihood(npar, gin, f, u, flag);
}

////////////////////////////////////////////////////////////////////////////////
/// Minimization function for Graphs using a Chisquare method.
/// In case of a TGraphErrors object, ex, the error along x,  is projected
/// along the y-direction by calculating the function at the points x-exlow and
/// x+exhigh.
///
/// The chisquare is computed as the sum of the quantity below at each point:
///
///                     (y - f(x))**2
///         -----------------------------------
///         ey**2 + (0.5*(exl + exh)*f'(x))**2
///
/// where x and y are the point coordinates and f'(x) is the derivative of function f(x).
/// This method to approximate the uncertainty in y because of the errors in x, is called
/// "effective variance" method.
/// The improvement, compared to the previously used  method (f(x+ exhigh) - f(x-exlow))/2
/// is of (error of x)**2 order.
///
///  NOTE:
///
///  1. By using the "effective variance" method a simple linear regression
///      becomes a non-linear case , which takes several iterations
///      instead of 0 as in the linear case .
///
///  2. The effective variance technique assumes that there is no correlation
///      between the x and y coordinate .
///
/// In case the function lies below (above) the data point, ey is ey_low (ey_high).

void GraphFitChisquareFumili(Int_t &npar, Double_t * gin, Double_t &f,
                       Double_t *u, Int_t flag)
{
   Double_t cu,eu,exl,exh,ey,eux,fu,fsum;
   Double_t x[1];
   Int_t i, bin, npfits=0;

   TFumili *grFitter = (TFumili*)TVirtualFitter::GetFitter();
   TGraph *gr     = (TGraph*)grFitter->GetObjectFit();
   TF1 *f1   = (TF1*)grFitter->GetUserFunc();
   Foption_t fitOption = grFitter->GetFitOption();

   Int_t n        = gr->GetN();
   Double_t *gx   = gr->GetX();
   Double_t *gy   = gr->GetY();
   npar           = f1->GetNpar();

   grFitter->SetParNumber(npar);

   if(flag == 9) return;
   Double_t *zik = grFitter->GetZ();
   Double_t *pl0 = grFitter->GetPL0();
   Double_t *df  = new Double_t[npar];


   f1->InitArgs(x,u);
   f      = 0;
   for (bin=0;bin<n;bin++) {
      x[0] = gx[bin];
      if (!f1->IsInside(x)) continue;
      cu   = gy[bin];
      TF1::RejectPoint(kFALSE);
      fu   = f1->EvalPar(x,u);
      if (TF1::RejectedPoint()) continue;
      npfits++;
      Double_t eusq=1.;
      if (fitOption.W1) {
         eu = 1.;
      } else {
         exh  = gr->GetErrorXhigh(bin);
         exl  = gr->GetErrorXlow(bin);
         ey  = gr->GetErrorY(bin);
         if (exl < 0) exl = 0;
         if (exh < 0) exh = 0;
         if (ey < 0)  ey  = 0;
         if (exh > 0 && exl > 0) {
//          "Effective variance" method for projecting errors
            eux = 0.5*(exl + exh)*f1->Derivative(x[0], u);
         } else
            eux = 0.;
         eu = ey*ey+eux*eux;
         if (eu <= 0) eu = 1;
         eusq = TMath::Sqrt(eu);
      }
      grFitter->Derivatives(df,x);
      n = 0;
      fsum = (fu-cu)/eusq;
      for (i=0;i<npar;i++) {
         if (pl0[i]>0){
            df[n] = df[i]/eusq;
            // left only non-fixed param derivatives / by Sigma
            gin[i] += df[n]*fsum;
            n++;
         }
      }
      Int_t l = 0;
      for (i=0;i<n;i++)
         for (Int_t j=0;j<=i;j++)
            zik[l++] += df[i]*df[j];
      f += .5*fsum*fsum;

   }
   delete [] df;
   f1->SetNumberFitPoints(npfits);
}

