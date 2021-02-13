// @(#)root/hist:$Id$
// Authors: Lorenzo Moneta, Aurélie Flandi  27/08/14
//
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  ROOT  Team, CERN/PH-SFT                        *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "TROOT.h"
#include "TClass.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1NormSum.h"
#include "Math/WrappedFunction.h"
#include "Math/WrappedTF1.h"

ClassImp(TF1NormSum);

/** \class TF1NormSum
    \ingroup Hist
Class adding two functions: c1*f1+c2*f2
*/

////////////////////////////////////////////////////////////////////////////////
/// Function to find and rename duplicate parameters with the same name

template<class Iterator>
void FixDuplicateNames(Iterator begin, Iterator end) {

   // make a map of values

   std::multimap<TString, int > parMap;
   for (Iterator it = begin; it !=  end; ++it) {
      parMap.insert( std::make_pair( *it, std::distance(begin,it) ) );
   }
   for ( auto & elem : parMap) {
      TString name = elem.first;
      int n = parMap.count( name);
      if (n > 1 ) {
          std::pair <std::multimap<TString,int>::iterator, std::multimap<TString,int>::iterator> ret;
          ret = parMap.equal_range(name);
          int i = 0;
          for (std::multimap<TString,int>::iterator it=ret.first; it!=ret.second; ++it) {
             *(begin+it->second) = TString::Format("%s%d",name.Data(),++i);
          }
      }
   }

}

////////////////////////////////////////////////////////////////////////////////

void TF1NormSum::InitializeDataMembers(const std::vector<TF1 *> &functions, const std::vector<Double_t> &coeffs,
                                       Double_t scale)
{

   fScale           = scale;
   fCoeffs          = coeffs;
   fNOfFunctions    = functions.size();
   fCstIndexes      = std::vector < Int_t     > (fNOfFunctions);
   fParNames        = std::vector<TString> (fNOfFunctions);
   fParNames.reserve(3*fNOfFunctions);  // enlarge capacity for function parameters

   // fill fFunctions with unique_ptr's
   fFunctions = std::vector<std::unique_ptr<TF1>>(functions.size());
   for (unsigned int n = 0; n < fNOfFunctions; n++) {
      // use TF1::Copy and not clone to copy the TF1 pointers
      // and use IsA()::New() in case we have base class pointers 
      TF1 * f = (TF1*) functions[n]->IsA()->New();
      functions[n]->Copy(*f);
      fFunctions[n] = std::unique_ptr<TF1>(f);
     

      if (!fFunctions[n])
         Fatal("InitializeDataMembers", "Invalid input function -- abort");

      fFunctions[n]->SetBit(TF1::kNotGlobal, kTRUE);
   }

   for (unsigned int n=0; n < fNOfFunctions; n++)
   {
      int npar = fFunctions[n] -> GetNpar();
      fCstIndexes[n]      = fFunctions[n] -> GetParNumber("Constant");//return -1 if there is no constant parameter
      fParNames[n] = TString::Format("Coeff%d",n);
      if (fCstIndexes[n]!= -1)                                        //if there exists a constant parameter
      {
         fFunctions[n] -> FixParameter(fCstIndexes[n], 1.); // fixes the parameters called "Constant" to 1
         int k = 0;                                         // index for the temp array, k wil go form 0 until fNofNonCstParameter
         for (int i=0; i<npar; i++)                         // go through all the parameter to
         {
            if (i==fCstIndexes[n])   continue;              // go to next step if this is the constant parameter
            fParNames.push_back(  fFunctions[n] -> GetParName(i) );
            k++;
         }
      }
      else {
         for (int i=0; i < npar; i++)                        //go through all the parameter to
         {
            fParNames.push_back( fFunctions[n] -> GetParName(i) );
         }
      }
      //normalize the functions if it is not already done (do at the end so constant parameter is not zero)
      if (!fFunctions[n] -> IsEvalNormalized())  fFunctions[n]  -> SetNormalized(true);
   }

   // Set range
   if (fNOfFunctions == 0) {
      fXmin = 0.;
      fXmax = 1.;
      // Info("InitializeDataMembers", "Initializing empty TF1NormSum with default [0,1] range");
   } else {
      fFunctions[0]->GetRange(fXmin, fXmax);
      if (fXmin >= fXmax) {
         fXmin = 0.;
         fXmax = 1.;
         // Info("InitializeDataMembers", "Initializing empty TF1NormSum with default [0,1] range");
      }
      for (unsigned int n = 1; n < fNOfFunctions; n++) {
         fFunctions[n]->SetRange(fXmin, fXmax);
         fFunctions[n]->Update();
      }
   }

   FixDuplicateNames(fParNames.begin() + fNOfFunctions, fParNames.end());
}

////////////////////////////////////////////////////////////////////////////////

TF1NormSum::TF1NormSum()
{
   fNOfFunctions  = 0;
   fScale         = 1.;
   fFunctions = std::vector<std::unique_ptr<TF1>>(0);         // Vector of size fNOfFunctions containing TF1 functions
   fCoeffs        = std::vector < Double_t  >(0) ;            // Vector of size fNOfFunctions containing coefficients in front of each function
   fCstIndexes = std::vector < Int_t     > (0);
   fXmin = 0; // Dummy values of xmin and xmax
   fXmax = 1;
}

////////////////////////////////////////////////////////////////////////////////

TF1NormSum::TF1NormSum(const std::vector <TF1*> &functions, const std::vector <Double_t> &coeffs, Double_t scale)
{
   InitializeDataMembers(functions, coeffs, scale);
}

////////////////////////////////////////////////////////////////////////////////
/// TF1NormSum constructor taking 2 functions, and 2 coefficients (if not equal to 1)

TF1NormSum::TF1NormSum(TF1* function1, TF1* function2, Double_t coeff1, Double_t coeff2, Double_t scale)
{
   std::vector<TF1 *> functions(2);
   std::vector < Double_t > coeffs(2);

   functions = {function1, function2};
   coeffs = {coeff1, coeff2};

   InitializeDataMembers(functions, coeffs,scale);
}

////////////////////////////////////////////////////////////////////////////////
/// TF1NormSum constructor taking 3 functions, and 3 coefficients (if not equal to 1)

TF1NormSum::TF1NormSum(TF1* function1, TF1* function2, TF1* function3, Double_t coeff1, Double_t coeff2, Double_t coeff3, Double_t scale)
{
   std::vector<TF1 *> functions(3);
   std::vector < Double_t > coeffs(3);

   functions = {function1, function2, function3};
   coeffs = {coeff1, coeff2, coeff3};

   InitializeDataMembers(functions, coeffs,scale);
}

////////////////////////////////////////////////////////////////////////////////
/// TF1NormSum constructor taking any addition of formulas with coefficient or not
///
/// - example 1 : 2.*expo + gauss + 0.5* gauss
/// - example 2 : expo + 0.3*f1 if f1 is defined in the list of functions

TF1NormSum::TF1NormSum(const TString &formula, Double_t xmin, Double_t xmax)
{
   TF1::InitStandardFunctions();

   TObjArray *arrayall    = formula.Tokenize("*+");
   TObjArray *arraytimes  = formula.Tokenize("*") ;
   Int_t noffunctions     = (formula.Tokenize("+")) -> GetEntries();
   Int_t nofobj           = arrayall  -> GetEntries();
   Int_t nofcoeffs        = nofobj - noffunctions;

   std::vector<TF1 *> functions(noffunctions);
   std::vector < Double_t > coeffs(noffunctions);
   std::vector < TString  > funcstringall(nofobj);
   std::vector < Int_t    > indexsizetimes(nofcoeffs+1);
   std::vector < Bool_t   > isacoeff(nofobj);//1 is it is a coeff, 0 if it is a functions

   for (int i=0; i<nofobj; i++)
   {
      funcstringall[i] = ((TObjString*)((*arrayall)[i])) -> GetString();
      funcstringall[i].ReplaceAll(" ","");
   }
   //algorithm to determine which object is a coefficient and which is a function
   //uses the fact that the last item of funcstringtimes[i].Tokenize("+") is always a coeff.
   Int_t j = 0;
   Int_t k = 1;
   for (int i=0; i<nofcoeffs+1; i++)
   {
      indexsizetimes[i] = ( ( ( (TObjString*)(*arraytimes)[i] ) -> GetString() ).Tokenize("+") ) -> GetEntries();
      while (k < indexsizetimes[i])
      {
         isacoeff[k+j-1] = 0;
         k++;
      }
      j = j+indexsizetimes[i];
      if (j==nofobj)    isacoeff[j-1] = 0;    //the last one is never a coeff
      else              isacoeff[j-1] = 1;
      k = 1;
   }

   Double_t old_xmin = 0.0, old_xmax = 0.0;
   k = 0; // index of term in funcstringall
   for (int i=0; i<noffunctions; i++)
   {
      // first, handle coefficient
      if (isacoeff[k]) {
         coeffs[i] = funcstringall[k].Atof();
         k++;
      } else {
         coeffs[i] = 1.;
      }

      // then, handle function
      functions[i] = (TF1 *)(gROOT->GetListOfFunctions()->FindObject(funcstringall[k]));
      if (!functions[i])
         Error("TF1NormSum", "Function %s does not exist", funcstringall[k].Data());
      // (set range for first function, which determines range of whole TF1NormSum)
      if (i == 0) {
         functions[i]->GetRange(old_xmin, old_xmax);
         functions[i]->SetRange(xmin, xmax);
      }

      k++;
   }
   InitializeDataMembers(functions, coeffs,1.);

   // Set range of first function back to original state
   if (noffunctions > 0 && functions[0])
      functions[0]->SetRange(old_xmin, old_xmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor (necessary to hold unique_ptr as member variable)

TF1NormSum::TF1NormSum(const TF1NormSum &nsum)
{
   nsum.Copy((TObject &)*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TF1NormSum &TF1NormSum::operator=(const TF1NormSum &rhs)
{
   if (this != &rhs)
      rhs.Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Overload the parenthesis to add the functions

double TF1NormSum::operator()(const Double_t *x, const Double_t *p)
{
   // first refresh the parameters
   if (p != 0)
      SetParameters(p);

   Double_t sum = 0.;
   for (unsigned int n=0; n<fNOfFunctions; n++)
      sum += fCoeffs[n]*(fFunctions[n] -> EvalPar(x,0));

   // normalize by a scale parameter (typically the bin width)
   return fScale * sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return array of parameters

std::vector<double>  TF1NormSum::GetParameters() const {
   std::vector<double> params(GetNpar() );
   int offset = 0;
   int nOfNonCstParams = 0;
   for (unsigned int n=0; n<fNOfFunctions; n++)
   {
      params[n] = fCoeffs[n];   // copy the coefficients
      offset += nOfNonCstParams;           // offset to go along the list of parameters
      int k = 0;
      for (int j = 0; j < fFunctions[n]->GetNpar(); ++j) {
         if (j != fCstIndexes[n]) {
            params[k+fNOfFunctions+offset] = fFunctions[n]->GetParameter(j);
            k++;
         }
      }
      nOfNonCstParams = k;
   }
   return params;
}
////////////////////////////////////////////////////////////////////////////////
/// Initialize array of all parameters.
///
/// double *params must contains first an array of the coefficients, then an array of the parameters.

void TF1NormSum::SetParameters(const Double_t *params) // params should have the size [fNOfFunctions][fNOfNonCstParams]
{
   for (unsigned int n=0; n<fNOfFunctions; n++)                         //initialization of the coefficients
   {
      fCoeffs[n] = params[n];
   }
   Int_t    offset     = 0;
   int k = 0;  // k indicates the number of non-constant parameter per function
   for (unsigned int n=0; n<fNOfFunctions; n++)
   {
      bool equalParams = true;
      Double_t * funcParams = fFunctions[n]->GetParameters();
      int npar = fFunctions[n]->GetNpar();
      offset += k;      // offset to go along the list of parameters
      k = 0; // reset k value for next function
      for (int i = 0; i < npar; ++i) {
         // constant parameters can be only one
         if (i != fCstIndexes[n])
         {
            // check if they are equal
            equalParams &= (funcParams[i] == params[k+fNOfFunctions+offset] );
            funcParams[i] = params[k+fNOfFunctions+offset];
            k++;
         }
      }
      // update function integral if not equal
      if (!equalParams) fFunctions[n]->Update();

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize array of all parameters.
///
/// Overload the TF1::SetParameters() method.
/// A maximum of 10 parameters must be used, with first the coefficients, then the parameters

void TF1NormSum::SetParameters(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4,
                               Double_t p5, Double_t p6, Double_t p7, Double_t p8, Double_t p9, Double_t p10)
{
   const double params[] = {p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10};
   TF1NormSum::SetParameters(params);

}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of (non constant) parameters including the coefficients: for 2 functions: c1,c2,p0,p1,p2,p3...

Int_t TF1NormSum::GetNpar() const
{
   Int_t nofparams = 0;
   for (unsigned int n=0; n<fNOfFunctions; ++n)
   {
      nofparams += fFunctions[n]->GetNpar();
      if (fCstIndexes[n] >= 0) nofparams -= 1;
   }
   return nofparams + fNOfFunctions;  //fNOfFunctions for the coefficients
}

////////////////////////////////////////////////////////////////////////////////

void TF1NormSum::SetRange(Double_t a, Double_t b)
{
   if (a >= b) {
      Warning("SetRange", "Invalid range: %f >= %f", a, b);
      return;
   }

   fXmin = a;
   fXmax = b;

   for (unsigned int n = 0; n < fNOfFunctions; n++) {
      fFunctions[n]->SetRange(a, b);
      fFunctions[n]->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////

void TF1NormSum::GetRange(Double_t &a, Double_t &b) const
{
   a = fXmin;
   b = fXmax;
}

////////////////////////////////////////////////////////////////////////////////
///   Update the component functions of the normalized sum

void TF1NormSum::Update()
{
   for (unsigned int n = 0; n < fNOfFunctions; n++)
      fFunctions[n]->Update();
}

////////////////////////////////////////////////////////////////////////////////

void TF1NormSum::Copy(TObject &obj) const
{
   ((TF1NormSum &)obj).fNOfFunctions = fNOfFunctions;
   ((TF1NormSum &)obj).fScale = fScale;
   ((TF1NormSum &)obj).fXmin = fXmin;
   ((TF1NormSum &)obj).fXmax = fXmax;
   ((TF1NormSum &)obj).fCoeffs = fCoeffs;
   ((TF1NormSum &)obj).fCstIndexes = fCstIndexes;
   ((TF1NormSum &)obj).fParNames = fParNames;

   // Clone objects in unique_ptr's
   ((TF1NormSum &)obj).fFunctions = std::vector<std::unique_ptr<TF1>>(fNOfFunctions);
   for (unsigned int n = 0; n < fNOfFunctions; n++) {
      TF1 * f = (TF1*) fFunctions[n]->IsA()->New();   
      fFunctions[n]->Copy(*f);
      ((TF1NormSum &)obj).fFunctions[n] = std::unique_ptr<TF1>(f);
   }
}
