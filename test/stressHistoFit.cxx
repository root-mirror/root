// @(#)root/test:$name:  $:$id: stressHistoFit.cxx,v 1.15 2002/10/25 10:47:51 rdm exp $
// Authors: David Gonzalez Maline November 2008

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*//
//                                                                               //
//                                                                               //
// Set of tests for different minimization algorithms and for                    //
// different objects. The tests are divided into three types:                    //
//                                                                               //
// 1. 1D and 2D Objects, including 1D and 2D histograms, 1D and 2D               //
//    histograms with variable bins, TGraph, TGraphErrors, TGraph2D,             //
//    TGraph2DErrors                                                             //
// 2. Same as before, but trying linear fitters.                                 //
// 3. Unbinned fits with trees of different dimensions.                          //
//                                                                               //
// Each test will performed fits with different functions and                    //
// different minimization algorithms selected. There is an error                 //
// tolerance for each one of them. There is also the possibility to              //
// inspect each one of the test individually changing the                        //
// defaultOptions variable.                                                      //
//                                                                               //
//                                                                               //
// An example of output when all the tests run OK is shown below:                //
// ****************************************************************************
// *  Starting  stress  H I S T O F I T                                       *
// ****************************************************************************

// Test 1D and 2D objects

// Test   1:  'Histogram 1D Variable' with 'GAUS'...................OK
// Test   2:  'Histogram 1D' with 'GAUS'............................OK
// Test   3:  'TGraph 1D' with 'GAUS'...............................OK
// Test   4:  'TGraphErrors 1D' with 'GAUS'.........................OK
// Test   5:  'THnSparse 1D' with 'GAUS'............................OK
// Test   6:  'Histogram 1D Variable' with 'Polynomial'.............OK
// Test   7:  'Histogram 1D' with 'Polynomial'......................OK
// Test   8:  'TGraph 1D' with 'Polynomial'.........................OK
// Test   9:  'TGraphErrors 1D' with 'Polynomial'...................OK
// Test  10:  'THnSparse 1D' with 'Polynomial'......................OK
// Test  11:  'Histogram 2D Variable' with 'gaus2D'.................OK
// Test  12:  'Histogram 2D' with 'gaus2D'..........................OK
// Test  13:  'TGraph 2D' with 'gaus2D'.............................OK
// Test  14:  'TGraphErrors 2DGE' with 'gaus2D'.....................OK
// Test  15:  'THnSparse 2D' with 'gaus2D'..........................OK

// Test Linear fits

// Test  16:  'Histogram 1D Variable' with 'Polynomial'.............OK
// Test  17:  'Histogram 1D' with 'Polynomial'......................OK
// Test  18:  'TGraph 1D' with 'Polynomial'.........................OK
// Test  19:  'TGraphErrors 1D' with 'Polynomial'...................OK
// Test  20:  'THnSparse 1D' with 'Polynomial'......................OK
// Test  21:  'Histogram 2D Variable' with 'Poly2D'.................OK
// Test  22:  'Histogram 2D' with 'Poly2D'..........................OK
// Test  23:  'TGraph 2D' with 'Poly2D'.............................OK
// Test  24:  'TGraphErrors 2DGE' with 'Poly2D'.....................OK
// Test  25:  'THnSparse 2D' with 'Poly2D'..........................OK

// Test unbinned fits

// Test  26:  'tree' with 'gausn'...................................OK
// Test  27:  'tree' with 'gaus2Dn'.................................OK
// Test  28:  'tree' with 'gausND'..................................OK

// ****************************************************************************
// stressHistoFit: Real Time =  37.49 seconds Cpu Time =  37.24 seconds
//  ROOTMARKS = 2663.8 ROOT version: 5.27/01 trunk@32822
// ****************************************************************************
//
//                                                                               //
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*//

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
#include "snprintf.h"

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
#include "TMath.h"
#include "Math/DistFunc.h"

#include "TUnuran.h"
#include "TUnuranMultiContDist.h"
#include "Math/MinimizerOptions.h"
#include "Math/IntegratorOptions.h"
#include "TBackCompFitter.h"
#include "TVirtualFitter.h"

#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"

#include "TRandom3.h"

#include "TROOT.h"
//#include "RConfigure.h"
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TApplication.h"

#include <vector>
#include <string>
#include <cassert>
#include <cmath>

#ifdef R__WIN32
#ifndef __CLING__
#define FOREGROUND_BLUE      1
#define FOREGROUND_GREEN     2
#define FOREGROUND_RED       4
#define FOREGROUND_INTENSITY 8
#define FOREGROUND_GREY FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE
extern "C" {
   void *__stdcall GetStdHandle(unsigned long);
   bool __stdcall SetConsoleTextAttribute(void *, unsigned int);
}
#pragma comment(lib, "Kernel32.lib")
#endif
#endif

#include "Riostream.h"
using namespace std;
// Next line should not exist. It is now there for testing
// pourpuses.
#undef R__HAS_MATHMORE

unsigned int __DRAW__ = 0;

int gSelectedTest = 0;

bool gEnableMT = kFALSE;

// set a small tolerance for the tests
// The default of 10*-2 make sometimes Simplex do not converge
const double gDefaultTolerance = 1.E-4;

TRandom3 rndm;

enum cmpOpts {
   cmpNone = 0,
   cmpPars = 1,
   cmpDist = 2,
   cmpChi2 = 4,
   cmpErr  = 8,
};

// Reference structure to compare the fitting results
struct RefValue {
   const double* pars;
   const double  chi;
   RefValue(const double* p = 0, const double c = 0.0): pars(p), chi(c) {};
};

// Class that keeps a reference structure and some tolerance values to
// make a comparision between the reference and the result of a
// fit. The options define what has to be compared.
class CompareResult {
public:
   struct RefValue* refValue;
   int opts;
   double tolPar;
   double tolChi2;
   CompareResult(int _opts = cmpPars, double _tolPar = 3, double _tolChi2 = 0.01):
      refValue(0), opts(_opts), tolPar(_tolPar), tolChi2(_tolChi2) {};

   // use default copy-ctor and assignment operator
   
   void setRefValue(struct RefValue* _refValue)
   {
      refValue = _refValue;
   };

   int parameters(int npar, double val, double ref) const
   {
      int ret = 0;
      if ( refValue && (opts & cmpPars) )
      {
         ret = compareResult(val, refValue->pars[npar], tolPar*ref);
//          printf("[TOL:%f]", ref);
      }
      return ret;
   };

   int chi2(double val) const
   { return ( refValue && (opts & cmpChi2) ) ? compareResult(val, refValue->chi, tolChi2) : 0; };

public:
   // Compares two doubles with a given tolerence
   int compareResult(double v1, double v2, double tol = 0.01) const {
      if (std::abs(v1-v2) > tol ) return 1;
      return 0;
   }
};

// Create a variable range in a vector (to be passed to the histogram
// constructor
void FillVariableRange(Double_t v[], Int_t numberOfBins, Double_t minRange, Double_t maxRange)
{
   Double_t minLimit = (maxRange-minRange)  / (numberOfBins*4);
   Double_t maxLimit = (maxRange-minRange)*4/ (numberOfBins);
   v[0] = 0;
   for ( Int_t i = 1; i < numberOfBins + 1; ++i )
   {
      Double_t limit = rndm.Uniform(minLimit, maxLimit);
      v[i] = v[i-1] + limit;
   }

   Double_t k = (maxRange-minRange)/v[numberOfBins];
   for ( Int_t i = 0; i < numberOfBins + 1; ++i )
   {
      v[i] = v[i] * k + minRange;
   }
}

// Class defining the different algorithms. It contains the library,
// the particular algorithm and the options which will be used to
// invoke the algorithm. It also contains a CompareResult to indicate
// what sort of checking has to be done once the algorithm has been
// used.
class algoType {
public:
   const char* type;
   const char* algo;
   const char* opts;
   CompareResult cmpResult;

   algoType(): type(0), algo(0), opts(0), cmpResult(0) {}

   algoType(const char* s1, const char* s2, const char* s3,
            CompareResult _cmpResult):
      type(s1), algo(s2), opts(s3), cmpResult(_cmpResult) {}
};

// Different vectors containing the list of algorithms to be used.
vector<algoType> commonAlgos;
vector<algoType> simplexAlgos;
vector<algoType> specialAlgos;
vector<algoType> noGraphAlgos;
vector<algoType> noGraphErrorAlgos;
vector<algoType> graphErrorAlgos;
vector<algoType> histGaus2D;
vector<algoType> linearAlgos;

vector< vector<algoType> > listTH1DAlgos;
vector< vector<algoType> > listAlgosTGraph;
vector< vector<algoType> > listAlgosTGraphError;

vector< vector<algoType> > listLinearAlgos;

vector< vector<algoType> > listTH2DAlgos;
vector< vector<algoType> > listAlgosTGraph2D;
vector< vector<algoType> > listAlgosTGraph2DError;


// Class defining the limits in the parameters of a function.
class ParLimit {
public:
   int npar;
   double min;
   double max;
   ParLimit(int _npar = 0, double _min = 0, double _max = 0): npar(_npar), min(_min), max(_max) {};
};

// Set the limits of a function given a vector of ParLimit
void SetParsLimits(vector<ParLimit>& v, TF1* func)
{
   for ( vector<ParLimit>::iterator it = v.begin();
         it !=  v.end(); ++it ) {
//       printf("Setting parameters: %d, %f, %f\n", (*it)->npar, (*it)->min, (*it)->max);
      func->SetParLimits( it->npar, it->min, it->max);
   }
}

// Class that defined a fitting function. It will contain:
//     The name of the function
//     A pointer to the method that implements the function
//     origPars is the original parameters used to fill the histogram/object
//     fitPars parameters used right before fitting.
//     parLimits limits of the parameters to be set before fitting
class fitFunctions {
public:
   TString name;
   double (*func)(double*, double*) = nullptr;
   unsigned int npars = 0;
   vector<double> origPars;
   vector<double> fitPars;
   vector<ParLimit> parLimits;

   fitFunctions() {}

   fitFunctions(const char* s1, double (*f)(double*, double*),
                unsigned int n,
                double* v1, double* v2,
                vector<ParLimit>& limits):
      name(s1), func(f), npars(n),
      origPars(npars), fitPars(npars), parLimits(limits.size())
   {
      copy(v1, v1 + npars, origPars.begin());
      copy(v2, v2 + npars, fitPars.begin());
      copy(limits.begin(), limits.end(), parLimits.begin());
   }
};

// List of functions that will be used in the test
vector<fitFunctions> l1DFunctions;
vector<fitFunctions> l2DFunctions;
vector<fitFunctions> treeFunctions;
vector<fitFunctions> l1DLinearFunctions;
vector<fitFunctions> l2DLinearFunctions;

// Gaus 1D implementation
Double_t gaus1DImpl(Double_t* x, Double_t* p)
{
   return p[2]*TMath::Gaus(x[0], p[0], p[1]);
}

// 1D Polynomial implementation
Double_t poly1DImpl(Double_t *x, Double_t *p)
{
   Double_t xx = x[0];
   return p[0]*xx*xx*xx+p[1]*xx*xx+p[2]*xx+p[3];
}

// 2D Polynomial implementation
Double_t poly2DImpl(Double_t *x, Double_t *p)
{
   Double_t xx = x[0];
   Double_t yy = x[1];
   return p[0]*xx*xx*xx+p[1]*xx*xx+p[2]*xx +
          p[3]*yy*yy*yy+p[4]*yy*yy+p[5]*yy +
          p[6];
}

// Gaus 2D Implementation
Double_t gaus2DImpl(Double_t *x, Double_t *p)
{
   return p[0]*TMath::Gaus(x[0], p[1], p[2])*TMath::Gaus(x[1], p[3], p[4]);
}

// Gaus 1D Normalized Implementation
double gausNormal(Double_t* x, Double_t* p)
{
   return p[2]*TMath::Gaus(x[0],p[0],p[1],1);
}

// Gaus 2D Normalized Implementation
double gaus2dnormal(double *x, double *p) {

   double mu_x = p[0];
   double sigma_x = p[1];
   double mu_y = p[2];
   double sigma_y = p[3];
   double rho = p[4];
   double u = (x[0] - mu_x)/ sigma_x ;
   double v = (x[1] - mu_y)/ sigma_y ;
   double c = 1 - rho*rho ;
   double result = (1 / (2 * TMath::Pi() * sigma_x * sigma_y * sqrt(c)))
      * exp (-(u * u - 2 * rho * u * v + v * v ) / (2 * c));
   return result;
}

// N-dimensional Gaus
double gausNd(double *x, double *p) {

   double f = gaus2dnormal(x,p);
   f *= ROOT::Math::normal_pdf(x[2],p[6],p[5]);
   f *= ROOT::Math::normal_pdf(x[3],p[8],p[7]);
   f *= ROOT::Math::normal_pdf(x[4],p[10],p[9]);
   f *= ROOT::Math::normal_pdf(x[5],p[12],p[11]);

   if (f <= 0) {
      std::cout << "invalid f value " << f << " for x ";
      for (int i = 0; i < 6; ++i) std::cout << "  " << x[i];
      std::cout << "\t P = ";
      for (int i = 0; i < 11; ++i) std::cout << "  " << p[i];
      std::cout << "\n\n ";
      return 1.E-300;
   }
   else if (f > 0) return f;

   std::cout << " f is a nan " << f << std::endl;
   for (int i = 0; i < 6; ++i) std::cout << "  " << x[i];
   std::cout << "\t P = ";
   for (int i = 0; i < 11; ++i) std::cout << "  " << p[i];
   std::cout << "\n\n ";
   Error("gausNd","f is a nan");
   assert(1);
   return 0;
}

const double minX = -5.;
const double maxX = +5.;
const double minY = -5.;
const double maxY = +5.;
const int nbinsX = 30;
const int nbinsY = 30;

// Options to indicate how the test has to be run
enum testOpt {
   testOptPars  = 1,  // Check parameters
   testOptChi   = 2,  // Check Chi2 Test
   testOptErr   = 4,  // Show the errors
   testOptColor = 8,  // Show wrong output in color
   testOptDebug = 16, // Print out debug version
   testOptCheck = 32, // Make the checkings
   testOptFitDbg = 64 // Make fit verbose
};

// Default options that all tests will have
int defaultOptions = testOptColor | testOptCheck;// | testOptDebug;

// Object to manage the fitter depending on the optiones used
template <typename T>
class ObjectWrapper {
public:
   T object;
   ObjectWrapper(T _obj): object(_obj) {};
   template <typename F>
   Int_t Fit(F func, const char* opts)
   {
      if ( opts[0] == 'G' )
      {
         ROOT::Fit::BinData d;
         ROOT::Fit::FillData(d,object,func);
//         ROOT::Math::WrappedTF1 f(*func);
         ROOT::Math::WrappedMultiTF1 f(*func);
//          f->SetDerivPrecision(10e-6);
         ROOT::Fit::Fitter fitter;
//          printf("Gradient? FIT?!?\n");
         fitter.Fit(d, f);
         const ROOT::Fit::FitResult & fitResult = fitter.Result();
         // one could set directly the fit result in TF1
         Int_t iret = fitResult.Status();
         if (!fitResult.IsEmpty() ) {
            // set in f1 the result of the fit
            func->SetChisquare( fitResult.Chi2() );
            func->SetNDF( fitResult.Ndf() );
            func->SetNumberFitPoints( d.Size() );

            func->SetParameters( &(fitResult.Parameters().front()) );
            if ( int( fitResult.Errors().size()) >= func->GetNpar() )
               func->SetParErrors( &(fitResult.Errors().front()) );

         }
         // Next line only for debug
//          fitResult.Print(std::cout);
         return iret;
      } else {
//          printf("Normal FIT\n");
         return object->Fit(func, opts);
      }
   };
   const char* GetName() { return object->GetName(); }
};

// Print the Name of the test
int gTestIndex = 0;
template <typename T>
void printTestName(T* object, TF1* func)
{
   string str = "Test  ";
   if (gTestIndex < 10) str += " ";  // add an extra space
   str += ROOT::Math::Util::ToString(gTestIndex);
   str += ":  '";
   str += object->GetName();
   str += "' with '";
   str += func->GetName();
   str += "'...";
   while ( str.length() != 65 )
      str += '.';
   printf("%s", str.c_str());
   fflush(stdout);
}

// In debug mode, prints the title of the debug table.
void printTitle(TF1* func)
{
   printf("\nMin Type    | Min Algo    | OPT  | PARAMETERS             ");
   int n = func->GetNpar();
   for ( int i = 1; i < n; ++i ) {
      printf("                       ");
   }
   printf(" | CHI2TEST        | ERRORS \n");
   fflush(stdout);
}

// In debug mode, separator for the different tests
void printSeparator()
{
   fflush(stdout);
   printf("*********************************************************************"
          "********************************************************************\n");
   fflush(stdout);
}

// Sets the color of the output to red or normal
void setColor(int red = 0)
{
#ifndef R__WIN32
   char command[13];
   if ( red )
      snprintf(command,13, "%c[%d;%d;%dm", 0x1B, 1, 1 + 30, 8 + 40);
   else
      snprintf(command,13, "%c[%dm", 0x1B, 0); // reset to default
   printf("%s", command);
#else
#ifndef __CLING__
   if ( red )
      SetConsoleTextAttribute(GetStdHandle((unsigned long)-11), FOREGROUND_RED );
   else
      SetConsoleTextAttribute(GetStdHandle((unsigned long)-11), FOREGROUND_GREY );
#endif
#endif
}

// Test a fit once it has been done:
//     @str1 Name of the library used
//     @str2 Name of the algorithm used
//     @str3 Options used when fitting
//     @func Fitted function
//     @cmpResult Object to compare the result. It contains all the reference
//               objects as well as the method to compare. It will know whether something has to be tested or not.
//     @opts Options of the test, to know what has to be printed or tested.
int testFit(const char* str1, const char* str2, const char* str3,
               TF1* func, CompareResult const& cmpResult, int opts)
{
   bool debug = opts & testOptDebug;
   // so far, status will just count the number of parameters wronly
   // calculated. There is no other test of the fitters
   int status = 0;
   int diff = 0;

   double chi2 = 0;
   if (  opts & testOptChi || opts & testOptCheck )
      chi2 = func->GetChisquare();

   fflush(stdout);
   if ( opts & testOptPars )
   {
      int n = func->GetNpar();
      double* values = func->GetParameters();
      if ( debug )
         printf("%-11s | %-11s | %-4s | ", str1, str2, str3);
      for ( int i = 0; i < n; ++i ) {
         if ( opts & testOptCheck )
            diff = cmpResult.parameters(i,
                                        values[i],
                                        std::max(std::sqrt(chi2/func->GetNDF()),1.0)*func->GetParError(i));
         status += diff;
         if ( opts & testOptColor )
            setColor ( diff );
         if ( debug )
            printf("%10.6f +/-(%-6.3f) ", values[i], func->GetParError(i));
         fflush(stdout);
      }
      setColor(0);
   }

   if ( opts & testOptChi )
   {
      if ( debug )
         printf(" | chi2: %9.4f | ",  chi2);
   }

   if ( opts & testOptErr )
   {
      // TVirtualFItter is not available in all case (e.g. when running with ROOT IMT)
      if (TVirtualFitter::GetFitter() != 0 ) {
         TBackCompFitter* fitter = dynamic_cast<TBackCompFitter*>( TVirtualFitter::GetFitter() );
         assert(fitter != 0);
         const ROOT::Fit::FitResult& fitResult = fitter->GetFitResult();
         if ( debug )
            printf("err: ");
         int n = func->GetNpar();
         for ( int i = 0; i < n; ++i ) {
            if ( debug )
               printf("%c ", (fitResult.LowerError(i) == fitResult.UpperError(i))?'E':'D');
         }
         if ( debug )
            printf("| ");
      }
   }

   if ( opts != 0 )
   {
      setColor(0);
      if ( debug )
         printf("\n");
   }
   fflush(stdout);

   return status;
}

// Makes all the tests combinations for:
//      @object The object to be fitted
//      @func The function to be used for the fitting
//      @listAlgos All the algorithms that should be tested
//      @fitFunction Parameters of the function used to fill the object
template <typename T, typename F>
int testFitters(T* object, F* func, vector< vector<algoType> >& listAlgos, fitFunctions const& fitFunction)
{
   // counts the number of parameters wronly calculated
   int status = 0;

   int numberOfTests = 0;
   const double* origpars = &(fitFunction.origPars[0]);
   const double* fitpars = &(fitFunction.fitPars[0]);

   func->SetParameters(fitpars);

   printTestName(object, func);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(commonAlgos[0].type, commonAlgos[0].algo);
   ROOT::Math::MinimizerOptions::SetDefaultTolerance(gDefaultTolerance);

   TString opt = "";
   if (! (defaultOptions & testOptFitDbg) ) opt += "Q";
   if (! __DRAW__) opt += "0";

   object->Fit(func, opt);
   if ( defaultOptions & testOptDebug ) printTitle(func);
   struct RefValue ref(origpars, func->GetChisquare());
   commonAlgos[0].cmpResult.setRefValue(&ref);
   int defMinOptions = testOptPars | testOptChi | testOptErr | defaultOptions;
   status += testFit(commonAlgos[0].type, commonAlgos[0].algo
                     , commonAlgos[0].opts, func
                     , commonAlgos[0].cmpResult, defMinOptions);
   numberOfTests += 1;

   if ( defaultOptions & testOptDebug )
   {
      printSeparator();
      func->SetParameters(origpars);
      status += testFit("Parameters", "Original", "", func, commonAlgos[0].cmpResult, testOptPars | testOptDebug);
      func->SetParameters(fitpars);
      status += testFit("Parameters", "Initial",  "", func, commonAlgos[0].cmpResult, testOptPars | testOptDebug);
      printSeparator();
   }

   for ( unsigned int j = 0; j < listAlgos.size(); ++j )
   {
      for ( unsigned int i = 0; i < listAlgos[j].size(); ++i )
      {
         int testFitOptions = testOptPars | testOptChi | testOptErr | defaultOptions;
         ROOT::Math::MinimizerOptions::SetDefaultMinimizer(listAlgos[j][i].type, listAlgos[j][i].algo);
         func->SetParameters(fitpars);
         fflush(stdout);
         object->Fit(func, listAlgos[j][i].opts);
         listAlgos[j][i].cmpResult.setRefValue(&ref);
         status += testFit(listAlgos[j][i].type, listAlgos[j][i].algo, listAlgos[j][i].opts
                           , func, listAlgos[j][i].cmpResult, testFitOptions);
         numberOfTests += 1;
         fflush(stdout);
      }
   }

   double percentageFailure = double( status * 100 ) / double( numberOfTests*func->GetNpar() );

   if ( defaultOptions & testOptDebug )
   {
      printSeparator();
      printf("Number of fails: %d Total Number of tests %d", status, numberOfTests);
      printf(" Percentage of failure: %f\n", percentageFailure );
   }

   // limit in the percentage of failure!
   return (percentageFailure < 4)?0:1;
}

// Test the diferent objects in 1D
int test1DObjects(vector< vector<algoType> >& listH,
                  vector< vector<algoType> >& listG,
                  vector< vector<algoType> >& listGE,
                  vector<fitFunctions>& listOfFunctions)
{
   // Counts how many tests failed.
   int globalStatus = 0;
   // To control if an individual test failed
   int status = 0;

   TF1* func = 0;
   TH1D* h1 = 0;
   TH1D* h2 = 0;
   THnSparse* s1 = 0;
   TGraph* g1 = 0;
   TGraphErrors* ge1 = 0;
   TCanvas *c0 = 0, *c1 = 0, *c2 = 0, *c3 = 0;
   for ( unsigned int j = 0; j < listOfFunctions.size(); ++j )
   {
      if ( func ) delete func;
      func = new TF1( listOfFunctions[j].name, listOfFunctions[j].func, minX, maxX, listOfFunctions[j].npars);
      func->SetParameters(&(listOfFunctions[j].origPars[0]));
      SetParsLimits(listOfFunctions[j].parLimits, func);

      // create here h1 since it is used to make TGrpahs's and THnsparse
      if (h1) delete h1;
      h1 = new TH1D("h1", "Histogram1D Equal Bins", nbinsX, minX, maxX);
      for (int i = 0; i <= h1->GetNbinsX() + 1; ++i)
         h1->SetBinContent(i, rndm.Poisson(func->Eval(h1->GetBinCenter(i))));

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         // fill equal bin 1D  histogram
         TString hname = TString::Format("H1D_%d", gTestIndex);
         h1->SetName(hname);

         if (c1 && !__DRAW__) delete c1;
         c1 = new TCanvas(TString::Format("c%d_H1D",gTestIndex), "Histogram1D");
         ObjectWrapper<TH1D*> owh1(h1);
         globalStatus += status = testFitters(&owh1, func, listH, listOfFunctions[j]);
         if (__DRAW__) {
            h1->DrawCopy();
            func->DrawCopy("SAME");
         }
         printf("%s\n", (status?"FAILED":"OK"));
      }

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         // variable bin test
         double v[nbinsX + 1];
         FillVariableRange(v, nbinsX, minX, maxX);
         if (h2) delete h2;
         TString hname = TString::Format("H1D_%d", gTestIndex);
         h2 = new TH1D(hname, "HIstogram1D Variable Bins", nbinsX, v);
         for (int i = 0; i <= h2->GetNbinsX() + 1; ++i)
            h2->SetBinContent(i, rndm.Poisson(func->Eval(h2->GetBinCenter(i))));

         if (c0 && !__DRAW__) delete c0;
         c0 = new TCanvas(TString::Format("c%d_H1D", gTestIndex), "Histogram1D Variable");
         ObjectWrapper<TH1D *> owh2(h2);
         globalStatus += status = testFitters(&owh2, func, listH, listOfFunctions[j]);
         printf("%s\n", (status ? "FAILED" : "OK"));
         if (__DRAW__) {
            h2->DrawCopy();
            func->DrawCopy("SAME");
         }
      }

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         delete g1; g1 = new TGraph(h1);
         g1->SetName("TGraph1D");   // no need for unique name of TGraphs
         g1->SetTitle("TGraph 1D");
         if (c2 && !__DRAW__) delete c2;
         c2 = new TCanvas(TString::Format("c%d_G1D", gTestIndex), "TGraph");
         ObjectWrapper<TGraph*> owg1(g1);
         globalStatus += status = testFitters(&owg1, func, listG, listOfFunctions[j]);
         printf("%s\n", (status?"FAILED":"OK"));
         if (__DRAW__) {
            g1->DrawClone("AB*");
            func->DrawCopy("SAME");
         }
      }

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         delete ge1; ge1 = new TGraphErrors(h1);
         ge1->SetName("TGraphErrors1D");
         ge1->SetTitle("TGraphErrors 1D");
         if (c3 && !__DRAW__) delete c3;
         c3 = new TCanvas(TString::Format("c%d_G1D", gTestIndex), "TGraphError");
         ObjectWrapper<TGraphErrors*> owge1(ge1);
         globalStatus += status = testFitters(&owge1, func, listGE, listOfFunctions[j]);
         printf("%s\n", (status?"FAILED":"OK"));
         if (__DRAW__) {
            ge1->DrawClone("AB*");
            func->DrawCopy("SAME");
         }
      }

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         delete s1; s1 = THnSparse::CreateSparse("THnSparse1D", "THnSparse 1D", h1);
         ObjectWrapper<THnSparse*> ows1(s1);
         globalStatus += status = testFitters(&ows1, func, listH, listOfFunctions[j]);
         printf("%s\n", (status?"FAILED":"OK"));
      }
   }

   if ( ! __DRAW__ )
   {
      delete func;
      delete h1;
      delete h2;
      delete g1;
      delete ge1;
      delete c0;
      delete c1;
      delete c2;
      delete c3;
   }

   return globalStatus;
}

// Test the different objects in 2S
int test2DObjects(vector< vector<algoType> >& listH,
                  vector< vector<algoType> >& listG,
                  vector< vector<algoType> >& listGE,
                  vector<fitFunctions>& listOfFunctions)
{
   // Counts how many tests failed.
   int globalStatus = 0;
   // To control if an individual test failed
   int status = 0;

   TF2* func = 0;
   TH2D* h1 = 0;
   TH2D* h2 = 0;
   THnSparse* s1 = 0;
   TGraph2D* g1 = 0;
   TGraph2DErrors* ge1 = 0;
   TCanvas *c0 = 0, *c1 = 0, *c2 = 0, *c3 = 0;
   for ( unsigned int h = 0; h < listOfFunctions.size(); ++h )
   {
      if ( func ) delete func;
      func = new TF2( listOfFunctions[h].name, listOfFunctions[h].func, minX, maxX, minY, maxY, listOfFunctions[h].npars);
      func->SetParameters(&(listOfFunctions[h].origPars[0]));
      SetParsLimits(listOfFunctions[h].parLimits, func);

      // fill histogram 2D
      if (h1) delete h1;
      h1 = new TH2D("h2d", "Histogram2D Equal Bins", nbinsX, minX, maxX, nbinsY, minY, maxY);
      if (ge1)
         delete ge1;
      ge1 = new TGraph2DErrors((h1->GetNbinsX() + 1) * (h1->GetNbinsY() + 1));
      ge1->SetName("Graph2DErrors");
      ge1->SetTitle("Graph2D with Errors");
      unsigned int counter = 0;
      for (int i = 0; i <= h1->GetNbinsX() + 1; ++i) {
         for (int j = 0; j <= h1->GetNbinsY() + 1; ++j) {
            double xc = h1->GetXaxis()->GetBinCenter(i);
            double yc = h1->GetYaxis()->GetBinCenter(j);
            double content = rndm.Poisson(func->Eval(xc, yc));
            h1->SetBinContent(i, j, content);
            ge1->SetPoint(counter, xc, yc, content);
            ge1->SetPointError(counter, h1->GetXaxis()->GetBinWidth(i) / 2, h1->GetYaxis()->GetBinWidth(j) / 2,
                               h1->GetBinError(i, j));
            counter += 1;
         }
      }

      // 2D Equal bins test
      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         TString hname = TString::Format("H2D_%d", gTestIndex);
         h1->SetName(hname);
         if ( c1 && ! __DRAW__) delete c1;
         c1 = new TCanvas(TString::Format("c%d_H2D", gTestIndex), "Histogram2D");
         ObjectWrapper<TH2D*> owh1(h1);
         globalStatus += status = testFitters(&owh1, func, listH, listOfFunctions[h]);
         printf("%s\n", (status?"FAILED":"OK"));
         if (__DRAW__) {
            h1->DrawCopy("COLZ");
            func->DrawCopy("SAME");
         }
      }

      // 2D Variable bins test
      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         // fill and test 2D variable bins histograms
         if (h2) delete h2;
         double x[nbinsX + 1];
         FillVariableRange(x, nbinsX, minX, maxX);
         double y[nbinsY + 1];
         FillVariableRange(y, nbinsY, minY, maxY);
         TString hname = TString::Format("H2D_%d", gTestIndex);
         h2 = new TH2D(hname, "Histogram2D Variable Bins", nbinsX, x, nbinsY, y);
         for (int i = 0; i <= h2->GetNbinsX() + 1; ++i) {
            for (int j = 0; j <= h2->GetNbinsY() + 1; ++j) {
               double xc = h2->GetXaxis()->GetBinCenter(i);
               double yc = h2->GetYaxis()->GetBinCenter(j);
               double content = rndm.Poisson(func->Eval(xc, yc));
               h2->SetBinContent(i, j, content);
            }
         }

         if (c0 && ! __DRAW__) delete c0;
         c0 = new TCanvas(TString::Format("c%d_H2D", gTestIndex), "Histogram2D Variable");
         ObjectWrapper<TH2D *> owh2(h2);
         globalStatus += status = testFitters(&owh2, func, listH, listOfFunctions[h]);
         printf("%s\n", (status ? "FAILED" : "OK"));
         if (__DRAW__) {
            h2->DrawCopy("COLZ");
            func->DrawCopy("SAME");
         }
      }

      // TGraph 2D test
      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         if (g1) delete g1;
         g1 = new TGraph2D(h1);
         g1->SetName("TGraph2D");
         g1->SetTitle("TGraph 2D");

         if ( c2 && !__DRAW__) delete c2;
         c2 = new TCanvas(TString::Format("c%d_G2D", gTestIndex), "TGraph");
         ObjectWrapper<TGraph2D*> owg1(g1);
         globalStatus += status = testFitters(&owg1, func, listG, listOfFunctions[h]);
         printf("%s\n", (status?"FAILED":"OK"));
         if (__DRAW__) {
            //g1->DrawClone("AB*");
            g1->DrawClone("surf1");
            func->DrawCopy("SAME");
         }
      }


      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         ge1->SetName("TGraphErrors2DGE");
         ge1->SetTitle("TGraphErrors 2D");

         if (c3 && !__DRAW__) delete c3;
         c3 = new TCanvas(TString::Format("c%d_G2DE", gTestIndex), "TGraphError");
         ObjectWrapper<TGraph2DErrors*> owge1(ge1);
         globalStatus += status = testFitters(&owge1, func, listGE, listOfFunctions[h]);
         printf("%s\n", (status?"FAILED":"OK"));
         if (__DRAW__) {
            ge1->DrawClone("AB*");
            func->DrawCopy("SAME");
         }
      }

      gTestIndex++;
      if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
         delete s1; s1 = THnSparse::CreateSparse("THnSparse2D", "THnSparse 2D", h1);
         ObjectWrapper<THnSparse*> ows1(s1);
         globalStatus += status = testFitters(&ows1, func, listH, listOfFunctions[h]);
         printf("%s\n", (status?"FAILED":"OK"));
      }
   }

   if ( ! __DRAW__ )
   {
      delete func;
      delete h1;
      delete h2;
      delete g1;
      delete ge1;
      delete c0;
      delete c1;
      delete c2;
      delete c3;
   }

   return globalStatus;
}

// Make a wrapper for the TTree, as the interface for fitting
// differs. This way, the same algorithms (testFit and testFitters)
// can be used for all the objects.
class TreeWrapper {
public:

   const char* vars;
   const char* cuts;
   TTree *tree;

   void set(TTree* t, const char* v, const char* c)
   {
      tree = t;
      vars = v;
      cuts = c;
   }

   const char* GetName() const {
      return tree->GetName();
   }

   Int_t Fit(TF1* f1, Option_t* option = "")
   {
      return tree->UnbinnedFit(f1->GetName(), vars, cuts, option);
   }
};

// Test the fittig algorithms for a TTree
int testUnBinnedFit(int n = 10000)
{
   // Counts how many tests failed.
   int globalStatus = 0;
   // To control if an individual test failed
   int status = 0;

   double origPars[13] = {1,2,3,0.5, 0.5, 0, 3, 0, 4, 0, 5, 1, 10 };
//   double fitPars[13] =  {1,1,1,  1, 0.1, 0, 2, 0, 3, 0, 4, 0,  9 };

   TF2 * func = new TF2("gaus2Dn",gaus2dnormal,-10,-10,-10,10,5);
   func->SetParameters(origPars);

   TUnuranMultiContDist dist(func);
   TUnuran unr(&rndm);

   // sampling with vnrou methods
   if (! unr.Init(dist,"vnrou")) {
         std::cerr << "error in init unuran " << std::endl; return -1;
   }

   TTree * tree =  new  TTree("tree","2 var gaus tree");
   double x,y,z,u,v,w;
   tree->Branch("x",&x,"x/D");
   tree->Branch("y",&y,"y/D");
   tree->Branch("z",&z,"z/D");
   tree->Branch("u",&u,"u/D");
   tree->Branch("v",&v,"v/D");
   tree->Branch("w",&w,"w/D");
   double xx[2];
   for (Int_t i=0;i<n;i++) {
      unr.SampleMulti(xx);
      x = xx[0];
      y = xx[1];
      z = rndm.Gaus(origPars[5],origPars[6]);
      u = rndm.Gaus(origPars[7],origPars[8]);
      v = rndm.Gaus(origPars[9],origPars[10]);
      w = rndm.Gaus(origPars[11],origPars[12]);

      tree->Fill();

   }

   delete func;

   vector< vector<algoType> > listAlgos(2);
   listAlgos[0] = commonAlgos;
   listAlgos[1] = simplexAlgos;

   TreeWrapper tw;
   TF1 * f1 = 0;
   TF2 * f2 = 0;
   TF1 * f4 = 0;

   gTestIndex++;
   if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
      f1 = new TF1(treeFunctions[0].name,treeFunctions[0].func,minX,maxY,treeFunctions[0].npars);
      f1->SetParameters( &(treeFunctions[0].fitPars[0]) );
      f1->FixParameter(2,1);
      tw.set(tree, "x", "");
      globalStatus += status = testFitters(&tw, f1, listAlgos, treeFunctions[0]);
      printf("%s\n", (status?"FAILED":"OK"));
   }

   vector<algoType> noCompareInTree;
   // exclude Simplex in tree
   //noCompareInTree.push_back(algoType( "Minuit2",     "Simplex",     "Q0", CompareResult(0)));

   vector< vector<algoType> > listAlgosND(2);
   listAlgosND[0] = commonAlgos;
   listAlgosND[1] = noCompareInTree;

   gTestIndex++;
   if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
      f2 = new TF2(treeFunctions[1].name,treeFunctions[1].func,minX,maxX,minY,maxY,treeFunctions[1].npars);
      f2->SetParameters( &(treeFunctions[1].fitPars[0]) );
      tw.set(tree, "x:y", "");
      globalStatus += status = testFitters(&tw, f2, listAlgosND, treeFunctions[1]);
      printf("%s\n", (status?"FAILED":"OK"));
   }

   gTestIndex++;
   if (gSelectedTest == 0 || gSelectedTest == gTestIndex) {
      f4 = new TF1("gausND",gausNd,0,1,13);
      f4->SetParameters(&(treeFunctions[2].fitPars[0]));
      tw.set(tree, "x:y:z:u:v:w", "");
      globalStatus += status = testFitters(&tw, f4, listAlgosND, treeFunctions[2]);
      printf("%s\n", (status?"FAILED":"OK"));
   }

   delete tree;
   if (f1) delete f1;
   if (f2) delete f2;
   if (f4) delete f4;

   return globalStatus;
}

// Initialize the data for the tests: List of different algorithms and
// fitting functions.
void init_structures()
{
   commonAlgos.push_back( algoType( "Minuit",      "Migrad",      "Q0X", CompareResult())  );
   commonAlgos.push_back( algoType( "Minuit",      "Minimize",    "Q0X", CompareResult())  );
   commonAlgos.push_back( algoType( "Minuit",      "Scan",        "Q0X", CompareResult(0)) );
   commonAlgos.push_back( algoType( "Minuit",      "Seek",        "Q0X", CompareResult())  );
   commonAlgos.push_back( algoType( "Minuit2",     "Migrad",      "Q0X", CompareResult())  );
   commonAlgos.push_back( algoType( "Minuit2",     "Minimize",    "Q0X", CompareResult())  );
   commonAlgos.push_back( algoType( "Minuit2",     "Scan",        "Q0X", CompareResult(0)) );
   commonAlgos.push_back( algoType( "Minuit2",     "Fumili2",     "Q0X", CompareResult())  );
#ifdef R__HAS_MATHMORE
   commonAlgos.push_back( algoType( "GSLMultiMin", "conjugatefr", "Q0X", CompareResult()) );
   commonAlgos.push_back( algoType( "GSLMultiMin", "conjugatepr", "Q0X", CompareResult()) );
   commonAlgos.push_back( algoType( "GSLMultiMin", "bfgs2",       "Q0X", CompareResult()) );
   commonAlgos.push_back( algoType( "GSLSimAn",    "",            "Q0X", CompareResult()) );
#endif

// simplex
   simplexAlgos.push_back( algoType( "Minuit",      "Simplex",     "Q0X", CompareResult()) );
   //simplex MInuit2 does not work well (needs to be checked)
   // simplexAlgos.push_back( algoType( "Minuit2",     "Simplex",     "Q0", CompareResult())  );

   specialAlgos.push_back( algoType( "Minuit",      "Migrad",      "QE0X", CompareResult()) );
   specialAlgos.push_back( algoType( "Minuit",      "Migrad",      "QW0X", CompareResult()) );

   noGraphAlgos.push_back( algoType( "Minuit",      "Migrad","      Q0IX", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Minuit",      "Migrad","      QL0", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Minuit",      "Migrad","     QLI0", CompareResult()) );

   // Gradient algorithms
   // No Minuit algorithms to use with the 'G' options until some stuff is fixed.
   noGraphAlgos.push_back( algoType( "Minuit",      "Migrad",      "GQ0", CompareResult()) );
//    noGraphAlgos.push_back( algoType( "Minuit",      "Minimize",    "GQ0", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Minuit2",     "Migrad",      "GQ0", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Minuit2",     "Minimize",    "GQ0", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Fumili",      "Fumili",      "GQ0", CompareResult()) );
   noGraphAlgos.push_back( algoType( "Minuit2",     "Fumili",      "GQ0", CompareResult()) );
//    noGraphAlgos.push_back( algoType( "Minuit",      "Migrad",      "GQE0", CompareResult()) );
    noGraphAlgos.push_back( algoType( "Minuit",      "Migrad",      "GQL0", CompareResult()) );

   noGraphErrorAlgos.push_back( algoType( "Fumili",      "Fumili",      "Q0", CompareResult()) );
#ifdef R__HAS_MATHMORE
   noGraphErrorAlgos.push_back( algoType( "GSLMultiFit", "",            "Q0", CompareResult()) ); // Not in TGraphError
#endif

   // Same as TH1D (but different comparision scheme!): commonAlgos,
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "Q0",   CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit",      "Minimize",    "Q0",   CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit",      "Scan",        "Q0",   CompareResult(0))         );
   histGaus2D.push_back( algoType( "Minuit",      "Seek",        "Q0",   CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit2",     "Migrad",      "Q0",   CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit2",     "Minimize",    "Q0",   CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit2",     "Scan",        "Q0",   CompareResult(0))         );
   histGaus2D.push_back( algoType( "Minuit2",     "Fumili2",     "Q0",   CompareResult(cmpPars,6)) );
#ifdef R__HAS_MATHMORE
   histGaus2D.push_back( algoType( "GSLMultiMin", "conjugatefr", "Q0", CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "GSLMultiMin", "conjugatepr", "Q0", CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "GSLMultiMin", "bfgs2",       "Q0", CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "GSLSimAn",    "",            "Q0", CompareResult(cmpPars,6)) );
#endif   // treeFail
   histGaus2D.push_back( algoType( "Minuit",      "Simplex",     "Q0",   CompareResult(cmpPars,6)) );
   // minuit2 simplex fails in 2d
   //histGaus2D.push_back( algoType( "Minuit2",     "Simplex",     "Q0",   CompareResult(cmpPars,6)) );
   // special algos
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "QE0",  CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "QW0",  CompareResult())          );
   // noGraphAlgos
   histGaus2D.push_back( algoType( "Minuit",    "Migrad", "       Q0I",  CompareResult(cmpPars,6)) );
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "QL0",  CompareResult())          );
   histGaus2D.push_back( algoType( "Minuit",    "Migrad", "      QLI0", CompareResult())          );

// Gradient algorithms
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "GQ0", CompareResult()) );
   histGaus2D.push_back( algoType( "Minuit",      "Minimize",    "GQ0", CompareResult()) );
   histGaus2D.push_back( algoType( "Minuit2",     "Migrad",      "GQ0", CompareResult()) );
   histGaus2D.push_back( algoType( "Minuit2",     "Minimize",    "GQ0", CompareResult()) );
   histGaus2D.push_back( algoType( "Fumili",      "Fumili",      "GQ0", CompareResult()) );
   histGaus2D.push_back( algoType( "Minuit2",     "Fumili",      "GQ0", CompareResult()) );
//    histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "GQE0", CompareResult()) );
   histGaus2D.push_back( algoType( "Minuit",      "Migrad",      "GQL0", CompareResult()) );

   // noGraphErrorAlgos
   histGaus2D.push_back( algoType( "Fumili",      "Fumili",      "Q0",   CompareResult(cmpPars,6)) );
#ifdef R__HAS_MATHMORE
   histGaus2D.push_back( algoType( "GSLMultiFit", "",            "Q0",   CompareResult(cmpPars,6)) );
#endif

   graphErrorAlgos.push_back( algoType( "Minuit",      "Migrad",      "Q0EX0", CompareResult()) );
   graphErrorAlgos.push_back( algoType( "Minuit2",      "Migrad",      "Q0EX0", CompareResult()) );


   // For testing the linear fitter we can force the use by setting Linear the default minimizer and use
   // the G option. In this case the fit is linearized using the gradient as the linear components
   // Option "G" has not to be set as first option character to avoid using Fitter class in
   // the test program
   // Use option "X" to force Chi2 calculations
   linearAlgos.push_back( algoType( "Linear",      "",            "Q0XG", CompareResult()) );
   listLinearAlgos.push_back( linearAlgos );

   listTH1DAlgos.push_back( commonAlgos );
   listTH1DAlgos.push_back( simplexAlgos );
   listTH1DAlgos.push_back( specialAlgos );
   listTH1DAlgos.push_back( noGraphAlgos );
   listTH1DAlgos.push_back( noGraphErrorAlgos );

   listAlgosTGraph.push_back( commonAlgos );
   listAlgosTGraph.push_back( simplexAlgos );
   listAlgosTGraph.push_back( specialAlgos );
   listAlgosTGraph.push_back( noGraphErrorAlgos );

   listAlgosTGraphError.push_back( commonAlgos );
   listAlgosTGraphError.push_back( simplexAlgos );
   listAlgosTGraphError.push_back( specialAlgos );
   listAlgosTGraphError.push_back( graphErrorAlgos );

   listTH2DAlgos.push_back( histGaus2D );

   listAlgosTGraph2D.push_back( commonAlgos );
   listAlgosTGraph2D.push_back( specialAlgos );
   listAlgosTGraph2D.push_back( noGraphErrorAlgos );

   listAlgosTGraph2DError.push_back( commonAlgos );
   listAlgosTGraph2DError.push_back( specialAlgos );
   listAlgosTGraph2DError.push_back( graphErrorAlgos );

   vector<ParLimit> emptyLimits(0);

   double gausOrig[] = {  0.,  3., 200.};
   double gausFit[] =  {0.5, 3.7,  250.};
   vector<ParLimit> gaus1DLimits;
   gaus1DLimits.push_back( ParLimit(1, 0, 5) );
   l1DFunctions.push_back( fitFunctions("GAUS",       gaus1DImpl, 3, gausOrig,  gausFit, gaus1DLimits) );
   double poly1DOrig[] = { 2, 3, 4, 200};
   double poly1DFit[] = { 6.4, -2.3, 15.4, 210.5};
   l1DFunctions.push_back( fitFunctions("Polynomial", poly1DImpl, 4, poly1DOrig, poly1DFit, emptyLimits) );

   // range os -5,5
   double gaus2DOrig[] = { 500., +.5, 2.7, -.5, 3.0 };
   double gaus2DFit[] = { 510., .0, 1.8, -1.0, 1.6};
   l2DFunctions.push_back( fitFunctions("gaus2D", gaus2DImpl, 5, gaus2DOrig, gaus2DFit, emptyLimits) );

   double gausnOrig[3] = {1,2,1};
   double treeOrig[13] = {1,2,3,0.5, 0.5, 0, 3, 0, 4, 0, 5, 1, 10 };
   double treeFit[13]  = {1,1,1,  1, 0.1, 0, 2, 0, 3, 0, 4, 0,  9 };
   treeFunctions.push_back( fitFunctions("gausn", gausNormal, 3, gausnOrig, treeFit, emptyLimits ));
   treeFunctions.push_back( fitFunctions("gaus2Dn", gaus2dnormal, 5, treeOrig, treeFit, emptyLimits));
   treeFunctions.push_back( fitFunctions("gausND", gausNd, 13, treeOrig, treeFit, emptyLimits));

   l1DLinearFunctions.push_back( fitFunctions("Polynomial", poly1DImpl, 4, poly1DOrig, poly1DFit, emptyLimits) );

   double poly2DOrig[] = { 2, 3, 4, 5, 6, 7, 200, };
   double poly2DFit[] = { 6.4, -2.3, 15.4, 3, 10, -3, 210.5};
   l2DLinearFunctions.push_back( fitFunctions("Poly2D", poly2DImpl, 7, poly2DOrig, poly2DFit, emptyLimits) );
}

int stressHistoFit()
{
   rndm.SetSeed(10);

   if (gEnableMT) ROOT::EnableImplicitMT();

   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");

   init_structures();

   int iret = 0;

   TBenchmark bm;
   bm.Start("stressHistoFit");

   cout << "****************************************************************************" <<endl;
   cout << "*  Starting  stress  H I S T O F I T                                       *" <<endl;
   cout << "****************************************************************************" <<endl;

   std::cout << "\nTest 1D and 2D objects\n\n";
   iret += test1DObjects(listTH1DAlgos, listAlgosTGraph, listAlgosTGraphError, l1DFunctions);
   iret += test2DObjects(listTH2DAlgos, listAlgosTGraph2D, listAlgosTGraph2DError, l2DFunctions);

   std::cout << "\nTest Linear fits\n\n";
   iret += test1DObjects(listLinearAlgos, listLinearAlgos, listLinearAlgos, l1DLinearFunctions);
   iret += test2DObjects(listLinearAlgos, listLinearAlgos, listLinearAlgos, l2DLinearFunctions);
   //defaultOptions = testOptColor | testOptCheck;
   // tree test
   std::cout << "\nTest unbinned fits\n\n";
   iret += testUnBinnedFit(2000);  // reduce statistics

   bm.Stop("stressHistoFit");
   std::cout <<"\n****************************************************************************\n";
   bm.Print("stressHistoFit");
   const double reftime = 124; // ref time on  pcbrun4
   double rootmarks = 800 * reftime / bm.GetCpuTime("stressHistoFit");
   std::cout << " ROOTMARKS = " << rootmarks << " ROOT version: " << gROOT->GetVersion() << "\t"
             << gROOT->GetGitBranch() << "@" << gROOT->GetGitCommit() << std::endl;
   std::cout <<"****************************************************************************\n";

   return iret;
}

int main(int argc, char** argv)
{

   TApplication* theApp = 0;

   Int_t  verbose     =      0;
   Int_t testNumber   =      0;
   Bool_t doDraw      = kFALSE;

   // Parse command line arguments
   for (Int_t i = 1 ;  i < argc ; i++) {

      string arg = argv[i] ;

      if (arg == "-v") {
         cout << "stressHistoFit: running in verbose mode" << endl;
         verbose = 1;
      } else if (arg == "-vv") {
         cout << "stressHistoFit: running in very verbose mode" << endl;
         verbose = 2;
      } else if (arg == "-a") {
         cout << "stressHistoFit: deploying full suite of tests" << endl;
      } else if (arg == "-n") {
         cout << "stressHistoFit: running single test" << endl;
         testNumber = atoi(argv[++i]);
      } else if (arg == "-d") {
         cout << "stressHistoFit: setting gDebug to " << argv[i + 1] << endl;
         gDebug = atoi(argv[++i]);
      } else if (arg == "-g") {
         cout << "stressHistoFit: running in graphics mode " << endl;
         doDraw = kTRUE;
      } else if (arg == "-t") {
         cout << "stressHistoFit: running in multi-thread mode " << endl;
         gEnableMT = kTRUE;
      } else if (arg == "-h") {
         cout << "usage: stressHistoFit [ options ] " << endl;
         cout << "" << endl;
         cout << "       -n N      : only run test with sequential number N" << endl;
         cout << "       -a        : run full suite of tests (default is basic suite); this overrides the -n single test option" << endl;
         cout << "       -g        : create a TApplication and produce plots" << endl;
         cout << "       -v/-vv    : set verbose mode (show result of each regression test) or very verbose mode (show all roofit output as well)" << endl;
         cout << "       -d N      : set ROOT gDebug flag to N" << endl ;
         cout << "       -t        : set ROOT to run in Multi-thread mode" << endl;
         cout << " " << endl ;
         return 0 ;
      }
   }

   __DRAW__ = ( doDraw || verbose == 2);

   if (verbose > 0) {
      defaultOptions |= testOptDebug;   // debug mode (print test results)
      if (verbose > 1) defaultOptions |= testOptFitDbg;   // very debug (print also fit outputs)
   }

   gSelectedTest = testNumber;  // number of selected test

   if ( __DRAW__ )
      theApp = new TApplication("App",&argc,argv);

   int ret = stressHistoFit();

   if ( __DRAW__ ) {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return ret;
}
