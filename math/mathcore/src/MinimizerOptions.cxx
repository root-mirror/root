// @(#)root/mathcore:$Id$ 
// Author: L. Moneta Fri Aug 15 2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2008  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Math/MinimizerOptions.h"

#include "Math/GenAlgoOptions.h"

// case of using ROOT plug-in manager
#ifndef MATH_NO_PLUGIN_MANAGER
#include "TEnv.h"
#endif 

#include <iomanip>

namespace ROOT { 
   

namespace Math { 

   namespace Minim { 
      static std::string gDefaultMinimizer = ""; // take from /etc/system.rootrc in ROOT Fitter
      static std::string gDefaultMinimAlgo = "Migrad";
      static double gDefaultErrorDef = 1.;
      static double gDefaultTolerance = 1.E-2; 
      static double gDefaultPrecision = -1; // value <= 0 means left to minimizer
      static int  gDefaultMaxCalls = 0; // 0 means leave default values Deaf
      static int  gDefaultMaxIter  = 0; 
      static int  gDefaultStrategy  = 1; 
      static int  gDefaultPrintLevel  = 0; 
   }


void MinimizerOptions::SetDefaultMinimizer(const char * type, const char * algo ) {   
   // set the default minimizer type and algorithm
   if (type) Minim::gDefaultMinimizer = std::string(type); 
   if (algo) Minim::gDefaultMinimAlgo = std::string(algo);
}
void MinimizerOptions::SetDefaultErrorDef(double up) {
   // set the default error definition 
   Minim::gDefaultErrorDef = up; 
}
void MinimizerOptions::SetDefaultTolerance(double tol) {
   // set the defult tolerance
   Minim::gDefaultTolerance = tol; 
}
void MinimizerOptions::SetDefaultPrecision(double prec) {
   // set the defult precision
   Minim::gDefaultPrecision = prec; 
}
void MinimizerOptions::SetDefaultMaxFunctionCalls(int maxcall) {
   // set the default maximum number of function calls
   Minim::gDefaultMaxCalls = maxcall; 
}
void MinimizerOptions::SetDefaultMaxIterations(int maxiter) {
   // set the default maximum number of iterations
   Minim::gDefaultMaxIter = maxiter; 
}
void MinimizerOptions::SetDefaultStrategy(int stra) {
   // set the default minimization strategy
   Minim::gDefaultStrategy = stra; 
}
void MinimizerOptions::SetDefaultPrintLevel(int level) {
   // set the default printing level 
   Minim::gDefaultPrintLevel = level; 
}

const std::string & MinimizerOptions::DefaultMinimizerAlgo() { return Minim::gDefaultMinimAlgo; }
double MinimizerOptions::DefaultErrorDef()         { return Minim::gDefaultErrorDef; }
double MinimizerOptions::DefaultTolerance()        { return Minim::gDefaultTolerance; }
double MinimizerOptions::DefaultPrecision()        { return Minim::gDefaultPrecision; }
int    MinimizerOptions::DefaultMaxFunctionCalls() { return Minim::gDefaultMaxCalls; }
int    MinimizerOptions::DefaultMaxIterations()    { return Minim::gDefaultMaxIter; }
int    MinimizerOptions::DefaultStrategy()         { return Minim::gDefaultStrategy; }
int    MinimizerOptions::DefaultPrintLevel()       { return Minim::gDefaultPrintLevel; }

const std::string & MinimizerOptions::DefaultMinimizerType() 
{ 
   // return default minimizer
   // if is "" (no default is set) read from etc/system.rootrc

   if (Minim::gDefaultMinimizer.size() == 0) { 
#ifndef MATH_NO_PLUGIN_MANAGER
   // use value defined in etc/system.rootrc  (if not found Minuit is used) 
      if (gEnv) 
         Minim::gDefaultMinimizer = gEnv->GetValue("Root.Fitter","Minuit");   
#else
      Minim::gDefaultMinimizer = "Minuit2";  // in case no PM exists 
#endif
   }

   return Minim::gDefaultMinimizer; 
}


MinimizerOptions::MinimizerOptions(IOptions * extraOpts): 
   fLevel( Minim::gDefaultPrintLevel),
   fMaxCalls( Minim::gDefaultMaxCalls ), 
   fMaxIter( Minim::gDefaultMaxIter ), 
   fStrategy( Minim::gDefaultStrategy ), 
   fErrorDef(  Minim::gDefaultErrorDef ), 
   fTolerance( Minim::gDefaultTolerance ),
   fPrecision( Minim::gDefaultPrecision ),
   fExtraOptions(extraOpts)
{
   // constructor using  the default options

   fMinimType = MinimizerOptions::DefaultMinimizerType();

   fAlgoType =  Minim::gDefaultMinimAlgo;

   // case of Fumili2 and TMinuit
   if (fMinimType == "TMinuit") fMinimType = "Minuit";
   else if (fMinimType == "Fumili2") { 
      fMinimType = "Minuit2";
      fAlgoType = "Fumili";
   }
   else if (fMinimType.find("cmaes")!=std::string::npos
	    || fMinimType.find("ipop")!=std::string::npos)
     {
       fAlgoType = fMinimType;
       fMinimType = "cmaes";
     }
   else if (fMinimType == "GSLMultiMin" && fAlgoType == "Migrad") 
      fAlgoType = "BFGS2";
   
   
   // check if extra options exists (copy them if needed)
   if (!fExtraOptions) { 
      IOptions * gopts = FindDefault( fMinimType.c_str() );
      if (gopts) fExtraOptions = gopts->Clone();
   }
}


MinimizerOptions::MinimizerOptions(const MinimizerOptions & opt) : fExtraOptions(0) {  
   // copy constructor 
   (*this) = opt; 
}

MinimizerOptions & MinimizerOptions::operator=(const MinimizerOptions & opt) {  
   // assignment operator 
   if (this == &opt) return *this; // self assignment
   fLevel = opt.fLevel;
   fMaxCalls = opt.fMaxCalls; 
   fMaxIter = opt.fMaxIter; 
   fStrategy = opt.fStrategy; 
   fErrorDef = opt.fErrorDef;
   fTolerance = opt.fTolerance;
   fPrecision = opt.fPrecision; 
   fMinimType = opt.fMinimType; 
   fAlgoType = opt.fAlgoType; 

   if (fExtraOptions) delete fExtraOptions; 
   fExtraOptions = 0; 
   if (opt.fExtraOptions)  fExtraOptions =  (opt.fExtraOptions)->Clone();
   return *this;
}

MinimizerOptions::~MinimizerOptions() { 
   if (fExtraOptions) delete fExtraOptions; 
}

void MinimizerOptions::SetExtraOptions(const IOptions & opt) {  
   // set extra options (clone the passed one)
   if (fExtraOptions) delete fExtraOptions; 
   fExtraOptions = opt.Clone(); 
}

void MinimizerOptions::Print(std::ostream & os) const {
   //print all the options
   os << std::setw(25) << "Minimizer Type"        << " : " << std::setw(15) << fMinimType << std::endl;
   os << std::setw(25) << "Minimizer Algorithm"   << " : " << std::setw(15) << fAlgoType << std::endl;
   os << std::setw(25) << "Strategy"              << " : " << std::setw(15) << fStrategy << std::endl;
   os << std::setw(25) << "Tolerance"              << " : " << std::setw(15) << fTolerance << std::endl;
   os << std::setw(25) << "Max func calls"         << " : " << std::setw(15) << fMaxCalls << std::endl;
   os << std::setw(25) << "Max iterations"         << " : " << std::setw(15) << fMaxIter << std::endl;
   os << std::setw(25) << "Func Precision"         << " : " << std::setw(15) << fPrecision << std::endl;
   os << std::setw(25) << "Error definition"       << " : " << std::setw(15) << fErrorDef << std::endl;
   os << std::setw(25) << "Print Level"            << " : " << std::setw(15) << fLevel << std::endl;
   
   if (ExtraOptions()) { 
      os << fMinimType << " specific options :"  << std::endl;
      ExtraOptions()->Print(os);
   }
}

IOptions & MinimizerOptions::Default(const char * name) { 
   // create default extra options for the given algorithm type 
   return GenAlgoOptions::Default(name);
}

IOptions * MinimizerOptions::FindDefault(const char * name) { 
   // find extra options for the given algorithm type 
   return GenAlgoOptions::FindDefault(name);
}

void MinimizerOptions::PrintDefault(const char * name, std::ostream & os) {
   //print default options
   MinimizerOptions tmp;
   tmp.Print(os);
   if (!tmp.ExtraOptions() ) {
      IOptions * opt = FindDefault(name);
      os << "Specific options for "  << name << std::endl;
      if (opt) opt->Print(os);
   }
}




} // end namespace Math

} // end namespace ROOT

