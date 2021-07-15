// @(#)root/mathcore:$Id$
// Author: L. Moneta Fri Aug 15 2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2008  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#include "Math/DistSamplerOptions.h"

#include "Math/GenAlgoOptions.h"


#include <iomanip>

namespace ROOT {


namespace Math {

   namespace Sampler {
      static std::string gDefaultSampler = "Unuran";
      static std::string gDefaultAlgorithm1D = "auto";
      static std::string gDefaultAlgorithmND = "vnrou";
      static int  gDefaultPrintLevel  = 0;
   }


void DistSamplerOptions::SetDefaultSampler(const char * type ) {
   // set the default minimizer type and algorithm
   if (type) Sampler::gDefaultSampler = std::string(type);
}
void DistSamplerOptions::SetDefaultAlgorithm1D(const char * algo ) {
   // set the default minimizer type and algorithm
   if (algo) Sampler::gDefaultAlgorithm1D = std::string(algo);
}
void DistSamplerOptions::SetDefaultAlgorithmND(const char * algo ) {
   // set the default minimizer type and algorithm
   if (algo) Sampler::gDefaultAlgorithmND = std::string(algo);
}
void DistSamplerOptions::SetDefaultPrintLevel(int level) {
   // set the default printing level
   Sampler::gDefaultPrintLevel = level;
}

const std::string & DistSamplerOptions::DefaultAlgorithm1D() { return Sampler::gDefaultAlgorithm1D; }
const std::string & DistSamplerOptions::DefaultAlgorithmND() { return Sampler::gDefaultAlgorithmND; }
int    DistSamplerOptions::DefaultPrintLevel()       { return Sampler::gDefaultPrintLevel; }

const std::string & DistSamplerOptions::DefaultSampler()
{
   // return default minimizer
   // if is "" (no default is set) read from etc/system.rootrc
   // use form /etc/ ??

   return Sampler::gDefaultSampler;
}


DistSamplerOptions::DistSamplerOptions(int dim):
   fLevel( Sampler::gDefaultPrintLevel),
   fExtraOptions(0)
{
   // constructor using  the default options

   fSamplerType = DistSamplerOptions::DefaultSampler();

   if (dim == 1)
      fAlgoType =  DistSamplerOptions::DefaultAlgorithm1D();
   else if (dim >1)
      fAlgoType =  DistSamplerOptions::DefaultAlgorithmND();
   else
      // not specified - keep null string
      fAlgoType = std::string();

   // check if extra options exists (copy them if needed)
   if (!fExtraOptions) {
      IOptions * gopts = FindDefault( fSamplerType.c_str() );
      if (gopts) fExtraOptions = gopts->Clone();
   }
}


DistSamplerOptions::DistSamplerOptions(const DistSamplerOptions & opt) : fExtraOptions(0) {
   // copy constructor
   (*this) = opt;
}

DistSamplerOptions & DistSamplerOptions::operator=(const DistSamplerOptions & opt) {
   // assignment operator
   if (this == &opt) return *this; // self assignment
   fLevel = opt.fLevel;
   fSamplerType = opt.fSamplerType;
   fAlgoType = opt.fAlgoType;

   if (fExtraOptions) delete fExtraOptions;
   fExtraOptions = 0;
   if (opt.fExtraOptions)  fExtraOptions =  (opt.fExtraOptions)->Clone();
   return *this;
}

DistSamplerOptions::~DistSamplerOptions() {
   if (fExtraOptions) delete fExtraOptions;
}

void DistSamplerOptions::SetExtraOptions(const IOptions & opt) {
   // set extra options (clone the passed one)
   if (fExtraOptions) delete fExtraOptions;
   fExtraOptions = opt.Clone();
}

void DistSamplerOptions::CreateExtraOptions() { 
   if (fExtraOptions) return;
   fExtraOptions = new ROOT::Math::GenAlgoOptions();
}

void DistSamplerOptions::Print(std::ostream & os) const {
   //print all the options
   os << std::setw(25) << "DistSampler Type"        << " : " << std::setw(15) << fSamplerType << std::endl;
   os << std::setw(25) << "DistSampler Algorithm"   << " : " << std::setw(15) << fAlgoType << std::endl;
   os << std::setw(25) << "Print Level"            << " : " << std::setw(15) << fLevel << std::endl;

   if (ExtraOptions()) {
      os << fSamplerType << " specific options :"  << std::endl;
      ExtraOptions()->Print(os);
   }
}

IOptions & DistSamplerOptions::Default(const char * name) {
   // create default extra options for the given algorithm type
   return GenAlgoOptions::Default(name);
}

IOptions * DistSamplerOptions::FindDefault(const char * name) {
   // find extra options for the given algorithm type
   return GenAlgoOptions::FindDefault(name);
}

void DistSamplerOptions::PrintDefault(const char * name, std::ostream & os) {
   //print default options
   os << "Default DistSampler options " << std::endl;
   os << std::setw(25) << "Default  Type"        << " : " << std::setw(15) << DistSamplerOptions::DefaultSampler() << std::endl;
   os << std::setw(25) << "Default Algorithm 1D"   << " : " << std::setw(15) << DistSamplerOptions::DefaultAlgorithm1D() << std::endl;
   os << std::setw(25) << "Default Algorithm ND"   << " : " << std::setw(15) << DistSamplerOptions::DefaultAlgorithmND() << std::endl;
   os << std::setw(25) << "Default Print Level"    << " : " << std::setw(15) << DistSamplerOptions::DefaultPrintLevel() << std::endl;
   IOptions * opt = FindDefault(name);
   if (opt) {
      os << "Specific default options for "  << name << std::endl;
      opt->Print(os);
   }
}




} // end namespace Math

} // end namespace ROOT

