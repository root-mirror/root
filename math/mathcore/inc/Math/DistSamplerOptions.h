// @(#)root/mathcore:$Id$
// Author: L. Moneta Fri Aug 15 2008

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2008  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Math_DistSamplerOptions
#define ROOT_Math_DistSamplerOptions

#include <string>

#include <iostream>
#include "Math/IOptions.h"

namespace ROOT {


namespace Math {


//_______________________________________________________________________________
/**
    DistSampler options class

    @ingroup NumAlgo
*/
class DistSamplerOptions {

public:

   // static methods for setting and retrieving the default options

   static void SetDefaultSampler(const char * type);
   static void SetDefaultAlgorithm1D(const char * algo );
   static void SetDefaultAlgorithmND(const char * algo );
   static void SetDefaultPrintLevel(int level);

   static const std::string & DefaultSampler();
   static const std::string & DefaultAlgorithm1D();
   static const std::string & DefaultAlgorithmND();
   static int DefaultPrintLevel();

   /// retrieve extra options - if not existing create a IOptions
   static ROOT::Math::IOptions & Default(const char * name);

   // find extra options - return 0 if not existing
   static ROOT::Math::IOptions * FindDefault(const char * name);

   /// print all the default options for the name given
   static void PrintDefault(const char * name, std::ostream & os = std::cout);

public:

   // constructor using the default options
   // pass optionally a pointer to the additional options
   // otherwise look if they exist for this default minimizer
   // and in that case they are copied in the constructed instance
   // constructor takes dimension since a different default algorithm
   // is used if the dimension is 1 or greater than 1
   DistSamplerOptions(int dim = 0);

   // destructor
   ~DistSamplerOptions();

   // copy constructor
   DistSamplerOptions(const DistSamplerOptions & opt);

   /// assignment operators
   DistSamplerOptions & operator=(const DistSamplerOptions & opt);

   /** non-static methods for  retrivieng options */

   /// set print level
   int PrintLevel() const { return fLevel; }

   /// return extra options (NULL pointer if they are not present)
   IOptions * ExtraOptions() const { return fExtraOptions; }

   /// type of minimizer
   const std::string & Sampler() const { return fSamplerType; }

   /// type of algorithm (method)
   const std::string & Algorithm() const { return fAlgoType; }

   /// print all the options
   void Print(std::ostream & os = std::cout) const;

   /** non-static methods for setting options */

   /// set print level
   void SetPrintLevel(int level) { fLevel = level; }

   /// set minimizer type
   void SetSampler(const char * type) { fSamplerType = type; }

   /// set minimizer algorithm
   void SetAlgorithm(const char *type) { fAlgoType = type; }

   /// set extra options (in this case pointer is cloned)
   void  SetExtraOptions(const IOptions & opt);

   /// set a specific algorithm option
   template <class T>
   void SetAlgoOption(const char * name, T value) {
      if (!fExtraOptions) CreateExtraOptions();
      fExtraOptions->SetValue(name, value);
   }

protected:

   void CreateExtraOptions();

private:

   int fLevel;               // debug print level
   std::string fSamplerType;   // DistSampler type (Unuran, Foam, etc...)xs
   std::string fAlgoType;    // DistSampler algorithmic specification (for Unuran only)

   // extra options
   ROOT::Math::IOptions *   fExtraOptions;  // extra options

};

   } // end namespace Math

} // end namespace ROOT

#endif
