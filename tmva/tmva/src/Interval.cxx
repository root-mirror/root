// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::Interval                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 * File and Version Information:                                                  *
 **********************************************************************************/

/*! \class TMVA::Interval
\ingroup TMVA

The TMVA::Interval Class

Interval definition, continuous and discrete

  - Interval(min,max)  : a continous interval [min,max]
  - Interval(min,max,n): a "discrete interval" [min,max], i.e the n numbers:
         min, min+step, min+2*step,...., min+(n-1)*step, min+n*step=max

  e.g.:

  - Interval(1,5,5) = 1,2,3,4,5
  - Interval(.5,1.,6) = .5, .6., .7, .8, .9, 1.0

 Note: **bin** counting starts from ZERO unlike in ROOT histograms

  - Interval definition, continuous and discrete

    - Interval(min,max)  : a continous interval [min,max]
    - Interval(min,max,n): a "discrete interval" [min,max], i.e the n numbers:

      min, min+step, min+2*step,...., min+(n-1)*step=max

      e.g.:

      - Interval(1,5,5)=1,2,3,4,5                    <br>
      - Interval(.5,1.,6)= .5, .6., .7, .8, .9, 1.0        <br>

~~~ {.cpp}
   Example:   Interval(.5,1.,6)

             [ min                           max ]
         -----------------------------------------------
                |     |     |     |     |     |
               .5    .6    .7    .8    .9    1.0

         bin    0     1     2     3     4     5
~~~
*/

#include "TRandom3.h"
#include "ThreadLocalStorage.h"

#include "TMVA/Interval.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

ClassImp(TMVA::Interval);

////////////////////////////////////////////////////////////////////////////////
/// defines minimum and maximum of an interval
///  - when nbins > 0, interval describes a discrete distribution (equally distributed in the interval)
///  - when nbins == 0, interval describes a continous interval

TMVA::Interval::Interval( Double_t min, Double_t max, Int_t nbins ) :
fMin(min),
   fMax(max),
   fNbins(nbins)
{
   if (fMax - fMin < 0) Log() << kFATAL << "maximum lower than minimum" << Endl;
   if (nbins < 0) {
      Log() << kFATAL << "nbins < 0" << Endl;
      return;
   }
   else if (nbins == 1) {
      Log() << kFATAL << "interval has to have at least 2 bins if discrete" << Endl;
      return;
   }
}

TMVA::Interval::Interval( const Interval& other ) :
   fMin  ( other.fMin ),
   fMax  ( other.fMax ),
   fNbins( other.fNbins )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Interval::~Interval()
{
}

////////////////////////////////////////////////////////////////////////////////
/// calculates the value of the "number" bin in a discrete interval.
/// Parameters:
///        Double_t position
///

Double_t TMVA::Interval::GetElement( Int_t bin ) const
{
   if (fNbins <= 0) {
      Log() << kFATAL << "GetElement only defined for discrete value Intervals" << Endl;
      return 0.0;
   }
   else if (bin < 0 || bin >= fNbins) {
      Log() << kFATAL << "bin " << bin << " out of range: interval *bins* count from 0 to " << fNbins-1  << Endl;
      return 0.0;
   }
   return fMin + ( (Double_t(bin)/(fNbins-1)) *(fMax - fMin) );
}

////////////////////////////////////////////////////////////////////////////////
/// returns the step size between the numbers of a "discrete Interval"

Double_t TMVA::Interval::GetStepSize( Int_t iBin )  const
{
   if (fNbins <= 0) {
      Log() << kFATAL << "GetElement only defined for discrete value Intervals" << Endl;
   }
   if (iBin<0) {
      Log() << kFATAL << "You asked for iBin=" << iBin
            <<" in interval .. and.. sorry, I cannot let this happen.."<<Endl;
   }
   return (fMax-fMin)/(Double_t)(fNbins-1);
}

////////////////////////////////////////////////////////////////////////////////
/// get uniformly distributed number within interval

Double_t TMVA::Interval::GetRndm( TRandom3& rnd )  const
{
   return rnd.Rndm()*(fMax - fMin) + fMin;
}

Double_t TMVA::Interval::GetWidth() const
{
   return fMax - fMin;
}
Double_t TMVA::Interval::GetMean()  const
{
   return (fMax + fMin)/2;
}

void TMVA::Interval::Print(std::ostream &os) const
{
   for (Int_t i=0; i<GetNbins(); i++){
      os << "| " << GetElement(i)<<" |" ;
   }
}

TMVA::MsgLogger& TMVA::Interval::Log() const {
   TTHREAD_TLS_DECL_ARG(MsgLogger,logger,"Interval");   // message logger
   return logger;
}
