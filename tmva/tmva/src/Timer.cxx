// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Timer                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::Timer
\ingroup TMVA
Timing information for training and evaluation of MVA methods

Usage:

~~~ {.cpp}
   TMVA::gConfig().SetDrawProgressBar(true);

   TMVA::Timer timer( Nloops, "MyClassName" );
   for (Int_t i=0; i<Nloops; i++) {
     ... // some code

     // now, print progress bar:
     timer.DrawProgressBar( i );

     // **OR** text output of left time (never both !)
     fLogger << " time left: " << timer.GetLeftTime( i ) << Endl;

   }
   fLogger << "MyClassName" << ": elapsed time: " << timer.GetElapsedTime()
           << Endl;
~~~

Remark: in batch mode, the progress bar is quite ugly; you may
        want to use the text output then

Note that by default in TMVA::Config the drawing of the
progress bar is switched off. To have the progress bar visible you need 
to enable it by calling TMVA::gConfig().SetDrawProgressBar(true)

*/

#include "TMVA/Timer.h"

#include "TMVA/Config.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"

#include "TStopwatch.h"

#ifdef _MSC_VER
#include <io.h>
#define isatty _isatty
#define STDERR_FILENO 2
#else
#include <unistd.h>
#endif

const TString TMVA::Timer::fgClassName = "Timer";
const Int_t   TMVA::Timer::fgNbins     = 16;

ClassImp(TMVA::Timer);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Timer::Timer( const char* prefix, Bool_t colourfulOutput )
   : Timer(0, prefix, colourfulOutput)
{
}

////////////////////////////////////////////////////////////////////////////////
/// standard constructor: ncounts gives the total number of counts that
/// the loop will iterate through. At each call of the timer, the current
/// number of counts is provided by the user, so that the timer can obtain
/// the due time from linearly interpolating the spent time.

TMVA::Timer::Timer( Int_t ncounts, const char* prefix, Bool_t colourfulOutput  )
   : fNcounts        ( ncounts ),
     fPrefix         ( strcmp(prefix,"")==0?Timer::fgClassName:TString(prefix) ),
     fColourfulOutput( colourfulOutput ),
     fPreviousProgress(-1),
     fOutputToFile(!isatty(STDERR_FILENO)),
     fProgressBarStringLength (0),
     fLogger         ( new MsgLogger( fPrefix.Data() ) )
{
   fColourfulOutput = fColourfulOutput && !fOutputToFile;
   Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Timer::~Timer( void )
{
   delete fLogger;
}

void TMVA::Timer::Init( Int_t ncounts )
{
   // timer initialisation
   fNcounts = ncounts;
   Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// resets timer

void TMVA::Timer::Reset( void )
{
   TStopwatch::Start( kTRUE );
   fPreviousProgress = -1;
   fPreviousTimeEstimate.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// computes elapsed tim in seconds

Double_t TMVA::Timer::ElapsedSeconds( void )
{
   Double_t rt = TStopwatch::RealTime(); TStopwatch::Start( kFALSE );
   return rt;
}

////////////////////////////////////////////////////////////////////////////////
/// returns pretty string with elapsed time

TString TMVA::Timer::GetElapsedTime( Bool_t Scientific )
{
   return SecToText( ElapsedSeconds(), Scientific );
}

////////////////////////////////////////////////////////////////////////////////
/// returns pretty string with time left

TString TMVA::Timer::GetLeftTime( Int_t icounts )
{
   Double_t leftTime = ( icounts <= 0 ? -1 :
                         icounts > fNcounts ? -1 :
                         Double_t(fNcounts - icounts)/Double_t(icounts)*ElapsedSeconds() );

   return SecToText( leftTime, kFALSE );
}

////////////////////////////////////////////////////////////////////////////////
/// draws the progressbar

void TMVA::Timer::DrawProgressBar()
{
   fProgressBarStringLength = 0;
   fNcounts++;
   if (fNcounts == 1) {
      std::clog << fLogger->GetPrintedSource();
      std::clog << "Please wait ";
   }

   std::clog << "." << std::flush;
}

////////////////////////////////////////////////////////////////////////////////
/// draws a string in the progress bar

void TMVA::Timer::DrawProgressBar( TString theString )
{

   std::clog << fLogger->GetPrintedSource();

   std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << "[" << gTools().Color("reset");

   std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << theString << gTools().Color("reset");

   std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << "]" << gTools().Color("reset");

   for (int i = fProgressBarStringLength; i < theString.Length (); ++i)
      std::cout << " ";
   std::clog << "\r" << std::flush;
   fProgressBarStringLength = theString.Length ();
}

////////////////////////////////////////////////////////////////////////////////
/// draws progress bar in color or B&W
/// caution:

void TMVA::Timer::DrawProgressBar( Int_t icounts, const TString& comment  )
{
   if (!gConfig().DrawProgressBar()) return;

   // sanity check:
   if (icounts > fNcounts-1) icounts = fNcounts-1;
   if (icounts < 0         ) icounts = 0;
   Int_t ic = Int_t(Float_t(icounts)/Float_t(fNcounts)*fgNbins);

   auto timeLeft = this->GetLeftTime( icounts );

   // do not redraw progress bar when neither time not ticks are different
   if (ic == fPreviousProgress && timeLeft == fPreviousTimeEstimate && icounts != fNcounts-1) return;
   // check if we are redirected to a file
   if (fOutputToFile) {
       if (ic != fPreviousProgress) {
           std::clog << Int_t((100*(icounts+1))/Float_t(fNcounts)) << "%, time left: " << timeLeft << std::endl;
           fPreviousProgress = ic;
       }
       return;
   }
   fPreviousProgress = ic;
   fPreviousTimeEstimate = timeLeft;

   std::clog << fLogger->GetPrintedSource();
   if (fColourfulOutput) std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << "[" << gTools().Color("reset");
   else                  std::clog << "[";
   for (Int_t i=0; i<ic; i++) {
      if (fColourfulOutput) std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << ">" << gTools().Color("reset");
      else                  std::clog << ">";
   }
   for (Int_t i=ic+1; i<fgNbins; i++) {
      if (fColourfulOutput) std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << "." << gTools().Color("reset");
      else                  std::clog << ".";
   }
   if (fColourfulOutput) std::clog << gTools().Color("white_on_green") << gTools().Color("dyellow") << "]" << gTools().Color("reset");
   else                  std::clog << "]" ;

   // timing information
   if (fColourfulOutput) {
      std::clog << gTools().Color("reset") << " " ;
      std::clog << "(" << gTools().Color("red") << Int_t((100*(icounts+1))/Float_t(fNcounts)) << "%" << gTools().Color("reset")
                << ", "
                << "time left: "
                << timeLeft << gTools().Color("reset") << ") ";
   }
   else {
      std::clog << "] " ;
      std::clog << "(" << Int_t((100*(icounts+1))/Float_t(fNcounts)) << "%" 
                << ", " << "time left: " << timeLeft << ") ";
   }
   if (comment != "") {
      std::clog << "[" << comment << "]  ";
   }
   std::clog << "\r" << std::flush;
}

////////////////////////////////////////////////////////////////////////////////
/// pretty string output

TString TMVA::Timer::SecToText( Double_t seconds, Bool_t Scientific ) const
{
   TString out = "";
   if      (Scientific    ) out = Form( "%.3g sec", seconds );
   else if (seconds <  0  ) out = "unknown";
   else if (seconds <= 300) out = Form( "%i sec", Int_t(seconds) );
   else {
      if (seconds > 3600) {
         Int_t h = Int_t(seconds/3600);
         if (h <= 1) out = Form( "%i hr : ", h );
         else        out = Form( "%i hrs : ", h );

         seconds = Int_t(seconds)%3600;
      }
      Int_t m = Int_t(seconds/60);
      if (m <= 1) out += Form( "%i min", m );
      else        out += Form( "%i mins", m );
   }

   return (fColourfulOutput) ? gTools().Color("red") + out + gTools().Color("reset") : out;
}
