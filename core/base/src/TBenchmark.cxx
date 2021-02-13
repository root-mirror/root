
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBenchmark.h"
#include "TStopwatch.h"

TBenchmark *gBenchmark = nullptr;

ClassImp(TBenchmark);

/** \class TBenchmark
\ingroup Base

This class is a ROOT utility to help benchmarking applications
*/

////////////////////////////////////////////////////////////////////////////////
/// Benchmark default constructor

TBenchmark::TBenchmark(): TNamed()
{
   fNbench   = 0;
   fNmax     = 20;
   fNames    = 0;
   fRealTime = 0;
   fCpuTime  = 0;
   fTimer    = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TBenchmark::TBenchmark(const TBenchmark& bm) :
  TNamed(bm),
  fNbench(bm.fNbench),
  fNmax(bm.fNmax),
  fNames(0),
  fRealTime(0),
  fCpuTime(0),
  fTimer(0)
{
   fNames    = new TString[fNmax];
   fRealTime = new Float_t[fNmax];
   fCpuTime  = new Float_t[fNmax];
   fTimer    = new TStopwatch[fNmax];

   for(Int_t i = 0; i<fNmax; ++i) {
      fNames[i] = bm.fNames[i];
      fRealTime[i] = bm.fRealTime[i];
      fCpuTime[i] = bm.fCpuTime[i];
      fTimer[i] = bm.fTimer[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TBenchmark& TBenchmark::operator=(const TBenchmark& bm)
{
   if (this!=&bm) {
      TNamed::operator=(bm);
      fNbench=bm.fNbench;
      fNmax=bm.fNmax;

      delete [] fNames;
      delete [] fRealTime;
      delete [] fCpuTime;
      delete [] fTimer;

      fNames    = new TString[fNmax];
      fRealTime = new Float_t[fNmax];
      fCpuTime  = new Float_t[fNmax];
      fTimer    = new TStopwatch[fNmax];

      for(Int_t i = 0; i<fNmax; ++i) {
         fNames[i] = bm.fNames[i];
         fRealTime[i] = bm.fRealTime[i];
         fCpuTime[i] = bm.fCpuTime[i];
         fTimer[i] = bm.fTimer[i];
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Benchmark destructor.

TBenchmark::~TBenchmark()
{
   fNbench   = 0;
   if (fNames)    { delete [] fNames;    fNames  = 0;}
   if (fRealTime) { delete [] fRealTime; fRealTime  = 0;}
   if (fCpuTime)  { delete [] fCpuTime;  fCpuTime  = 0;}
   if (fTimer  )  { delete [] fTimer;    fTimer  = 0;}
}

////////////////////////////////////////////////////////////////////////////////
/// Returns index of Benchmark name.

Int_t TBenchmark::GetBench(const char *name) const
{
   for (Int_t i=0;i<fNbench;i++) {
      if (!strcmp(name,(const char*)fNames[i])) return i;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Cpu time used by Benchmark name.

Float_t TBenchmark::GetCpuTime(const char *name)
{
   Int_t bench = GetBench(name);
   if (bench >= 0) return fCpuTime[bench];
   else            return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Realtime used by Benchmark name.

Float_t TBenchmark::GetRealTime(const char *name)
{
   Int_t bench = GetBench(name);
   if (bench >= 0) return fRealTime[bench];
   else            return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints parameters of Benchmark name.

void TBenchmark::Print(const char *name) const
{
   Int_t bench = GetBench(name);
   if (bench < 0) return;
   Printf("%-10s: Real Time = %6.2f seconds Cpu Time = %6.2f seconds",name,fRealTime[bench],fCpuTime[bench]);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset all Benchmarks

void TBenchmark::Reset()
{
   fNbench = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stops Benchmark name and Prints results

void TBenchmark::Show(const char *name)
{
   Stop(name);
   Print((char*)name);
}

////////////////////////////////////////////////////////////////////////////////
/// Starts Benchmark with the specified name.
///
/// An independent timer (see class TStopwatch) is started.
/// The name of the benchmark is entered into the list of benchmarks.
/// Benchmark can be stopped via TBenchmark::Stop().
/// Results can be printed via TBenchmark::Print().
/// TBenchmark::Show() can be used to stop benchmark and print results.
/// If name is an already existing benchmark, timing will resume.
/// A summary of all benchmarks can be seen via TBenchmark::Summary().

void TBenchmark::Start(const char *name)
{
   if (!fNames) {
      fNames    = new TString[fNmax];
      fRealTime = new Float_t[fNmax];
      fCpuTime  = new Float_t[fNmax];
      fTimer    = new TStopwatch[fNmax];
   }
   Int_t bench = GetBench(name);
   if (bench < 0 && fNbench < fNmax ) {
      // define a new benchmark to Start
      fNames[fNbench] = name;
      bench = fNbench;
      fNbench++;
      fTimer[bench].Reset();
      fTimer[bench].Start();
      fRealTime[bench] = 0;
      fCpuTime[bench]  = 0;
   } else if (bench >= 0) {
      // Resume the existing benchmark
      fTimer[bench].Continue();
   }
   else
      Warning("Start","too many benchmarks");
}

////////////////////////////////////////////////////////////////////////////////
/// Terminates Benchmark with specified name.

void TBenchmark::Stop(const char *name)
{
   Int_t bench = GetBench(name);
   if (bench < 0) return;

   fTimer[bench].Stop();
   fRealTime[bench] = fTimer[bench].RealTime();
   fCpuTime[bench]  = fTimer[bench].CpuTime();
}

////////////////////////////////////////////////////////////////////////////////
/// Prints a summary of all benchmarks.

void TBenchmark::Summary(Float_t &rt, Float_t &cp)
{
   rt = 0;
   cp = 0;
   for (Int_t i=0;i<fNbench;i++) {
      Printf("%-10s: Real Time = %6.2f seconds Cpu Time = %6.2f seconds",(const char*)fNames[i],fRealTime[i],fCpuTime[i]);
      rt += fRealTime[i];
      cp += fCpuTime[i];
   }
   Printf("%-10s: Real Time = %6.2f seconds Cpu Time = %6.2f seconds","TOTAL",rt,cp);
}
