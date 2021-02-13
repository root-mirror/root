// @(#)root/hist:$Id$
// Author: Christophe.Delaere@cern.ch   21/08/2002

///////////////////////////////////////////////////////////////////////////
/** \class TLimitDataSource
  This class serves as interface to feed data into the TLimit routines
*/
///////////////////////////////////////////////////////////////////////////

#include "TLimitDataSource.h"
#include "TH1.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "snprintf.h"

ClassImp(TLimitDataSource);

TLimitDataSource::TLimitDataSource()
{
   // Default constructor
   fDummyTA.SetOwner();
   fDummyIds.SetOwner();
}

TLimitDataSource::TLimitDataSource(TH1 * s, TH1 * b, TH1 * d)
{
   // Another constructor, directly adds one channel
   // with signal, background and data given as input.
   fDummyTA.SetOwner();
   fDummyIds.SetOwner();
   AddChannel(s, b, d);
}

TLimitDataSource::TLimitDataSource(TH1 * s, TH1 * b, TH1 * d,
                                   TVectorD * es, TVectorD * eb, TObjArray * names)
{
   // Another constructor, directly adds one channel
   // with signal, background and data given as input.
   fDummyTA.SetOwner();
   fDummyIds.SetOwner();
   AddChannel(s, b, d, es, eb, names);
}

void TLimitDataSource::AddChannel(TH1 * s, TH1 * b, TH1 * d)
{
   // Adds a channel with signal, background and data given as input.

   TVectorD *empty;
   TRandom3 generator;
   fSignal.AddLast(s);
   fBackground.AddLast(b);
   fCandidates.AddLast(d);
   char rndname[20];
   snprintf(rndname,20, "rndname%f", generator.Rndm());
   empty = new TVectorD(1);
   fErrorOnSignal.AddLast(empty);
   fDummyTA.AddLast(empty);
   snprintf(rndname,20, "rndname%f", generator.Rndm());
   empty = new TVectorD(1);
   fErrorOnBackground.AddLast(empty);
   fDummyTA.AddLast(empty);
   TObjArray *dummy = new TObjArray(0);
   fIds.AddLast(dummy);
   fDummyIds.AddLast(dummy);
}

void TLimitDataSource::AddChannel(TH1 * s, TH1 * b, TH1 * d, TVectorD * es,
                                  TVectorD * eb, TObjArray * names)
{
   // Adds a channel with signal, background and data given as input.
   // In addition, error sources are defined.
   // TH1 are here used for convenience: each bin has to be seen as
   // an error source (relative).
   // names is an array of strings containing the names of the sources.
   // Sources with the same name are correlated.

   fSignal.AddLast(s);
   fBackground.AddLast(b);
   fCandidates.AddLast(d);
   fErrorOnSignal.AddLast(es);
   fErrorOnBackground.AddLast(eb);
   fIds.AddLast(names);
}

void TLimitDataSource::SetOwner(bool swtch)
{
   // Gives to the TLimitDataSource the ownership of the various objects
   // given as input.
   // Objects are then deleted by the TLimitDataSource destructor.

   fSignal.SetOwner(swtch);
   fBackground.SetOwner(swtch);
   fCandidates.SetOwner(swtch);
   fErrorOnSignal.SetOwner(swtch);
   fErrorOnBackground.SetOwner(swtch);
   fIds.SetOwner(swtch);
   fDummyTA.SetOwner(!swtch);
   fDummyIds.SetOwner(!swtch);
}

