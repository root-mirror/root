// @(#)root/test:$Id$
// Author: Rene Brun   19/08/96

////////////////////////////////////////////////////////////////////////
//
//                       Event and Track classes
//                       =======================
//
//  The Event class is a naive/simple example of an event structure.
//     public:
//        char           fType[20];
//        char          *fEventName;         //run+event number in character format
//        Int_t          fNtrack;
//        Int_t          fNseg;
//        Int_t          fNvertex;
//        UInt_t         fFlag;
//        Double32_t     fTemperature;
//        Int_t          fMeasures[10];
//        Double32_t     fMatrix[4][4];
//        Double32_t    *fClosestDistance; //[fNvertex] indexed array!
//        EventHeader    fEvtHdr;
//        TClonesArray  *fTracks;
//        TRefArray     *fHighPt;            //array of High Pt tracks only
//        TRefArray     *fMuons;             //array of Muon tracks only
//        TRef           fLastTrack;         //pointer to last track
//        TRef           fHistoWeb;          //EXEC:GetHistoWeb reference to an histogram in a TWebFile
//        TH1F          *fH;
//        TBits          fTriggerBits;       //Bits triggered by this event.
//
//   The EventHeader class has 3 data members (integers):
//     public:
//        Int_t          fEvtNum;
//        Int_t          fRun;
//        Int_t          fDate;
//
//
//   The Event data member fTracks is a pointer to a TClonesArray.
//   It is an array of a variable number of tracks per event.
//   Each element of the array is an object of class Track with the members:
//     private:
//        Float_t      fPx;           //X component of the momentum
//        Float_t      fPy;           //Y component of the momentum
//        Float_t      fPz;           //Z component of the momentum
//        Float_t      fRandom;       //A random track quantity
//        Float_t      fMass2;        //The mass square of this particle
//        Float_t      fBx;           //X intercept at the vertex
//        Float_t      fBy;           //Y intercept at the vertex
//        Float_t      fMeanCharge;   //Mean charge deposition of all hits of this track
//        Float_t      fXfirst;       //X coordinate of the first point
//        Float_t      fXlast;        //X coordinate of the last point
//        Float_t      fYfirst;       //Y coordinate of the first point
//        Float_t      fYlast;        //Y coordinate of the last point
//        Float_t      fZfirst;       //Z coordinate of the first point
//        Float_t      fZlast;        //Z coordinate of the last point
//        Double32_t   fCharge;       //Charge of this track
//        Double32_t   fVertex[3];    //Track vertex position
//        Int_t        fNpoint;       //Number of points for this track
//        Short_t      fValid;        //Validity criterion
//        Int_t        fNsp;          //Number of points for this track with a special value
//        Double32_t  *fPointValue;   //[fNsp] a special quantity for some point.
//        TBits        fTriggerBits;  //Bits triggered by this track.
//
//   An example of a batch program to use the Event/Track classes is given
//   in this directory: MainEvent.
//   Look also in the same directory at the following macros:
//     - eventa.C  an example how to read the tree
//     - eventb.C  how to read events conditionally
//
//   During the processing of the event (optionally) also a large number
//   of histograms can be filled. The creation and handling of the
//   histograms is taken care of by the HistogramManager class.
//
//   Note:  This version of the class Event (see EventMT.h and EventMT.cxx
//   for an alternative) uses static variables to improve performance (by
//   reducing the number of memory allocations).  Consequently, only one
//   instance of the class Event should be in use at a time (a 2nd instance
//   would share the array of Tracks with the first instance).
//
////////////////////////////////////////////////////////////////////////

#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"
#include "snprintf.h"

#include "Event.h"


ClassImp(EventHeader);
ClassImp(Event);
ClassImp(Track);
ClassImp(HistogramManager);

TClonesArray *Event::fgTracks = 0;
TH1F *Event::fgHist = 0;

////////////////////////////////////////////////////////////////////////////////
/// Create an Event object.
/// When the constructor is invoked for the first time, the class static
/// variable fgTracks is 0 and the TClonesArray fgTracks is created.

Event::Event() : fIsValid(kFALSE)
{
   if (!fgTracks) fgTracks = new TClonesArray("Track", 1000);
   fTracks = fgTracks;
   fHighPt = new TRefArray;
   fMuons  = new TRefArray;
   fNtrack = 0;
   fH      = 0;
   Int_t i0,i1;
   for (i0 = 0; i0 < 4; i0++) {
      for (i1 = 0; i1 < 4; i1++) {
         fMatrix[i0][i1] = 0.0;
      }
   }
   for (i0 = 0; i0 <10; i0++) fMeasures[i0] = 0;
   for (i0 = 0; i0 <20; i0++) fType[i0] = 0;
   fClosestDistance = 0;
   fEventName = 0;
   fWebHistogram.SetAction(this);
}

////////////////////////////////////////////////////////////////////////////////

Event::~Event()
{
   Clear();
   if (fH == fgHist) fgHist = 0;
   delete fH; fH = 0;
   delete fHighPt; fHighPt = 0;
   delete fMuons;  fMuons = 0;
   delete [] fClosestDistance;
   if (fEventName) delete [] fEventName;
}

////////////////////////////////////////////////////////////////////////////////

void Event::Build(Int_t ev, Int_t arg5, Float_t ptmin) {
  fIsValid = kTRUE;
  char etype[20];
  Float_t sigmat, sigmas;
  gRandom->Rannor(sigmat,sigmas);
  Int_t ntrack   = Int_t(arg5 +arg5*sigmat/120.);
  Float_t random = gRandom->Rndm();

  //Save current Object count
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  Clear();
  fHighPt->Delete();
  fMuons->Delete();

  Int_t nch = 15;
  if (ev >= 100)   nch += 3;
  if (ev >= 10000) nch += 3;
  if (fEventName) delete [] fEventName;
  fEventName = new char[nch];
  snprintf(fEventName,nch,"Event%d_Run%d",ev,200);
  snprintf(etype,20,"type%d",ev%5);
  SetType(etype);
  SetHeader(ev, 200, 960312, random);
  SetNseg(Int_t(10*ntrack+20*sigmas));
  SetNvertex(Int_t(1+20*gRandom->Rndm()));
  SetFlag(UInt_t(random+0.5));
  SetTemperature(random+20.);

  for(UChar_t m = 0; m < 10; m++) {
     SetMeasure(m, Int_t(gRandom->Gaus(m,m+1)));
  }
  for(UChar_t i0 = 0; i0 < 4; i0++) {
    for(UChar_t i1 = 0; i1 < 4; i1++) {
       SetMatrix(i0,i1,gRandom->Gaus(i0*i1,1));
    }
  }

  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));

  //  Create and Fill the Track objects
  for (Int_t t = 0; t < ntrack; t++) AddTrack(random,ptmin);

  //Restore Object count
  //To save space in the table keeping track of all referenced objects
  //we assume that our events do not address each other. We reset the
  //object count to what it was at the beginning of the event.
  TProcessID::SetObjectCount(ObjectNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new track to the list of tracks for this event.
/// To avoid calling the very time consuming operator new for each track,
/// the standard but not well know C++ operator "new with placement"
/// is called. If tracks[i] is 0, a new Track object will be created
/// otherwise the previous Track[i] will be overwritten.

Track *Event::AddTrack(Float_t random, Float_t ptmin)
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,32,0)
   Track *track = (Track*)fTracks->ConstructedAt(fNtrack++);
   track->Set(random);
#else
   TClonesArray &tracks = *fTracks;
   Track *track = new(tracks[fNtrack++]) Track(random);
#endif
   //Save reference to last Track in the collection of Tracks
   fLastTrack = track;
   //Save reference in fHighPt if track is a high Pt track
   if (track->GetPt() > ptmin)   fHighPt->Add(track);
   //Save reference in fMuons if track is a muon candidate
   if (track->GetMass2() < 0.11) fMuons->Add(track);
   return track;
}

////////////////////////////////////////////////////////////////////////////////

void Event::Clear(Option_t * /*option*/)
{
   fTracks->Clear("C"); //will also call Track::Clear
   fHighPt->Delete();
   fMuons->Delete();
   fTriggerBits.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to reset all static objects for this event
///   fgTracks->Delete(option);

void Event::Reset(Option_t * /*option*/)
{
   delete fgTracks; fgTracks = 0;
   fgHist   = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Event::SetHeader(Int_t i, Int_t run, Int_t date, Float_t random)
{
   fNtrack = 0;
   fEvtHdr.Set(i, run, date);
   if (!fgHist) fgHist = new TH1F("hstat","Event Histogram",100,0,1);
   fH = fgHist;
   fH->Fill(random);
}

////////////////////////////////////////////////////////////////////////////////

void Event::SetMeasure(UChar_t which, Int_t what) {
   if (which<10) fMeasures[which] = what;
}

////////////////////////////////////////////////////////////////////////////////
/// This delete is to test the relocation of variable length array

void Event::SetRandomVertex() {
   if (fClosestDistance) delete [] fClosestDistance;
   if (!fNvertex) {
      fClosestDistance = 0;
      return;
   }
   fClosestDistance = new Double32_t[fNvertex];
   for (Int_t k = 0; k < fNvertex; k++ ) {
      fClosestDistance[k] = gRandom->Gaus(1,1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a track object

Track::Track(const Track &orig) : TObject(orig),fTriggerBits(orig.fTriggerBits)
{
   fPx = orig.fPx;
   fPy = orig.fPy;
   fPz = orig.fPx;
   fRandom = orig.fRandom;
   fMass2 = orig.fMass2;
   fBx = orig.fBx;
   fBy = orig.fBy;
   fMeanCharge = orig.fMeanCharge;
   fXfirst = orig.fXfirst;
   fXlast  = orig.fXlast;
   fYfirst = orig.fYfirst;
   fYlast  = orig.fYlast;
   fZfirst = orig.fZfirst;
   fZlast  = orig.fZlast;
   fCharge = orig.fCharge;

   fVertex[0] = orig.fVertex[0];
   fVertex[1] = orig.fVertex[1];
   fVertex[2] = orig.fVertex[2];
   fNpoint = orig.fNpoint;
   fNsp = orig.fNsp;
   if (fNsp) {
      fPointValue = new Double32_t[fNsp];
      for(int i=0; i<fNsp; i++) {
         fPointValue[i] = orig.fPointValue[i];
      }
   } else {
      fPointValue = 0;
   }
   fValid  = orig.fValid;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a track object.
/// Note that in this example, data members do not have any physical meaning.

Track::Track(Float_t random) : TObject(),fTriggerBits(64)
{
   Float_t a,b,px,py;
   gRandom->Rannor(px,py);
   fPx = px;
   fPy = py;
   fPz = TMath::Sqrt(px*px+py*py);
   fRandom = 1000*random;
   if (fRandom < 10) fMass2 = 0.106;
   else if (fRandom < 100) fMass2 = 0.8;
   else if (fRandom < 500) fMass2 = 4.5;
   else if (fRandom < 900) fMass2 = 8.9;
   else  fMass2 = 9.8;
   gRandom->Rannor(a,b);
   fBx = 0.1*a;
   fBy = 0.1*b;
   fMeanCharge = 0.01*gRandom->Rndm();
   gRandom->Rannor(a,b);
   fXfirst = a*10;
   fXlast  = b*10;
   gRandom->Rannor(a,b);
   fYfirst = a*12;
   fYlast  = b*16;
   gRandom->Rannor(a,b);
   fZfirst = 50 + 5*a;
   fZlast  = 200 + 10*b;
   fCharge = Double32_t(Int_t(3*gRandom->Rndm()) - 1);

   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));

   fVertex[0] = gRandom->Gaus(0,0.1);
   fVertex[1] = gRandom->Gaus(0,0.2);
   fVertex[2] = gRandom->Gaus(0,10);
   fNpoint = Int_t(60+10*gRandom->Rndm());
   fNsp = Int_t(3*gRandom->Rndm());
   if (fNsp) {
      fPointValue = new Double32_t[fNsp];
      for(int i=0; i<fNsp; i++) {
         fPointValue[i] = i+1;
      }
   } else {
      fPointValue = 0;
   }
   fValid  = Int_t(0.6+gRandom->Rndm());
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a track

Track &Track::operator=(const Track &orig)
{
   TObject::operator=(orig);
   fPx = orig.fPx;
   fPy = orig.fPy;
   fPz = orig.fPx;
   fRandom = orig.fRandom;
   fMass2 = orig.fMass2;
   fBx = orig.fBx;
   fBy = orig.fBy;
   fMeanCharge = orig.fMeanCharge;
   fXfirst = orig.fXfirst;
   fXlast  = orig.fXlast;
   fYfirst = orig.fYfirst;
   fYlast  = orig.fYlast;
   fZfirst = orig.fZfirst;
   fZlast  = orig.fZlast;
   fCharge = orig.fCharge;

   fVertex[0] = orig.fVertex[0];
   fVertex[1] = orig.fVertex[1];
   fVertex[2] = orig.fVertex[2];
   fNpoint = orig.fNpoint;
   if (fNsp > orig.fNsp) {
      fNsp = orig.fNsp;
      if (fNsp == 0) {
         delete [] fPointValue;
         fPointValue = 0;
      } else {
         for(int i=0; i<fNsp; i++) {
            fPointValue[i] = orig.fPointValue[i];
         }
      }
   } else {
      if (fNsp) {
         delete [] fPointValue;
      }
      fNsp = orig.fNsp;
      if (fNsp) {
         fPointValue = new Double32_t[fNsp];
         for(int i=0; i<fNsp; i++) {
            fPointValue[i] = orig.fPointValue[i];
         }
      } else {
         fPointValue = 0;
      }
   }
   fValid  = orig.fValid;

   fTriggerBits = orig.fTriggerBits;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Note that we intend on using TClonesArray::ConstructedAt, so we do not
/// need to delete any of the arrays.

void Track::Clear(Option_t * /*option*/)
{
   TObject::Clear();
   fTriggerBits.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the values of the Track data members.

void Track::Set(Float_t random)
{
   Float_t a,b,px,py;
   gRandom->Rannor(px,py);
   fPx = px;
   fPy = py;
   fPz = TMath::Sqrt(px*px+py*py);
   fRandom = 1000*random;
   if (fRandom < 10) fMass2 = 0.106;
   else if (fRandom < 100) fMass2 = 0.8;
   else if (fRandom < 500) fMass2 = 4.5;
   else if (fRandom < 900) fMass2 = 8.9;
   else  fMass2 = 9.8;
   gRandom->Rannor(a,b);
   fBx = 0.1*a;
   fBy = 0.1*b;
   fMeanCharge = 0.01*gRandom->Rndm();
   gRandom->Rannor(a,b);
   fXfirst = a*10;
   fXlast  = b*10;
   gRandom->Rannor(a,b);
   fYfirst = a*12;
   fYlast  = b*16;
   gRandom->Rannor(a,b);
   fZfirst = 50 + 5*a;
   fZlast  = 200 + 10*b;
   fCharge = Double32_t(Int_t(3*gRandom->Rndm()) - 1);

   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm()));

   fVertex[0] = gRandom->Gaus(0,0.1);
   fVertex[1] = gRandom->Gaus(0,0.2);
   fVertex[2] = gRandom->Gaus(0,10);
   fNpoint = Int_t(60+10*gRandom->Rndm());
   Int_t newNsp = Int_t(3*gRandom->Rndm());
   if (fNsp > newNsp) {
      fNsp = newNsp;
      if (fNsp == 0) {
         delete [] fPointValue;
         fPointValue = 0;
      } else {
         for(int i=0; i<fNsp; i++) {
            fPointValue[i] = i+1;
         }
      }

   } else {
      if (fNsp) {
         delete [] fPointValue;
      }
      fNsp = newNsp;
      if (fNsp) {
         fPointValue = new Double32_t[fNsp];
         for(int i=0; i<fNsp; i++) {
            fPointValue[i] = i+1;
         }
      } else {
         fPointValue = 0;
      }
   }
   fValid  = Int_t(0.6+gRandom->Rndm());
}

////////////////////////////////////////////////////////////////////////////////
/// Create histogram manager object. Histograms will be created
/// in the "dir" directory.

HistogramManager::HistogramManager(TDirectory *dir)
{
   // Save current directory and cd to "dir".
   TDirectory *saved = gDirectory;
   dir->cd();

   fNtrack      = new TH1F("hNtrack",    "Ntrack",100,575,625);
   fNseg        = new TH1F("hNseg",      "Nseg",100,5800,6200);
   fTemperature = new TH1F("hTemperature","Temperature",100,19.5,20.5);
   fPx          = new TH1F("hPx",        "Px",100,-4,4);
   fPy          = new TH1F("hPy",        "Py",100,-4,4);
   fPz          = new TH1F("hPz",        "Pz",100,0,5);
   fRandom      = new TH1F("hRandom",    "Random",100,0,1000);
   fMass2       = new TH1F("hMass2",     "Mass2",100,0,12);
   fBx          = new TH1F("hBx",        "Bx",100,-0.5,0.5);
   fBy          = new TH1F("hBy",        "By",100,-0.5,0.5);
   fMeanCharge  = new TH1F("hMeanCharge","MeanCharge",100,0,0.01);
   fXfirst      = new TH1F("hXfirst",    "Xfirst",100,-40,40);
   fXlast       = new TH1F("hXlast",     "Xlast",100,-40,40);
   fYfirst      = new TH1F("hYfirst",    "Yfirst",100,-40,40);
   fYlast       = new TH1F("hYlast",     "Ylast",100,-40,40);
   fZfirst      = new TH1F("hZfirst",    "Zfirst",100,0,80);
   fZlast       = new TH1F("hZlast",     "Zlast",100,0,250);
   fCharge      = new TH1F("hCharge",    "Charge",100,-1.5,1.5);
   fNpoint      = new TH1F("hNpoint",    "Npoint",100,50,80);
   fValid       = new TH1F("hValid",     "Valid",100,0,1.2);

   // cd back to original directory
   saved->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up all histograms.

HistogramManager::~HistogramManager()
{
   // Nothing to do. Histograms will be deleted when the directory
   // in which tey are stored is closed.
}

////////////////////////////////////////////////////////////////////////////////
/// Fill histograms.

void HistogramManager::Hfill(Event *event)
{
   fNtrack->Fill(event->GetNtrack());
   fNseg->Fill(event->GetNseg());
   fTemperature->Fill(event->GetTemperature());

   for (Int_t itrack = 0; itrack < event->GetNtrack(); itrack++) {
      Track *track = (Track*)event->GetTracks()->UncheckedAt(itrack);
      fPx->Fill(track->GetPx());
      fPy->Fill(track->GetPy());
      fPz->Fill(track->GetPz());
      fRandom->Fill(track->GetRandom());
      fMass2->Fill(track->GetMass2());
      fBx->Fill(track->GetBx());
      fBy->Fill(track->GetBy());
      fMeanCharge->Fill(track->GetMeanCharge());
      fXfirst->Fill(track->GetXfirst());
      fXlast->Fill(track->GetXlast());
      fYfirst->Fill(track->GetYfirst());
      fYlast->Fill(track->GetYlast());
      fZfirst->Fill(track->GetZfirst());
      fZlast->Fill(track->GetZlast());
      fCharge->Fill(track->GetCharge());
      fNpoint->Fill(track->GetNpoint());
      fValid->Fill(track->GetValid());
   }
}
