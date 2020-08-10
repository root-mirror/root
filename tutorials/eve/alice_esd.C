/// \file
/// \ingroup tutorial_eve
/// Complex example showing ALICE ESD track visualization.
///
///  alice_esd.C - a simple event-display for ALICE ESD tracks and clusters
///
///
///  Only standard ROOT is used to process the ALICE ESD files.
///
///  No ALICE code is needed, only four simple coordinate-transformation
///  functions declared in this macro.
///
///  A simple geometry of 10KB, extracted from the full TGeo-geometry, is
///  used to outline the central detectors of ALICE.
///
///  All files are access from the web by using the "CACHEREAD" option.
///
///
///  ### Automatic building of ALICE ESD class declarations and dictionaries.
///
///  ALICE ESD is a TTree containing tracks and other event-related
///  information with one entry per event. All these classes are part of
///  the AliROOT offline framework and are not available to standard
///  ROOT.
///
///  To be able to access the event data in a natural way, by using
///  data-members of classes and object containers, the header files and
///  class dictionaries are automatically generated from the
///  TStreamerInfo classes stored in the ESD file by using the
///  TFile::MakeProject() function. The header files and a shared library
///  is created in the aliesd/ directory and can be loaded dynamically
///  into the ROOT session.
///
///  See the run_alice_esd.C macro.
///
///
///  ### Creation of simple GUI for event navigation.
///
///  Most common use of the event-display is to browse through a
///  collection of events. Thus a simple GUI allowing this is created in
///  the function make_gui().
///
///  Eve uses the configurable ROOT-browser as its main window and so we
///  create an extra tab in the left working area of the browser and
///  provide backward/forward buttons.
///
///
///  ### Event-navigation functions.
///
///  As this is a simple macro, we store the information about the
///  current event in the global variable 'Int_t esd_event_id'. The
///  functions for event-navigation simply modify this variable and call
///  the load_event() function which does the following:
///  1. drop the old visualization objects;
///  2. retrieve given event from the ESD tree;
///  3. call alice_esd_read() function to create visualization objects
///     for the new event.
///
///
///  ### Reading of ALICE data and creation of visualization objects.
///
///  This is performed in alice_esd_read() function, with the following
///  steps:
///  1. create the track container object - TEveTrackList;
///  2. iterate over the ESD tracks, create TEveTrack objects and append
///     them to the container;
///  3. instruct the container to extrapolate the tracks and set their
///     visual attributes.
///
/// \image html eve_alice_esd.png
/// \macro_code
///
/// \author Matevz Tadel

#ifndef __RUN_ALICE_ESD__

void alice_esd()
{
   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("alice_esd.C","");
   dir.ReplaceAll("/./","/");
   gROOT->LoadMacro(dir +"MultiView.C+");
   const char* esd_file_name = "http://root.cern/files/alice_ESDs.root";
   TFile::SetCacheFileDir(".");
   TString lib(Form("aliesd/aliesd.%s", gSystem->GetSoExt()));

   if (gSystem->AccessPathName(lib, kReadPermission)) {
      TFile* f = TFile::Open(esd_file_name, "CACHEREAD");
      if (f == 0) return;
      TTree *tree = (TTree*) f->Get("esdTree");
      tree->SetBranchStatus ("ESDfriend*", 1);
      f->MakeProject("aliesd", "*", "++");
      f->Close();
      delete f;
   }
   gSystem->Load(lib);
   gROOT->ProcessLine("#define __RUN_ALICE_ESD__ 1");
   gROOT->ProcessLine("#include \"alice_esd.C\"");
   gROOT->ProcessLine("run_alice_esd()");
   gROOT->ProcessLine("#undef __RUN_ALICE_ESD__");
}

#else

#include "alice_esd_html_summary.C"

// Forward declarations.

class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class AliExternalTrackParam;

void       make_gui();
void       load_event();

void       alice_esd_read();
TEveTrack* esd_make_track(TEveTrackPropagator* trkProp, Int_t index,
                          AliESDtrack* at,
                          AliExternalTrackParam* tp=0);
Bool_t     trackIsOn(AliESDtrack* t, Int_t mask);
void       trackGetPos(AliExternalTrackParam* tp, Double_t r[3]);
void       trackGetMomentum(AliExternalTrackParam* tp, Double_t p[3]);
Double_t   trackGetP(AliExternalTrackParam* tp);


// Configuration and global variables.

const char* esd_file_name = "http://root.cern/files/alice_ESDs.root";
// Temporarily disable reading of ESD friend.
// There seems to be no way to get it working without AliRoot.
// const char* esd_friends_file_name =
//       "http://root.cern/files/alice_ESDfriends.root";
const char* esd_friends_file_name = 0;

const char* esd_geom_file_name =
   "http://root.cern/files/alice_ESDgeometry.root";

// For testing
// const char* esd_file_name         = "AliESDs.root";
// const char* esd_friends_file_name = "AliESDfriends.root";

TFile *esd_file          = 0;
TFile *esd_friends_file  = 0;

TTree *esd_tree          = 0;

AliESDEvent  *esd        = 0;
TList        *esd_objs   = 0;
AliESDfriend *esd_friend = 0;

Int_t esd_event_id       = 0; // Current event id.

TEveTrackList *gTrackList = 0;

TEveGeoShape *gGeomGentle = 0;

// Implemented in MultiView.C
class MultiView;
MultiView* gMultiView = 0;

/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

//______________________________________________________________________________
void run_alice_esd()
{
   // Main function, initializes the application.
   //
   // 1. Load the auto-generated library holding ESD classes and
   //    ESD dictionaries.
   // 2. Open ESD data-files.
   // 3. Load cartoon geometry.
   // 4. Spawn simple GUI.
   // 5. Load first event.

   const TString weh("alice_esd()");

   TFile::SetCacheFileDir(".");

   printf("*** Opening ESD ***\n");
   esd_file = TFile::Open(esd_file_name, "CACHEREAD");
   if (!esd_file)
      return;

   esd_tree = (TTree*)       esd_file->Get("esdTree");
   esd      = (AliESDEvent*) esd_tree->GetUserInfo()->FindObject("AliESDEvent");
   esd_objs = esd->fESDObjects;

   if (esd_friends_file_name != 0)
   {
      printf("*** Opening ESD-friends ***\n");
      esd_friends_file = TFile::Open(esd_friends_file_name, "CACHEREAD");
      if (!esd_friends_file)
         return;

      esd_tree->SetBranchStatus ("ESDfriend*", 1);
   }

   // Set the branch addresses.
   {
      TIter next(esd_objs);
      TObject *el;
      while ((el = (TNamed*)next()))
      {
         TString bname(el->GetName());
         if (bname == "AliESDfriend")
         {
            // AliESDfriend needs special treatment.
            TBranch *br = esd_tree->GetBranch("ESDfriend.");
            br->SetAddress(esd_objs->GetObjectRef(el));
         }
         else
         {
            TBranch *br = esd_tree->GetBranch(bname);
            if (br)
            {
               br->SetAddress(esd_objs->GetObjectRef(el));
            }
            else
            {
               br = esd_tree->GetBranch(bname + ".");
               if (br)
               {
                  br->SetAddress(esd_objs->GetObjectRef(el));
               }
               else
               {
                  Warning("AliESDEvent::ReadFromTree() "
                          "No Branch found with Name '%s' or '%s.'.",
                          bname.Data(),bname.Data());
               }
            }
         }
      }
   }


   TEveManager::Create();

   { // Simple geometry
      TFile* geom = TFile::Open(esd_geom_file_name, "CACHEREAD");
      if (!geom)
         return;
      TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) geom->Get("Gentle");
      gGeomGentle = TEveGeoShape::ImportShapeExtract(gse, 0);
      geom->Close();
      delete geom;
      gEve->AddGlobalElement(gGeomGentle);
   }


   // Standard multi-view
   //=====================

   gMultiView = new MultiView;

   gMultiView->ImportGeomRPhi(gGeomGentle);
   gMultiView->ImportGeomRhoZ(gGeomGentle);


   // HTML summary view
   //===================

   fgHtmlSummary = new HtmlSummary("Alice Event Display Summary Table");
   slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
   fgHtml = new TGHtml(0, 100, 100);
   TEveWindowFrame *wf = slot->MakeFrame(fgHtml);
   fgHtml->MapSubwindows();
   wf->SetElementName("Summary");


   // Final stuff
   //=============

   gEve->GetBrowser()->GetTabRight()->SetTab(1);

   make_gui();

   load_event();

   gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
}

//______________________________________________________________________________
void load_event()
{
   // Load event specified in global esd_event_id.
   // The contents of previous event are removed.

   printf("Loading event %d.\n", esd_event_id);

   gEve->GetViewers()->DeleteAnnotations();

   if (gTrackList)
      gTrackList->DestroyElements();

   esd_tree->GetEntry(esd_event_id);
   // esd_tree->Show();

   alice_esd_read();

   TEveElement* top = gEve->GetCurrentEvent();

   gMultiView->DestroyEventRPhi();
   gMultiView->ImportEventRPhi(top);

   gMultiView->DestroyEventRhoZ();
   gMultiView->ImportEventRhoZ(top);

   update_html_summary();

   gEve->Redraw3D(kFALSE, kTRUE);
}


/******************************************************************************/
// GUI
/******************************************************************************/

//______________________________________________________________________________
//
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{
public:
   void Fwd()
   {
      if (esd_event_id < esd_tree->GetEntries() - 1) {
         ++esd_event_id;
         load_event();
      } else {
         printf("Already at last event.\n");
      }
   }
   void Bck()
   {
      if (esd_event_id > 0) {
         --esd_event_id;
         load_event();
      } else {
         printf("Already at first event.\n");
      }
   }
};

//______________________________________________________________________________
void make_gui()
{
   // Create minimal GUI for event navigation.

   TEveBrowser* browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("XX GUI");
   frmMain->SetCleanup(kDeepCleanup);

   TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
   {

      TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
      TGPictureButton* b = 0;
      EvNavHandler    *fh = new EvNavHandler;

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
   }
   frmMain->AddFrame(hf);

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);
}


/******************************************************************************/
// Code for reading AliESD and creating visualization objects
/******************************************************************************/

enum ESDTrackFlags {
   kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
   kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
   kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
   kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
   kHMPIDpid=0x20000,
   kEMCALmatch=0x40000,
   kTRDbackup=0x80000,
   kTRDStop=0x20000000,
   kESDpid=0x40000000,
   kTIME=0x80000000
};

//______________________________________________________________________________
void alice_esd_read()
{
   // Read tracks and associated clusters from current event.

   AliESDRun    *esdrun = (AliESDRun*)    esd_objs->FindObject("AliESDRun");
   TClonesArray *tracks = (TClonesArray*) esd_objs->FindObject("Tracks");

   // This needs further investigation. Clusters not shown.
   // esd_friend = (AliESDfriend*) esd_objs->FindObject("AliESDfriend");
   // printf("Friend %p, n_tracks:%d\n",
   //        esd_friend,
   //        esd_friend->fTracks.GetEntries());

   if (gTrackList == 0)
   {
      gTrackList = new TEveTrackList("ESD Tracks");
      gTrackList->SetMainColor(6);
      gTrackList->SetMarkerColor(kYellow);
      gTrackList->SetMarkerStyle(4);
      gTrackList->SetMarkerSize(0.5);

      gEve->AddElement(gTrackList);
   }

   TEveTrackPropagator* trkProp = gTrackList->GetPropagator();
   trkProp->SetMagField( 0.1 * esdrun->fMagneticField ); // kGaus to Tesla

   for (Int_t n=0; n<tracks->GetEntriesFast(); ++n)
   {
      AliESDtrack* at = (AliESDtrack*) tracks->At(n);

      // If ITS refit failed, take track parameters at inner TPC radius.
      AliExternalTrackParam* tp = at;
      if (! trackIsOn(at, kITSrefit)) {
         tp = at->fIp;
      }

      TEveTrack* track = esd_make_track(trkProp, n, at, tp);
      track->SetAttLineAttMarker(gTrackList);
      gTrackList->AddElement(track);

      // This needs further investigation. Clusters not shown.
      // if (frnd)
      // {
      //     AliESDfriendTrack* ft = (AliESDfriendTrack*) frnd->fTracks->At(n);
      //     printf("%d friend = %p\n", ft);
      // }
   }

   gTrackList->MakeTracks();
}

//______________________________________________________________________________
TEveTrack* esd_make_track(TEveTrackPropagator*   trkProp,
                          Int_t                  index,
                          AliESDtrack*           at,
                          AliExternalTrackParam* tp)
{
   // Helper function creating TEveTrack from AliESDtrack.
   //
   // Optionally specific track-parameters (e.g. at TPC entry point)
   // can be specified via the tp argument.

   Double_t      pbuf[3], vbuf[3];
   TEveRecTrack  rt;

   if (tp == 0) tp = at;

   rt.fLabel  = at->fLabel;
   rt.fIndex  = index;
   rt.fStatus = (Int_t) at->fFlags;
   rt.fSign   = (tp->fP[4] > 0) ? 1 : -1;

   trackGetPos(tp, vbuf);      rt.fV.Set(vbuf);
   trackGetMomentum(tp, pbuf); rt.fP.Set(pbuf);

   Double_t ep = trackGetP(at);
   Double_t mc = 0.138; // at->GetMass(); - Complicated function, requiring PID.

   rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

   TEveTrack* track = new TEveTrack(&rt, trkProp);
   track->SetName(Form("TEveTrack %d", rt.fIndex));
   track->SetStdTitle();

   return track;
}

//______________________________________________________________________________
Bool_t trackIsOn(AliESDtrack* t, Int_t mask)
{
   // Check is track-flag specified by mask are set.

   return (t->fFlags & mask) > 0;
}

//______________________________________________________________________________
void trackGetPos(AliExternalTrackParam* tp, Double_t r[3])
{
   // Get global position of starting point of tp.

  r[0] = tp->fX; r[1] = tp->fP[0]; r[2] = tp->fP[1];

  Double_t cs=TMath::Cos(tp->fAlpha), sn=TMath::Sin(tp->fAlpha), x=r[0];
  r[0] = x*cs - r[1]*sn; r[1] = x*sn + r[1]*cs;
}

//______________________________________________________________________________
void trackGetMomentum(AliExternalTrackParam* tp, Double_t p[3])
{
   // Return global momentum vector of starting point of tp.

   p[0] = tp->fP[4]; p[1] = tp->fP[2]; p[2] = tp->fP[3];

   Double_t pt=1./TMath::Abs(p[0]);
   Double_t cs=TMath::Cos(tp->fAlpha), sn=TMath::Sin(tp->fAlpha);
   Double_t r=TMath::Sqrt(1 - p[1]*p[1]);
   p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];
}

//______________________________________________________________________________
Double_t trackGetP(AliExternalTrackParam* tp)
{
   // Return magnitude of momentum of tp.

   return TMath::Sqrt(1.+ tp->fP[3]*tp->fP[3])/TMath::Abs(tp->fP[4]);
}

#endif

