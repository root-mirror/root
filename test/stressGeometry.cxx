// Program to check a TGeo geometry
// The first time you run this program, the geometry files will be taken
// from http://root.cern/files
//
//    How the program works
// If the file <geom_name>_ref.root does not exist, it is generated. The file
// contains a TTree with Npoints (default=100000) obtained with the following
// algorithm:
//   -a point is generated with a uniform distribution x,y,z in the master volume
//   -a direction theta, phi is generated uniformly in -2pi<phi<2pi and 0<theta<pi
//   -gGeoManager finds the geometry path for the point
//   -the number of boundaries (nbound), total length (length), safety distance
// from the starting point (safe) and number of radiation lengths (rad) from x,y,z
// is calculated to the exit of the detector. The total number of crossings, detector
// weight and total number of radiation lengths for all tracks are stored as user info in the tree.
//
//  Using the file <geom_name>_ref.root (generated typically with a previous version
//  of the TGeo classes), the Npoints in the Tree are used to perform the
//  same operation with the new version.
//  In case of a disagreement, an error message is reported.
//
//  The ReadRef case is also used as a benchmark
//  The ROOTMARKS reported are relative to a Linux/P IV 2.8 GHz gcc3.2.3 machine
//  normalized at 800 ROOTMARKS when running with CINT.
//
// To run this script, do
//   stressGeometry
// or  stressGeometry *
// or  stressGeometry alice
// or from the ROOT command line
// root > .L stressGeometry.cxx  or .L stressGeometry.cxx+
// root > stressGeometry(exp_name); // where exp_name is the geometry file name without .root
// OR simply: stressGeometry(); to run tests for a set of geometries
//
// Authors: Rene Brun, Andrei Gheata, 22 march 2005

#include "TStopwatch.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoBBox.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TError.h"
#include "TApplication.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVirtualGeoConverter.h"

// Total and reference times
Double_t tpstot = 0;
Double_t tpsref = 112.1; //time including the generation of the ref files
Bool_t testfailed = kFALSE;
#ifndef __CINT__
void stressGeometry(const char*, Bool_t, Bool_t);

int main(int argc, char **argv)
{
   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   Bool_t vecgeom = kFALSE;
   TString geom = "*";
   if (argc > 1) geom = argv[1];
   geom.ToLower();
   if (geom == "all") geom = "*";
   printf("geom: %s\n", geom.Data());

   if (argc > 1) {
       for (Int_t iarg=1; iarg<argc; ++iarg) {
          if (!strcmp(argv[iarg], "vecgeom")) vecgeom = kTRUE;
       }
   }
   stressGeometry(geom,kFALSE,vecgeom);
   return 0;
}

#endif

// data structure for one point
typedef struct {
   Double_t x,y,z,theta,phi;      // Initial track position and direction
   Int_t    nbound;               // Number of boundaries crossed until exit
   Float_t  length;               // Total length up to exit
   Float_t  safe;                 // Safety distance for the initial location
   Float_t  rad;                  // Number of radiation lengths up to exit
} p_t;
p_t p;

const Int_t NG = 33;
const char *exps[NG] = {"aleph",
                        "barres",
                        "felix",
                        "phenix",
                        "chambers",
                        "p326",
                        "bes",
                        "dubna",
                        "ganil",
                        "e907",
                        "phobos2",
                        "hermes",
                        "na35",
                        "na47",
                        "na49",
                        "wa91",
                        "sdc",
                        "integral",
                        "ams",
                        "brahms",
                        "gem",
                        "tesla",
                        "btev",
                        "cdf",
                        "hades2",
                        "lhcbfull",
                        "star",
                        "sld",
                        "cms",
                        "alice3",
                        "babar2",
                        "belle",
                        "atlas"
};
const Int_t versions[NG] =  {5, //aleph
                             3, //barres
                             3, //felix
                             3, //phenix
                             3, //chambers
                             4, //p326
                             3, //bes
                             3, //dubna
                             3, //ganil
                             3, //e907
                             4, //phobos2
                             3, //hermes
                             3, //na35
                             3, //na47
                             3, //na49
                             3, //wa91
                             3, //sdc
                             4, //integral
                             4, //ams
                             3, //brahms
                             5, //gem
                             4, //tesla
                             3, //btev
                             6, //cdf
                             4, //hades2
                             4, //lhcbfull
                             4, //star
                             4, //sld
                             4, //cms
                             6, //alice3
                             4, //babar2
                             3, //belle
                             6}; //atlas
// The timings below are on my machine PIV 3GHz
const Double_t cp_brun[NG] = {1.9,  //aleph
                              0.1,  //barres
                              0.12, //felix
                              0.62, //phenix
                              0.1,  //chambers
                              0.19, //p326
                              1.2,  //bes
                              0.12, //dubna
                              0.11, //ganil
                              0.17, //e907
                              0.22, //phobos2
                              0.24, //hermes
                              0.14, //na35
                              0.21, //na47
                              0.23, //na49
                              0.16, //wa91
                              0.17, //sdc
                              0.63, //integral
                              0.9,  //ams
                              1.1,  //brahms
                              1.8,  //gem
                              1.5,  //tesla
                              1.6,  //btev
                              2.2,  //cdf
                              1.2,  //hades2
                              1.6,  //lhcbfull
                              2.7,  //star
                              3.3,  //sld
                              7.5,  //cms
                              8.0,  //alice2
                             19.6,  //babar2
                             24.1,  //belle
                             26.7}; //atlas
// Bounding boxes for experiments
Double_t boxes[NG][3] = {{600,600,500},     // aleph
                         {100,100,220},     // barres
                         {200,200,12000},   // felix
                         {750,750,1000},    // phenix
                         {500,500,500},     // chambers
                         {201,201,26000},   // p326
                         {400,400,240},     // bes
                         {500,500,2000},    // dubna
                         {500,500,500},     // ganil
                         {250,250,2000},    // e907
                         {400,40,520},      // phobos2
                         {250,250,770},     // hermes
                         {310,160,1500},    // na35
                         {750,500,3000},    // na47
                         {600,200,2000},    // na49
                         {175,325,680},     // wa91
                         {1400,1400,2100},  // sdc
                         {100,100,200},     // integral
                         {200,200,200},     // ams
                         {50,50,50},        // brahms
                         {2000,2000,5000},  // gem
                         {1500,1500,1500},  // tesla
                         {600,475,1270},    // btev
                         {500,500,500},     // cdf
                         {250,250,200},     // hades2
                         {6700,5000,19000}, // lhcbfull
                         {350,350,350},     // star
                         {500,500,500},     // sld
                         {800,800,1000},    // cms
                         {400,400,400},     // alice2
                         {300,300,400},     // babar2
                         {440,440,538},     // belle
                         {1000,1000,1500}   // atlas
};

Int_t iexp[NG];
Bool_t gen_ref=kFALSE;
void FindRad(Double_t x, Double_t y, Double_t z,Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose=kFALSE);
void ReadRef(Int_t kexp);
void WriteRef(Int_t kexp);
void InspectRef(const char *exp="alice", Int_t vers=3);

void stressGeometry(const char *exp="*", Bool_t generate_ref=kFALSE, Bool_t vecgeom=kFALSE) {
   TGeoManager::SetVerboseLevel(0);
   gen_ref = generate_ref;
   gErrorIgnoreLevel = 10;

   fprintf(stderr,"******************************************************************\n");
   fprintf(stderr,"* STRESS GEOMETRY\n");
   TString opt = exp;
   opt.ToLower();
   Bool_t all = kFALSE;
   if (opt.Contains("*")) all = kTRUE;
   Int_t i;
   for (i=0; i<NG; i++) {
      if (all) {
         iexp[i] = 1;
         continue;
      }
      if (opt.Contains(exps[i])) iexp[i] = 1;
      else                       iexp[i] = 0;
   }
#if defined(linux) && !defined(__x86_64__)
   // 32bit linux: we have an error with ATLAS, see https://sft.its.cern.ch/jira/browse/ROOT-9893
   // Disable unless explicitly enabled.
   if (all) {
      printf("DISABLED ATLAS TEST due to known failure on Linux 32 bit!\n");
      iexp[32] = 0;
   }
#endif
   TFile::SetCacheFileDir(".");
   TString fname;
   for (i=0; i<NG; i++) {
      if (!iexp[i]) continue;
      fname = TString::Format("%s.root", exps[i]);
      if (gGeoManager) {
         delete gGeoManager;
         gGeoManager = 0;
      }
      TGeoManager::Import(Form("http://root.cern/files/%s",fname.Data()));
      if (!gGeoManager) return;
      if (vecgeom) TVirtualGeoConverter::Instance()->ConvertGeometry();
      
      fname = TString::Format("files/%s_ref_%d.root", exps[i],versions[i]);

      if (gen_ref || !TFile::Open(Form("http://root.cern/files/%s_ref_%d.root",exps[i],versions[i]),"CACHEREAD")) {
         if (!gen_ref) fprintf(stderr,"File: %s does not exist, generating it\n", fname.Data());
         else               fprintf(stderr,"Generating reference file %s\n", fname.Data());
         WriteRef(i);
      }

      ReadRef(i);
   }
   if (all && tpstot>0) {
      Float_t rootmarks = 800*tpsref/tpstot;
      Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
      if (UNIX) {
         TString sp = gSystem->GetFromPipe("uname -a");
         sp.Resize(60);
         printf("*  SYS: %s\n",sp.Data());
         if (strstr(gSystem->GetBuildNode(),"Linux")) {
            sp = gSystem->GetFromPipe("lsb_release -d -s");
            printf("*  SYS: %s\n",sp.Data());
         }
         if (strstr(gSystem->GetBuildNode(),"Darwin")) {
            sp  = gSystem->GetFromPipe("sw_vers -productVersion");
            sp += " Mac OS X ";
            printf("*  SYS: %s\n",sp.Data());
         }
      } else {
         const char *os = gSystem->Getenv("OS");
         if (!os) fprintf(stderr,"*  SYS: Windows 95\n");
         else     fprintf(stderr,"*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
      }
      fprintf(stderr,"******************************************************************\n");
      if (testfailed) fprintf(stderr,"*  stressGeometry found bad points ............. FAILED\n");
      else          fprintf(stderr,"*  stressGeometry .................................. OK\n");
      fprintf(stderr,"******************************************************************\n");
      fprintf(stderr,"*  CPU time in ReadRef = %6.2f seconds\n",tpstot);
      fprintf(stderr,"*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime());
   }
   fprintf(stderr,"******************************************************************\n");
}

void ReadRef(Int_t kexp) {
   TString fname;
   TFile *f = 0;
   //use ref_[version[i]] files
   if (!gen_ref)
      fname = TString::Format("http://root.cern/files/%s_ref_%d.root", exps[kexp],versions[kexp]);
   else
      fname.Format("files/%s_ref_%d.root", exps[kexp],versions[kexp]);

   f = TFile::Open(fname,"CACHEREAD");
   if (!f) {
      fprintf(stderr,"Reference file %s not found ! Skipping.\n", fname.Data());
      return;
   }
   // fprintf(stderr,"Reference file %s found\n", fname.Data());
   fname = TString::Format("%s_diff.root", exps[kexp]);
   TFile fdiff(fname,"RECREATE");
   TTree *TD = new TTree("TD","TGeo stress diff");
   TD->Branch("p",&p.x,"x/D:y/D:z/D:theta/D:phi/D:rad[4]/F");
   TTree *T = (TTree*)f->Get("T");
   T->SetBranchAddress("p",&p.x);
   Long64_t nentries = T->GetEntries();
   TVectorD *vref = (TVectorD *)T->GetUserInfo()->At(0);
   if (!vref) {
      fprintf(stderr," ERROR: User info not found, regenerate reference file\n");
      return;
   }
   TVectorD vect(4);
   TVectorD vect_ref = *vref;
   Int_t nbound;
   Float_t length, safe, rad;
   Float_t diff;
   Float_t diffmax = 0.01;  // percent of rad!
   Int_t nbad = 0;
   vect(0) = 0;//gGeoManager->Weight(0.01, "va");
   TStopwatch sw;
   for (Long64_t i=0;i<nentries;i++) {
      T->GetEntry(i);
      nbound = 0;
      length = 0.;
      safe = 0.;
      rad = 0.;
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound, length, safe, rad);
      vect(1) += Double_t(nbound);
      vect(2) += length;
      vect(3) += rad;
      diff = 0;
      diff += TMath::Abs(length-p.length);
      diff += TMath::Abs(safe-p.safe);
      diff += TMath::Abs(rad-p.rad);
      if (((p.rad>0) && (TMath::Abs(rad-p.rad)/p.rad)>diffmax) ||
           TMath::Abs(nbound-p.nbound)>100) {
         nbad++;
         if (nbad < 10) {
            fprintf(stderr," ==>Point %lld differs with diff = %g, x=%g, y=%g, z=%g\n",i,diff,p.x,p.y,p.z);
            fprintf(stderr,"    p.nbound=%d, p.length=%g, p.safe=%g, p.rad=%g\n",
                        p.nbound,p.length,p.safe,p.rad);
            fprintf(stderr,"      nbound=%d,   length=%g,   safe=%g,   rad=%g\n",
                        nbound,length,safe,rad);
         }
         TD->Fill();
         p.nbound = nbound;
         p.length = length;
         p.safe   = safe;
         p.rad    = rad;
         TD->Fill();
      }
   }
   diff = 0.;
   //for (Int_t j=1; j<4; j++) diff += TMath::Abs(vect_ref(j)-vect(j));
   diff += TMath::Abs(vect_ref(3)-vect(3))/vect_ref(3);
   if (diff > diffmax) {
//      fprintf(stderr,"Total weight=%g   ref=%g\n", vect(0), vect_ref(0));
      fprintf(stderr,"Total nbound=%g   ref=%g\n", vect(1), vect_ref(1));
      fprintf(stderr,"Total length=%g   ref=%g\n", vect(2), vect_ref(2));
      fprintf(stderr,"Total    rad=%g   ref=%g\n", vect(3), vect_ref(3));
      nbad++;
   }

   if (nbad) {
      testfailed = kTRUE;
      TD->AutoSave();
      TD->Print();
   }
   delete TD;
   delete f;

   Double_t cp = sw.CpuTime();
   tpstot += cp;
   if (nbad > 0) fprintf(stderr,"*     stress %-15s  found %5d bad points ............. failed\n",exps[kexp],nbad);
   else          fprintf(stderr,"*     stress %-15s: time/ref = %6.2f/%6.2f............ OK\n",exps[kexp],cp,cp_brun[kexp]);
}

void WriteRef(Int_t kexp) {
   TRandom3 r;
//   Double_t theta, phi;
   Double_t point[3];
   TVectorD vect(4);
   TGeoShape *top = gGeoManager->GetMasterVolume()->GetShape();
//   TGeoBBox *box = (TGeoBBox*)top;
   Double_t xmax = boxes[kexp][0]; //box->GetDX(); // 300;
   Double_t ymax = boxes[kexp][1]; //box->GetDY(); // 300;
   Double_t zmax = boxes[kexp][2]; //box->GetDZ(); // 500;
   TString fname(TString::Format("files/%s_ref_%d.root", exps[kexp], versions[kexp]));
   TFile f(fname,"recreate");
   TTree *T = new TTree("T","TGeo stress");
   T->Branch("p",&p.x,"x/D:y/D:z/D:theta/D:phi/D:nbound/I:length/F:safe/F:rad/F");
   T->GetUserInfo()->Add(&vect);
   Long64_t Npoints = 10000;
   Long64_t i = 0;
   vect(0) = 0; //gGeoManager->Weight(0.01, "va");
   while (i<Npoints) {
      p.x  = r.Uniform(-xmax,xmax);
      p.y  = r.Uniform(-ymax,ymax);
      p.z  = r.Uniform(-zmax,zmax);
      point[0] = p.x;
      point[1] = p.y;
      point[2] = p.z;
      if (top->Contains(point)) {
         p.phi   =  2*TMath::Pi()*r.Rndm();
         p.theta = TMath::ACos(1.-2.*r.Rndm());
         FindRad(p.x,p.y,p.z, p.theta, p.phi, p.nbound, p.length, p.safe, p.rad);
         vect(1) += Double_t(p.nbound);
         vect(2) += p.length;
         vect(3) += p.rad;
         T->Fill();
         i++;
      }
   }
   T->AutoSave();
   T->GetUserInfo()->Remove(&vect);
//   T->Print();
   delete T;
}

void FindRad(Double_t x, Double_t y, Double_t z,Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose) {
   Double_t xp  = TMath::Sin(theta)*TMath::Cos(phi);
   Double_t yp  = TMath::Sin(theta)*TMath::Sin(phi);
   Double_t zp  = TMath::Cos(theta);
   Double_t snext;
   TString path;
   Double_t pt[3];
   Double_t loc[3];
   Double_t epsil = 1.E-2;
   Double_t lastrad = 0.;
   Int_t ismall = 0;
   nbound = 0;
   length = 0.;
   safe   = 0.;
   rad    = 0.;
   TGeoMedium *med;
   TGeoShape *shape;
   TGeoNode *lastnode;
   gGeoManager->InitTrack(x,y,z,xp,yp,zp);
   if (verbose) {
      fprintf(stderr,"Track: (%15.10f,%15.10f,%15.10f,%15.10f,%15.10f,%15.10f)\n",
                       x,y,z,xp,yp,zp);
      path = gGeoManager->GetPath();
   }
   TGeoNode *nextnode = gGeoManager->GetCurrentNode();
   safe = gGeoManager->Safety();
   while (nextnode) {
      med = nextnode->GetVolume()->GetMedium();
      shape = nextnode->GetVolume()->GetShape();
      lastnode = nextnode;
      nextnode = gGeoManager->FindNextBoundaryAndStep();
      snext  = gGeoManager->GetStep();
      if (snext<1.e-8) {
         ismall++;
         if ((ismall<3) && (lastnode != nextnode)) {
            // First try to cross a very thin layer
            length += snext;
            nextnode = gGeoManager->FindNextBoundaryAndStep();
            snext  = gGeoManager->GetStep();
            if (snext<1.E-8) continue;
            // We managed to cross the layer
            ismall = 0;
         } else {
            // Relocate point
            if (ismall > 3) {
               fprintf(stderr,"ERROR: Small steps in: %s shape=%s\n",gGeoManager->GetPath(), shape->ClassName());
               return;
            }
            memcpy(pt,gGeoManager->GetCurrentPoint(),3*sizeof(Double_t));
            const Double_t *dir = gGeoManager->GetCurrentDirection();
            for (Int_t i=0;i<3;i++) pt[i] += epsil*dir[i];
            snext = epsil;
            length += snext;
            rad += lastrad*snext;
            gGeoManager->CdTop();
            nextnode = gGeoManager->FindNode(pt[0],pt[1],pt[2]);
            if (gGeoManager->IsOutside()) return;
            TGeoMatrix *mat = gGeoManager->GetCurrentMatrix();
            mat->MasterToLocal(pt,loc);
            if (!gGeoManager->GetCurrentVolume()->Contains(loc)) {
//            fprintf(stderr,"Woops - out\n");
               gGeoManager->CdUp();
               nextnode = gGeoManager->GetCurrentNode();
            }
            continue;
         }
      } else {
         ismall = 0;
      }
      nbound++;
      length += snext;
      if (med) {
         Double_t radlen = med->GetMaterial()->GetRadLen();
         if (radlen>1.e-5 && radlen<1.e10) {
            lastrad = med->GetMaterial()->GetDensity()/radlen;
            rad += lastrad*snext;
         } else {
            lastrad = 0.;
         }
         if (verbose) {
            fprintf(stderr," STEP #%d: %s\n",nbound, path.Data());
            fprintf(stderr,"    step=%g  length=%g  rad=%g %s\n", snext,length,
                   med->GetMaterial()->GetDensity()*snext/med->GetMaterial()->GetRadLen(),med->GetName());
            path =  gGeoManager->GetPath();
         }
      }
   }
}

void InspectDiff(const char* exp="alice",Long64_t ientry=-1) {
   Int_t nbound = 0;
   Float_t length = 0.;
   Float_t safe   = 0.;
   Float_t rad    = 0.;
   TString fname(TString::Format("%s.root",exp));
   if (gSystem->AccessPathName(fname)) {
      TGeoManager::Import(Form("http://root.cern/files/%s",fname.Data()));
   } else {
      TGeoManager::Import(fname);
   }
   fname = TString::Format("%s_diff.root",exp);
   TFile f(fname);
   if (f.IsZombie()) return;
   TTree *TD = (TTree*)f.Get("TD");
   TD->SetBranchAddress("p",&p.x);
   Long64_t nentries = TD->GetEntries();
   nentries = nentries>>1;
   if (ientry>=0 && ientry<nentries) {
      fprintf(stderr,"DIFFERENCE #%lld\n", ientry);
      TD->GetEntry(2*ientry);
      fprintf(stderr,"   NEW: nbound=%d  length=%g  safe=%g  rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      TD->GetEntry(2*ientry+1);
      fprintf(stderr,"   OLD: nbound=%d  length=%g  safe=%g  rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound,length,safe,rad, kTRUE);
      return;
   }
   for (Long64_t i=0;i<nentries;i++) {
      fprintf(stderr,"DIFFERENCE #%lld\n", i);
      TD->GetEntry(2*i);
      fprintf(stderr,"   NEW: nbound=%d  length=%g  safe=%g rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      TD->GetEntry(2*i+1);
      fprintf(stderr,"   OLD: nbound=%d  length=%g  safe=%g rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound,length,safe,rad, kTRUE);
   }
}

void InspectRef(const char *exp, Int_t vers) {
// Inspect current reference.
   TString fname(TString::Format("%s_ref_%d.root", exp, vers));
   if (gSystem->AccessPathName(fname)) {
      fprintf(stderr,"ERROR: file %s does not exist\n", fname.Data());
      return;
   }
   TFile f(fname);
   if (f.IsZombie()) return;
   TTree *T = (TTree*)f.Get("T");
   Long64_t nentries = T->GetEntries();
   fname.Format("Stress test for %s geometry", exp);
   TCanvas *c = new TCanvas("stress", fname,700,800);
   c->Divide(2,2,0.005,0.005);
   c->cd(1);
   gPad->SetLogy();
   T->Draw("p.nbound","","", nentries, 0);
   c->cd(2);
   gPad->SetLogy();
   T->Draw("p.length","","", nentries, 0);
   c->cd(3);
   gPad->SetLogy();
   T->Draw("p.safe","","", nentries, 0);
   c->cd(4);
   gPad->SetLogy();
   T->Draw("p.rad","","", nentries, 0);
   c->cd(0);
   c->SetFillColor(kYellow);
   TVectorD *vref = (TVectorD *)T->GetUserInfo()->At(0);
   TVectorD vect = *vref;
   fprintf(stderr,"=====================================\n");
//   fprintf(stderr,"Total weight:  %g [kg]\n", vect(0));
   fprintf(stderr,"Total nbound:  %g boundaries crossed\n", vect(1));
   fprintf(stderr,"Total length:  %g [m]\n", 0.01*vect(2));
   fprintf(stderr,"Total nradlen: %f\n", vect(3));
   fprintf(stderr,"=====================================\n");
}
