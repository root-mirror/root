// @(#)root/eg:$Id$
// Author: Ola Nordmann   21/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class  TGenerator
 \ingroup eg

The interface to various event generators

Is an base class, that defines the interface of ROOT to various      
event generators. Every event generator should inherit from          
TGenerator or its subclasses.                                        
                                                                     
Derived class can overload the member  function GenerateEvent        
to do the actual event generation (e.g., call PYEVNT or similar).    
                                                                     
The derived class should overload the member function                
ImportParticles (both types) to read the internal storage of the     
generated event into either the internal TObjArray or the passed     
TClonesArray of TParticles.                                          
                                                                     
If the generator code stores event data in the /HEPEVT/ common block 
Then the default implementation of ImportParticles should suffice.   
The common block /HEPEVT/ is structed like                           
                                                                     
\verbatim
  // C                                                            
  typedef struct {                                                   
     Int_t    nevhep;           // Event number                      
     Int_t    nhep;             // # of particles                    
     Int_t    isthep[4000];     // Status flag of i'th particle      
     Int_t    idhep[4000];      // PDG # of particle                 
     Int_t    jmohep[4000][2];  // 1st & 2nd mother particle #       
     Int_t    jdahep[4000][2];  // 1st & 2nd daughter particle #     
     Double_t phep[4000][5];    // 4-momentum and 1 word             
     Double_t vhep[4000][4];    // 4-position of production          
  } HEPEVT_DEF;                                                      
                                                                     
                                                                     
  C Fortran                                                          
        COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(4000),IDHEP(4000),          
      +    JMOHEP(2,4000),JDAHEP(2,4000),PHEP(5,4000),VHEP(4,4000)   
        INTEGER NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP               
        DOUBLE PRECISION PHEP,VHEP                                   
\endverbatim

The generic member functions SetParameter and GetParameter can be
overloaded to set and get parameters of the event generator.

Note, if the derived class interfaces a (set of) Fortran common
blocks (like TPythia, TVenus does), one better make the derived
class a singleton.  That is, something like
                         
\verbatim
    class MyGenerator : public TGenerator                            
    {                                                                
    public:                                                          
      static MyGenerator* Instance()                                 
      {                                                              
        if (!fgInstance) fgInstance = new MyGenerator;               
        return fgInstance;                                           
      }                                                              
      void  GenerateEvent() { ... }                                  
      void  ImportParticles(TClonesArray* a, Option_t opt="") {...}  
      Int_t ImportParticles(Option_t opt="") { ... }                 
      Int_t    SetParameter(const char* name, Double_t val) { ... }  
      Double_t GetParameter(const char* name) { ... }                
      virtual ~MyGenerator() { ... }                                 
    protected:                                                       
      MyGenerator() { ... }                                          
      MyGenerator(const MyGenerator& o) { ... }                      
      MyGenerator& operator=(const MyGenerator& o) { ... }           
      static MyGenerator* fgInstance;                                
      ClassDef(MyGenerator,0);                                       
    };                                                               
\endverbatim
                                                                     
Having multiple objects accessing the same common blocks is not      
safe.                                                                
                                                                     
Concrete TGenerator classes can be loaded in scripts and subseqent-  
ly used in compiled code:                                            
                                                                     
\verbatim
    // MyRun.h                                                       
    class MyRun : public TObject                                     
    {                                                                
    public:                                                          
      static MyRun* Instance() { ... }                               
      void SetGenerator(TGenerator* g) { fGenerator = g; }           
      void Run(Int_t n, Option_t* option="")                         
      {                                                              
        TFile*        file = TFile::Open("file.root","RECREATE");    
        TTree*        tree = new TTree("T","T");                     
        TClonesArray* p    = new TClonesArray("TParticles");         
        tree->Branch("particles", &p);                               
        for (Int_t event = 0; event < n; event++) {                  
          fGenerator->GenerateEvent();                               
          fGenerator->ImportParticles(p,option);                     
          tree->Fill();                                              
        }                                                            
        file->Write();                                               
        file->Close();                                               
      }                                                              
      ...                                                            
    protected:                                                       
      TGenerator* fGenerator;                                        
      ClassDef(MyRun,0);                                             
    };                                                               
                                                                     
    // Config.C                                                      
    void Config()                                                    
    {                                                                
       MyRun* run = MyRun::Instance();                               
       run->SetGenerator(MyGenerator::Instance());                   
    }                                                                
                                                                     
    // main.cxx                                                      
    int                                                              
    main(int argc, char** argv)                                      
    {                                                                
      TApplication app("", 0, 0);                                    
      gSystem->ProcessLine(".x Config.C");                           
      MyRun::Instance()->Run(10);                                    
      return 0;                                                      
    }                                                                
\endverbatim
                                                                     
This is especially useful for example with TVirtualMC or similar.
*/

#include "TROOT.h"
#include "TGenerator.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "Hepevt.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TText.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <iostream>


ClassImp(TGenerator);

////////////////////////////////////////////////////////////////////////////////
///  Event generator default constructor
///

TGenerator::TGenerator(const char *name,const char *title): TNamed(name,title)
{
   //  Initialize particles table
   TDatabasePDG::Instance();
   //TDatabasePDG *pdg = TDatabasePDG::Instance();
   //if (!pdg->ParticleList()) pdg->Init();

   fPtCut        = 0;
   fShowNeutrons = kTRUE;
   fParticles    =  new TObjArray(10000);
}

////////////////////////////////////////////////////////////////////////////////
///  Event generator default destructor
///

TGenerator::~TGenerator()
{
   //do nothing
   if (fParticles) {
      fParticles->Delete();
      delete fParticles;
      fParticles = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// must be implemented in concrete class (see eg TPythia6)

void TGenerator::GenerateEvent()
{
}

////////////////////////////////////////////////////////////////////////////////
///
///  It reads the /HEPEVT/ common block which has been filled by the
///  GenerateEvent method. If the event generator does not use the
///  HEPEVT common block, This routine has to be overloaded by the
///  subclasses.
///
///  The default action is to store only the stable particles (ISTHEP =
///  1) This can be demanded explicitly by setting the option = "Final"
///  If the option = "All", all the particles are stored.
///

TObjArray* TGenerator::ImportParticles(Option_t *option)
{
   fParticles->Clear();
   Int_t numpart = HEPEVT.nhep;
   if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i<numpart; i++) {
         if (HEPEVT.isthep[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
            TParticle *p = new TParticle(
                                   HEPEVT.idhep[i],
                                   HEPEVT.isthep[i],
                                   HEPEVT.jmohep[i][0]-1,
                                   HEPEVT.jmohep[i][1]-1,
                                   HEPEVT.jdahep[i][0]-1,
                                   HEPEVT.jdahep[i][1]-1,
                                   HEPEVT.phep[i][0],
                                   HEPEVT.phep[i][1],
                                   HEPEVT.phep[i][2],
                                   HEPEVT.phep[i][3],
                                   HEPEVT.vhep[i][0],
                                   HEPEVT.vhep[i][1],
                                   HEPEVT.vhep[i][2],
                                   HEPEVT.vhep[i][3]);
            fParticles->Add(p);
         }
      }
   } else if (!strcmp(option,"All")) {
      for (Int_t i = 0; i<numpart; i++) {
         TParticle *p = new TParticle(
                                   HEPEVT.idhep[i],
                                   HEPEVT.isthep[i],
                                   HEPEVT.jmohep[i][0]-1,
                                   HEPEVT.jmohep[i][1]-1,
                                   HEPEVT.jdahep[i][0]-1,
                                   HEPEVT.jdahep[i][1]-1,
                                   HEPEVT.phep[i][0],
                                   HEPEVT.phep[i][1],
                                   HEPEVT.phep[i][2],
                                   HEPEVT.phep[i][3],
                                   HEPEVT.vhep[i][0],
                                   HEPEVT.vhep[i][1],
                                   HEPEVT.vhep[i][2],
                                   HEPEVT.vhep[i][3]);
         fParticles->Add(p);
      }
   }
   return fParticles;
}

////////////////////////////////////////////////////////////////////////////////
///
///  It reads the /HEPEVT/ common block which has been filled by the
///  GenerateEvent method. If the event generator does not use the
///  HEPEVT common block, This routine has to be overloaded by the
///  subclasses.
///
///  The function loops on the generated particles and store them in
///  the TClonesArray pointed by the argument particles.  The default
///  action is to store only the stable particles (ISTHEP = 1) This can
///  be demanded explicitly by setting the option = "Final" If the
///  option = "All", all the particles are stored.
///

Int_t TGenerator::ImportParticles(TClonesArray *particles, Option_t *option)
{
   if (particles == 0) return 0;
   TClonesArray &clonesParticles = *particles;
   clonesParticles.Clear();
   Int_t numpart = HEPEVT.nhep;
   if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i<numpart; i++) {
         if (HEPEVT.isthep[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
            new(clonesParticles[i]) TParticle(
                                   HEPEVT.idhep[i],
                                   HEPEVT.isthep[i],
                                   HEPEVT.jmohep[i][0]-1,
                                   HEPEVT.jmohep[i][1]-1,
                                   HEPEVT.jdahep[i][0]-1,
                                   HEPEVT.jdahep[i][1]-1,
                                   HEPEVT.phep[i][0],
                                   HEPEVT.phep[i][1],
                                   HEPEVT.phep[i][2],
                                   HEPEVT.phep[i][3],
                                   HEPEVT.vhep[i][0],
                                   HEPEVT.vhep[i][1],
                                   HEPEVT.vhep[i][2],
                                   HEPEVT.vhep[i][3]);
         }
      }
   } else if (!strcmp(option,"All")) {
      for (Int_t i = 0; i<numpart; i++) {
         new(clonesParticles[i]) TParticle(
                                   HEPEVT.idhep[i],
                                   HEPEVT.isthep[i],
                                   HEPEVT.jmohep[i][0]-1,
                                   HEPEVT.jmohep[i][1]-1,
                                   HEPEVT.jdahep[i][0]-1,
                                   HEPEVT.jdahep[i][1]-1,
                                   HEPEVT.phep[i][0],
                                   HEPEVT.phep[i][1],
                                   HEPEVT.phep[i][2],
                                   HEPEVT.phep[i][3],
                                   HEPEVT.vhep[i][0],
                                   HEPEVT.vhep[i][1],
                                   HEPEVT.vhep[i][2],
                                   HEPEVT.vhep[i][3]);
      }
   }
   return numpart;
}

////////////////////////////////////////////////////////////////////////////////
///browse generator

void TGenerator::Browse(TBrowser *)
{
   Draw();
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to objects in event
///

Int_t TGenerator::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t big = 9999;
   const Int_t inview = 0;
   Int_t dist = big;
   if (px > 50 && py > 50) dist = inview;
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
///
///  Insert one event in the pad list
///

void TGenerator::Draw(Option_t *option)
{
   // Create a default canvas if a canvas does not exist
   if (!gPad) {
      gROOT->MakeDefCanvas();
      if (gPad->GetVirtCanvas())
         gPad->GetVirtCanvas()->SetFillColor(13);
   }

   static Float_t rbox = 1000;
   Float_t rmin[3],rmax[3];
   TView *view = gPad->GetView();
   if (!strstr(option,"same")) {
      if (view) { view->GetRange(rmin,rmax); rbox = rmax[2];}
      gPad->Clear();
   }

   AppendPad(option);

   view = gPad->GetView();
   //    compute 3D view
   if (view) {
      view->GetRange(rmin,rmax);
      rbox = rmax[2];
   } else {
      view = TView::CreateView(1,0,0);
      if (view) view->SetRange(-rbox,-rbox,-rbox, rbox,rbox,rbox );
   }
   const Int_t kColorProton    = 4;
   const Int_t kColorNeutron   = 5;
   const Int_t kColorAntiProton= 3;
   const Int_t kColorPionPlus  = 6;
   const Int_t kColorPionMinus = 2;
   const Int_t kColorKaons     = 7;
   const Int_t kColorElectrons = 0;
   const Int_t kColorGamma     = 18;

   Int_t nProtons    = 0;
   Int_t nNeutrons   = 0;
   Int_t nAntiProtons= 0;
   Int_t nPionPlus   = 0;
   Int_t nPionMinus  = 0;
   Int_t nKaons      = 0;
   Int_t nElectrons  = 0;
   Int_t nGammas     = 0;

   Int_t ntracks = fParticles->GetEntriesFast();
   Int_t i,lwidth,color,lstyle;
   TParticlePDG *ap;
   TParticle *p;
   const char *name;
   Double_t etot,vx,vy,vz;
   Int_t ninvol = 0;
   for (i=0;i<ntracks;i++) {
      p = (TParticle*)fParticles->UncheckedAt(i);
      if(!p) continue;
      ap = (TParticlePDG*)p->GetPDG();
      vx = p->Vx();
      vy = p->Vy();
      vz = p->Vz();
      if  (vx*vx+vy*vy+vz*vz > rbox*rbox) continue;
      Float_t pt = p->Pt();
      if (pt < fPtCut) continue;
      etot = p->Energy();
      if (etot > 0.1) lwidth = Int_t(6*TMath::Log10(etot));
      else lwidth = 1;
      if (lwidth < 1) lwidth = 1;
      lstyle = 1;
      color = 0;
      name = ap->GetName();
      if (!strcmp(name,"n"))     { if (!fShowNeutrons) continue;
                                   color = kColorNeutron;    nNeutrons++;}
      if (!strcmp(name,"p"))     { color = kColorProton;     nProtons++;}
      if (!strcmp(name,"p bar")) { color = kColorAntiProton; nAntiProtons++;}
      if (!strcmp(name,"pi+"))   { color = kColorPionPlus;   nPionPlus++;}
      if (!strcmp(name,"pi-"))   { color = kColorPionMinus;  nPionMinus++;}
      if (!strcmp(name,"e+"))    { color = kColorElectrons;  nElectrons++;}
      if (!strcmp(name,"e-"))    { color = kColorElectrons;  nElectrons++;}
      if (!strcmp(name,"gamma")) { color = kColorGamma;      nGammas++; lstyle = 3; }
      if ( strstr(name,"K"))     { color = kColorKaons;      nKaons++;}
      p->SetLineColor(color);
      p->SetLineStyle(lstyle);
      p->SetLineWidth(lwidth);
      p->AppendPad();
      ninvol++;
   }

   // event title
   TPaveText *pt = new TPaveText(-0.94,0.85,-0.25,0.98,"br");
   pt->AddText((char*)GetName());
   pt->AddText((char*)GetTitle());
   pt->SetFillColor(42);
   pt->Draw();

   // Annotate color codes
   Int_t tcolor = 5;
   if (gPad->GetFillColor() == 10) tcolor = 4;
   TText *text = new TText(-0.95,-0.47,"Particles");
   text->SetTextAlign(12);
   text->SetTextSize(0.025);
   text->SetTextColor(tcolor);
   text->Draw();
   text->SetTextColor(kColorGamma);      text->DrawText(-0.95,-0.52,"(on screen)");
   text->SetTextColor(kColorGamma);      text->DrawText(-0.95,-0.57,"Gamma");
   text->SetTextColor(kColorProton);     text->DrawText(-0.95,-0.62,"Proton");
   text->SetTextColor(kColorNeutron);    text->DrawText(-0.95,-0.67,"Neutron");
   text->SetTextColor(kColorAntiProton); text->DrawText(-0.95,-0.72,"AntiProton");
   text->SetTextColor(kColorPionPlus);   text->DrawText(-0.95,-0.77,"Pion +");
   text->SetTextColor(kColorPionMinus);  text->DrawText(-0.95,-0.82,"Pion -");
   text->SetTextColor(kColorKaons);      text->DrawText(-0.95,-0.87,"Kaons");
   text->SetTextColor(kColorElectrons);  text->DrawText(-0.95,-0.92,"Electrons,etc.");

   text->SetTextColor(tcolor);
   text->SetTextAlign(32);
   char tcount[32];
   snprintf(tcount,12,"%d",ntracks);      text->DrawText(-0.55,-0.47,tcount);
   snprintf(tcount,12,"%d",ninvol);       text->DrawText(-0.55,-0.52,tcount);
   snprintf(tcount,12,"%d",nGammas);      text->DrawText(-0.55,-0.57,tcount);
   snprintf(tcount,12,"%d",nProtons);     text->DrawText(-0.55,-0.62,tcount);
   snprintf(tcount,12,"%d",nNeutrons);    text->DrawText(-0.55,-0.67,tcount);
   snprintf(tcount,12,"%d",nAntiProtons); text->DrawText(-0.55,-0.72,tcount);
   snprintf(tcount,12,"%d",nPionPlus);    text->DrawText(-0.55,-0.77,tcount);
   snprintf(tcount,12,"%d",nPionMinus);   text->DrawText(-0.55,-0.82,tcount);
   snprintf(tcount,12,"%d",nKaons);       text->DrawText(-0.55,-0.87,tcount);
   snprintf(tcount,12,"%d",nElectrons);   text->DrawText(-0.55,-0.92,tcount);

   text->SetTextAlign(12);
   if (nPionPlus+nPionMinus) {
      snprintf(tcount,31,"Protons/Pions= %4f",Float_t(nProtons)/Float_t(nPionPlus+nPionMinus));
   } else {
      strlcpy(tcount,"Protons/Pions= inf",31);
   }
   text->DrawText(-0.45,-0.92,tcount);

   if (nPionPlus+nPionMinus) {
      snprintf(tcount,31,"Kaons/Pions= %4f",Float_t(nKaons)/Float_t(nPionPlus+nPionMinus));
   } else {
      strlcpy(tcount,"Kaons/Pions= inf",31);
   }
   text->DrawText(0.30,-0.92,tcount);
}


////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event
///

void TGenerator::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (gPad->GetView()) {
      gPad->GetView()->ExecuteRotateView(event, px, py);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of particles in the stack

Int_t TGenerator::GetNumberOfParticles() const
{
   return fParticles->GetLast()+1;
}

////////////////////////////////////////////////////////////////////////////////
///  Returns pointer to primary number i;
///

TParticle *TGenerator::GetParticle(Int_t i) const
{
   if (!fParticles) return 0;
   Int_t n = fParticles->GetLast();
   if (i < 0 || i > n) return 0;
   return (TParticle*)fParticles->UncheckedAt(i);
}

////////////////////////////////////////////////////////////////////////////////
///
///  Paint one event
///

void TGenerator::Paint(Option_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
///
///  Set Pt threshold below which primaries are not drawn
///

void TGenerator::SetPtCut(Float_t ptcut)
{
   fPtCut = ptcut;
   Draw();
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
///
///  Set lower and upper values of the view range
///

void TGenerator::SetViewRadius(Float_t rbox)
{
   SetViewRange(-rbox,-rbox,-rbox,rbox,rbox,rbox);
}

////////////////////////////////////////////////////////////////////////////////
///
///  Set lower and upper values of the view range
///

void TGenerator::SetViewRange(Float_t xmin, Float_t ymin, Float_t zmin, Float_t xmax, Float_t ymax, Float_t zmax)
{
   TView *view = gPad->GetView();
   if (!view) return;
   view->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);

   Draw();
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
///
///  Set flag to display or not neutrons
///

void TGenerator::ShowNeutrons(Bool_t show)
{
   fShowNeutrons = show;
   Draw();
   gPad->Update();
}
