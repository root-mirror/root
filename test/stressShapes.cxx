// Author : Mihaela Gheata   12-01-03
#ifndef __CINT__
#include <TRandom3.h>
#include <TROOT.h>
#include <snprintf.h>
#include <TH1.h>
#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoMatrix.h>
#include <TBenchmark.h>
#include <TApplication.h>

void stressShapes();

int main(int argc, char **argv)
{
   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   stressShapes();
   return 0;
}

#endif
//--- This macro creates a simple geometry based on all shapes known
//--- by TGeo. The first test generates 1 million random points inside
//--- the bounding box of each shape and computes the volume of the
//--- shape as Vbbox*Ninside/Ntotal.
//--- The second test tracks 100K random rays in the geometry, histogramming
//--- the length of all segments passing through each different shape.
//--- It computes mean, RMS and sum of lengths of all segments inside a
//--- given shape and compares with reference values.
//
// This test program is automatically created by $ROOTSYS/test/Makefile.
// To run it in batch, execute stressGeom.
// To run this test with interactive CINT, do
// root > .x stressShapes.cxx++
// or
// root > .x stressShapes.cxx

void sample_volume(Int_t ivol)
{
   const Double_t vshape[16] = {40000.0, 36028.3, 39978.70, 48001.3, 28481.2,
      8726.2, 42345.4, 9808.2, 12566.8, 64655.6, 37730.4, 23579.7,
      25559.5, 18418.3, 49960.2, 47707.9};
   gRandom = new TRandom3();
   TGeoVolume *vol = (TGeoVolume*)gGeoManager->GetListOfVolumes()->At(ivol);
   TGeoShape *shape = vol->GetShape();
   Double_t dx = ((TGeoBBox*)shape)->GetDX();
   Double_t dy = ((TGeoBBox*)shape)->GetDY();
   Double_t dz = ((TGeoBBox*)shape)->GetDZ();
   Double_t ox = (((TGeoBBox*)shape)->GetOrigin())[0];
   Double_t oy = (((TGeoBBox*)shape)->GetOrigin())[1];
   Double_t oz = (((TGeoBBox*)shape)->GetOrigin())[2];
   Double_t ratio;
   Double_t point[3];
   Double_t ngen=1000000;
   Double_t iin=0;
   Double_t i;
   for (i=0; i<ngen; i++) {
      point[0] = ox-dx+2*dx*gRandom->Rndm();
      point[1] = oy-dy+2*dy*gRandom->Rndm();
      point[2] = oz-dz+2*dz*gRandom->Rndm();
      if (shape->Contains(point)) iin++;
   }
   ratio = Double_t(iin)/Double_t(ngen);
   Double_t vbox = 8*dx*dy*dz;
   Double_t vv = vbox*ratio;
   Double_t dvv = TMath::Abs(vv-vshape[ivol-1]);
   Double_t sigma = vv/TMath::Sqrt(iin+1);
   char result[16];
   snprintf(result,16, "FAILED");
   if (dvv<2*sigma) snprintf(result,16, "OK");
   printf("---> testing %-4s ............... %s\n", vol->GetName(), result);
}

void length()
{
   const Double_t rms[16] = {6.284, 10.79, 9.545, 14.15, 11.45,
      5.871, 7.673, 5.935, 7.61, 5.334, 6.581, 4.954,
      7.718, 3.238, 19.09, 14.85};
   const Double_t mean[16] = {19.34, 22.53, 18.87, 21.95, 23.29,
      16.73, 15.09, 9.516, 12.68, 8.852, 9.518, 7.432,
      8.881, 6.489, 28.29, 22.74};
   TObjArray *vlist = gGeoManager->GetListOfVolumes();
   TGeoVolume *volume;
   Int_t nvolumes = vlist->GetEntriesFast();
   Double_t len[17];
   TList *hlist = new TList();
   TH1F *hist;
   Int_t i=0;
   memset(len, 0, nvolumes*sizeof(Double_t));
   for (i=0; i<nvolumes; i++) {
      volume = (TGeoVolume*)(vlist->At(i));
      hist = new TH1F(volume->GetName(), "lengths inside", 100, 0, 100);
      hist->SetCanExtend(TH1::kAllAxes);
      hlist->Add(hist);
   }
   Int_t nrays = 100000;

   Double_t dir[3];
   TGeoNode *startnode, *endnode;
   Int_t istep=0, icrt;
   Int_t itot=0;
   Int_t n10=nrays/10;

   Double_t theta,phi, step;
   while (itot<nrays) {
      itot++;
      if (n10) {
         if ((itot%n10) == 0) printf("    %i percent\n", Int_t(100*itot/nrays));
      }
      phi = 2*TMath::Pi()*gRandom->Rndm();
      theta= TMath::ACos(1.-2.*gRandom->Rndm());
      dir[0]=TMath::Sin(theta)*TMath::Cos(phi);
      dir[1]=TMath::Sin(theta)*TMath::Sin(phi);
      dir[2]=TMath::Cos(theta);
      gGeoManager->InitTrack(0,0,0, dir[0], dir[1], dir[2]);
      startnode = gGeoManager->GetCurrentNode();
      if (gGeoManager->IsOutside()) startnode=0;
      icrt = 0;
      if (startnode) icrt =vlist->IndexOf(startnode->GetVolume());
      // find where we end-up
      gGeoManager->FindNextBoundary();
      step = gGeoManager->GetStep();
      endnode = gGeoManager->Step();
      while (step<1E10) {
         while (!gGeoManager->IsEntering()) {
            istep++;
            if (istep>10000) break;
            gGeoManager->SetStep(1E-3);
            endnode = gGeoManager->Step();
            step += 1E-3;
         }
         if (istep>10000) break;
           len[icrt] += step;
         hist = (TH1F*)(hlist->At(icrt));
         hist->Fill(step);
         // now see if we can make an other step
         if (endnode==0 && step>1E10) break;
         istep = 0;
         // generate an extra step to cross boundary
         startnode = endnode;
         icrt = 0;
         if (startnode) icrt =vlist->IndexOf(startnode->GetVolume());
         gGeoManager->FindNextBoundary();
         step = gGeoManager->GetStep();
         endnode = gGeoManager->Step();
      }
   }
   // draw all segments
   Double_t drms, dmean;
   for (i=1; i<nvolumes; i++) {
      volume = (TGeoVolume*)(vlist->At(i));
      hist = (TH1F*)(hlist->At(i));
      char result[16];
      drms = TMath::Abs(rms[i-1]-hist->GetRMS());
      dmean = TMath::Abs(mean[i-1]-hist->GetMean());
      snprintf(result,16, "FAILED");
      if (dmean<0.01) {
         if (drms<0.01) snprintf(result,16,"OK");
      }
      printf("   %-4s : mean_len=%7.4g RMS=%7.4g total_len=%11.4g ... %s\n",
             volume->GetName(), hist->GetMean(), hist->GetRMS(), len[i],result);
   }
   hlist->Delete();
   delete hlist;
}

void stressShapes()
{
// New geometry test suite. Creates a geometry containing all shape
// types. Loop over all volumes and compute the following :
//  - generate 1 million random points and count how many are inside
//    each shape -> compute volume of each shape
//  - generate 10000 random directions and propagate from the center
//    of each volume -> compute total step length to exit current shape

#ifdef __CINT__
   gSystem->Load("libGeom");
#endif

   gBenchmark = new TBenchmark();
   gBenchmark->Start("stressShapes");

   TGeoManager *geom = new TGeoManager("stressShapes", "arbitrary shapes");
   TGeoMaterial *mat;
   TGeoMixture *mix;
   //---> create some materials
   mat = new TGeoMaterial("Vacuum",0,0,0);
   mat->SetUniqueID(0);
   mat = new TGeoMaterial("Be", 9.01,4,1.848);
   mat->SetUniqueID(1);
   mat = new TGeoMaterial("Al", 26.98,13,2.7);
   mat->SetUniqueID(2);
   mat = new TGeoMaterial("Fe", 55.85,26,7.87);
   mat->SetUniqueID(3);
   mat = new TGeoMaterial("Cu", 63.55,29,8.96);
   mat->SetUniqueID(4);
   mat = new TGeoMaterial("C",12.01,6,2.265);
   mat->SetUniqueID(5);
   mat = new TGeoMaterial("Pb",207.19,82,11.35);
   mat->SetUniqueID(6);
   mat = new TGeoMaterial("Si",28.09,14,2.33);
   mat->SetUniqueID(7);
   mix = new TGeoMixture("scint",2,   1.03200    );
      mix->DefineElement(0,1.008,1,0.7749078E-01);
      mix->DefineElement(1,12,6,0.9225092);
   mix->SetUniqueID(8);
   mat = new TGeoMaterial("Li",6.94,3,0.534);
   mat->SetUniqueID(9);
   mat = new TGeoMaterial("N",14.01,7,0.808);
   mat->SetUniqueID(10);
   mat = new TGeoMaterial("Ne",20.18,10,1.207);
   mat->SetUniqueID(11);
   mat = new TGeoMaterial("Ts",183.85,74,19.3);
   mat->SetUniqueID(12);
   mat = new TGeoMaterial("CFRP",12.01,6,2.3);
   mat->SetUniqueID(13);
   mix = new TGeoMixture("H2O",2,   1.00000    );
      mix->DefineElement(0,1.01,1,0.1120977);
      mix->DefineElement(1,16,8,0.8879023);
   mix->SetUniqueID(14);
   mix = new TGeoMixture("Ethane_Gas",2,  0.135600E-02);
      mix->DefineElement(0,12.01,6,0.7985373);
      mix->DefineElement(1,1.01,1,0.2014628);
   mix->SetUniqueID(15);
   mat = new TGeoMaterial("RHONEYCM",26.98,13,0.601);
   mat->SetUniqueID(16);
//---> create mediums
TGeoMedium *med0 = new TGeoMedium("Vacuum",0,0,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med1 = new TGeoMedium("Be",1,1,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med2 = new TGeoMedium("Al",2,2,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med3 = new TGeoMedium("Fe",3,3,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med4 = new TGeoMedium("Cu",4,4,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med5 = new TGeoMedium("C",5,5,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med6 = new TGeoMedium("Pb",6,6,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med7 = new TGeoMedium("Si",7,7,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med8 = new TGeoMedium("scint",8,8,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med9 = new TGeoMedium("Li",9,9,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med10 = new TGeoMedium("N",10,10,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med11 = new TGeoMedium("Ne",11,11,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med12 = new TGeoMedium("Ts",12,12,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med13 = new TGeoMedium("CFRP",13,13,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med14 = new TGeoMedium("H2O",14,14,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med15 = new TGeoMedium("Ethane gas",15,15,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);
TGeoMedium *med16 = new TGeoMedium("RHONEYCM",16,16,0,0,0,20,0.1000000E+11,0.212,0.1000000E-02,1.150551);

//---> create volumes
   TGeoVolume *TOP = geom->MakeBox("TOP", med0, 200, 200, 200);
   TGeoVolume *BOX = geom->MakeBox("BOX", med1, 10,20,25);
   TGeoVolume *TRD1 = geom->MakeTrd1("TRD1", med2, 20,10,20,15);
   TGeoVolume *TRD2 = geom->MakeTrd2("TRD2", med3, 20,5,10,25,25);
   TGeoVolume *PARA = geom->MakePara("PARA", med4, 10, 20, 30, 15, 15, 120);
   Double_t v[16];
   v[0]=-22; v[1]=-18; v[2]=-18; v[3]=22; v[4]=22; v[5]=18; v[6]=18; v[7]=-22;
   v[8]=-12; v[9]=-8; v[10]=-8; v[11]=12; v[12]=12; v[13]=8; v[14]=8; v[15]=-12;
   TGeoVolume *ARB8 = geom->MakeArb8("ARB8", med5, 15, v);
   TGeoVolume *SPHE = geom->MakeSphere("SPHE", med6, 5, 15, 45, 180, 0, 270);
   TGeoVolume *TUBE = geom->MakeTube("TUBE", med7, 20, 25, 30);
   TGeoVolume *TUBS = geom->MakeTubs("TUBS", med8, 10,15,20, 45, 270);
   TGeoVolume *ELTU = geom->MakeEltu("ELTU", med9, 10,20,10);
   TGeoVolume *CTUB = geom->MakeCtub("CTUB", med10, 25, 30, 50, 0,270, 0, 0.5, -0.5*TMath::Sqrt(3.),0.5,0,0.5*TMath::Sqrt(3.));
   TGeoVolume *CONE = geom->MakeCone("CONE", med11, 30, 25, 30, 10, 15);
   TGeoVolume *CONS = geom->MakeCons("CONS", med12, 30, 25, 30, 10, 15, -45, 180);
   TGeoVolume *PCON = geom->MakePcon("PCON", med13, 0, 360, 3);
   TGeoPcon *pcon = (TGeoPcon*)(PCON->GetShape());
   pcon->DefineSection(0,-25, 10, 15);
   pcon->DefineSection(1,0, 10, 15);
   pcon->DefineSection(2,25, 25, 30);
   TGeoVolume *PGON = geom->MakePgon("PGON", med14, 0, 270, 4, 3);
   pcon = (TGeoPcon*)(PGON->GetShape());
   pcon->DefineSection(0,-25, 10, 15);
   pcon->DefineSection(1,0, 10, 15);
   pcon->DefineSection(2,25, 15, 20);
   TGeoVolume *TRAP = geom->MakeTrap("TRAP", med15, 25, 15, 30, 20, 15,10,15, 20,15,10,15);
   TGeoVolume *GTRA = geom->MakeGtra("GTRA", med16, 25, 15, 30, 30, 20, 15,10,15, 20,15,10,15);
   //---> create nodes
   geom->SetTopVolume(TOP);
   TOP->AddNode(BOX, 1);
   TOP->AddNode(TRD1, 2,  new TGeoTranslation(100, 0, 0));
   TOP->AddNode(TRD2, 3,  new TGeoTranslation(-100, 0, 0));
   TOP->AddNode(PARA, 4,  new TGeoTranslation(0, 100, 0));
   TOP->AddNode(ARB8, 5,  new TGeoTranslation(0, -100, 0));
   TOP->AddNode(SPHE, 6,  new TGeoTranslation(0, 0, 100));
   TOP->AddNode(TUBE, 7,  new TGeoTranslation(0, 0, -100));
   TOP->AddNode(TUBS, 8,  new TGeoTranslation(100, 0, 100));
   TOP->AddNode(ELTU, 9,  new TGeoTranslation(100, 0, -100));
   TOP->AddNode(CTUB, 10, new TGeoTranslation(-100, 0, 100));
   TOP->AddNode(CONE, 11, new TGeoTranslation(-100, 0, -100));
   TOP->AddNode(CONS, 12, new TGeoTranslation(0, 100, 100));
   TOP->AddNode(PCON, 13, new TGeoTranslation(0, -100, 100));
   TOP->AddNode(PGON, 14, new TGeoTranslation(100, 100, 100));
   TOP->AddNode(TRAP, 15, new TGeoTranslation(-100, -100, -100));
   TOP->AddNode(GTRA, 16, new TGeoTranslation(100, 100, 0));
   //---> close geometry
   geom->CloseGeometry("d");
   geom->DefaultColors();

   TIter next(gGeoManager->GetListOfVolumes());
   TGeoVolume *vol = (TGeoVolume*)next();
   Int_t ivol=1;
   printf("=== testing shapes ...\n");
   while ((vol=(TGeoVolume*)next())) {
      sample_volume(ivol);
      ivol++;
   }
   printf("=== testing global tracking ...\n");
   length();

   // print ROOTMARKs
   printf("\n");
   gBenchmark->Stop("stressShapes");
   gBenchmark->Print("stressShapes");
   Float_t ct = gBenchmark->GetCpuTime("stressShapes");
   Float_t cp_brun   = 5.55;
   Float_t rootmarks = 860*cp_brun/ct;
   printf("*******************************************************************\n");
   printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d CP = %7.2fs\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime(),ct);
   printf("*******************************************************************\n");
}

