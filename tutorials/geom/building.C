/// \file
/// \ingroup tutorial_geom
/// Drawing a building where Dept. of Physics is, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x building.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_building.png width=800px
/// \macro_code
///
/// \author Hyung Ju Lee (laccalus@nate.com), Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void building()
{
   TGeoManager *geom = new TGeoManager("geom","My First 3D Geometry");

// Materials
   TGeoMaterial *Vacuum = new TGeoMaterial("vacuum",0,0,0);
   TGeoMaterial *Fe = new TGeoMaterial("Fe",55.845,26,7.87);

// Media
   TGeoMedium *Air = new TGeoMedium("Air",0,Vacuum);
   TGeoMedium *Iron = new TGeoMedium("Iron",0,Fe);

// Volume
   TGeoVolume *Phy_Building = geom->MakeBox("top",Air,150,150,150);
   geom->SetTopVolume(Phy_Building);
   geom->SetTopVisible(0);
   // If you want to see the boundary, please input the number, 1 instead of 0.
   // Like this, geom->SetTopVisible(1);


   TGeoVolume *mBlocks;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////         Front-Building        ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   int i = 0;
   int F = 0;           // Floor
   int N = 0;           // Block_no
   int nW = 8;          // Number of windows
   int nF = 3;          // Number of Floor
   char no_Block[100];  // Name of Block
   double  sP = 0;      // Starting Phi of Tubs
   double  hP = 21;     // Height of Tubs from Ground

   while (F<nF){
      N = 0; i = 0; sP = 0;

/////////////////////////// Front of Building
   while (i<nW){
      i++;
      sprintf(no_Block, "B1_F%d_%d", F, ++N);                  // Windows (6.25)
      mBlocks = geom->MakeTubs(no_Block,Iron,21,29,1.8,sP,sP+6.25);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

      if (i < nW) {
         sprintf(no_Block, "B1_F%d_%d", F, ++N);               // Walls (8)
         mBlocks = geom->MakeTubs(no_Block,Iron,21,30,1.8,sP,sP+2.5);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

         sprintf(no_Block, "B1_F%d_%d", F, ++N);
         mBlocks = geom->MakeTubs(no_Block,Iron,21,31,1.8,sP,sP+1);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

         sprintf(no_Block, "B1_F%d_%d", F, ++N);
         mBlocks = geom->MakeTubs(no_Block,Iron,21,30,1.8,sP,sP+1);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

         sprintf(no_Block, "B1_F%d_%d", F, ++N);
         mBlocks = geom->MakeTubs(no_Block,Iron,21,31,1.8,sP,sP+1);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

         sprintf(no_Block, "B1_F%d_%d", F, ++N);
         mBlocks = geom->MakeTubs(no_Block,Iron,21,30,1.8,sP,sP+2.5);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));
      }

      if (i>=nW) {
         sprintf(no_Block, "B1_F%d_%d", F, ++N);               // Walls
         mBlocks = geom->MakeTubs(no_Block,Iron,21,30,1.8,sP,103);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP+(8*F)));

      }
   }

   sprintf(no_Block, "B1_F%d", ++F);                  // No Windows Floor
   mBlocks = geom->MakeTubs(no_Block,Iron,21,30,2.2,0,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-12+(8*F)));

///////////////////////////////////////// Back of Building
   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.5,21,0.8,92,101);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-9.4+(8*F)));

   if(F<nF){
      sprintf(no_Block, "B1_F%d_%d", F, ++N);
      mBlocks = geom->MakeTubs(no_Block,Iron,18.5,21,3.2,92,102);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-5.4+(8*F)));

   }
   }

   sprintf(no_Block, "B1_F%d_%d", F, ++N);               // Walls
   mBlocks = geom->MakeTubs(no_Block,Iron,18.5,21,2,92,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-4));
   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.5,21,3.2,92,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-5.4+(8*F)));

   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,30,2,0,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-4.2+(8*F)));

   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,21,2,0,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-4.2+(8*F)));

   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,18.5,14,92,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,29));

//////////////////////// Front of Building
   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,29,2,0,97);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,13));

   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,32,2,37,97);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,9));

   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,29,1.95,0,37);
   mBlocks->SetLineColor(30);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,9.05));
   sprintf(no_Block, "B1_F%d_%d", F, ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,29,0.05,0,37);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,7.05));


//////////////////////// Rooftop
   sprintf(no_Block, "B1_RT%d", N = 0);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,29.5,0.2,0,102);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-2+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.5,21,0.2,0,101);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-2+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,30,0.7,102.9,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21.1,29.9,0.7,102,102.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21.1,21.5,0.5,98,102.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,21,21.1,0.7,98,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.6,21,0.7,101.9,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.6,21,0.7,101,101.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.8+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,29.5,29.9,0.5,0,102);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,29.9,30,0.5,0,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.1,18.5,0.5,-1,101.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,18.1,0.5,-0.5,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));

   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.1,18.4,0.5,101.9,102.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,18.1,0.5,102,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18.4,18.5,0.5,102,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));
   sprintf(no_Block, "B1_RT%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,18.5,0.5,102.9,103);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,hP-1.7+(8*F)));


/////////////////////////////// White Wall
   sprintf(no_Block, "B1_WW%d", N = 0);
   mBlocks = geom->MakeTubs(no_Block,Iron,20.8,31,19.5,sP,sP+1);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,26));

   sprintf(no_Block, "B1_WW%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,26.8,31,5,sP,sP+1);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,2));

   sprintf(no_Block, "B1_WW%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,23,24.3,5,sP,sP+1);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,2));

   sprintf(no_Block, "B1_WW%d", ++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,20.8,21.3,5,sP,sP+1);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,2));



////////////////////////// Zero Floor1
   sprintf(no_Block, "B1_ZF%d",N=0);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,21,9,0,92);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,6));

   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,21,7.5,0,92);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,31.5));

   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,18,21,4.5,0,92);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,19.5));

   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,18,0.2,0,101);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,18.6));
   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,18,1.7,0,100);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,16.7));
   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,18,1.2,101,101.9);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,19.6));
   sprintf(no_Block, "B1_ZF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,18,1.2,101.9,102);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,19.6));


////////////////////////// Zero Floor2
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6.5,7,2.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,10.75,13));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6.5,7,3);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,10.75,7.5));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,0.05,10);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,17.95,7));
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6.9,0.20,10);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,17.70,7));
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.1,0.20,10);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-13.9,17.70,7));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.05,7,3.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-13.95,10.5,13.5));
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.20,6.9,3.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-13.70,10.5,13.5));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.25,7,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-13.75,10.5,1));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,0.05,10);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,3.55,7));
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6.9,0.20,10);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7,3.8,7));
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.1,0.20,10);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-13.9,3.8,7));


////////////////////////// Zero Floor2
   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5,5,1);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-5,23,-2));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5,0.25,1.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-5,28.25,-1.5));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.25,5.5,1.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-10.25,23,-1.5));

   sprintf(no_Block, "B1_ZF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,3.5,5);
   mBlocks->SetLineColor(20);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-12.5,0,-4));


////////////////////////// Ground
   sprintf(no_Block, "B1_GRD%d",N=0);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,29,1,0,36.75);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-2));

   sprintf(no_Block, "B1_GRD%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,30.4,0.4,36.75,77.25);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-2.7));

   sprintf(no_Block, "B1_GRD%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,29.7,0.3,36.75,77.25);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-2));

   sprintf(no_Block, "B1_GRD%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,29,0.3,36.75,77.25);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-1.3));

   sprintf(no_Block, "B1_GRD%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,0,29,1,77.25,97);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-2));


///////////////////////////// Pillars & fences
   sprintf(no_Block, "B1_PF%d", N = 0);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(29,4.2,6, new TGeoRotation("r1",8.25,0,0)));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(24.2,16.5,6, new TGeoRotation("r1",34.25,0,0)));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(14.5,25.4,6, new TGeoRotation("r1",60.25,0,0)));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(1.9,29.2,6, new TGeoRotation("r1",86.25,0,0)));

   sprintf(no_Block, "B1_PF%d",++N);
   mBlocks = geom->MakeTubs(no_Block,Iron,29,30,2,0,36.75);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,0,-1));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3,2,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-0.75,29.3,-1));

   sprintf(no_Block, "B1_PF%d", ++N);      //장애인용
   mBlocks = geom->MakeBox(no_Block,Iron,0.25,4.3,1.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(6.5,30.6,-1.5));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.25,4.3,0.4);
   mBlocks->SetLineColor(10);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(1.125,30.6,-2.7));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,0.25,0.75);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(1.125,34.9,-2.25));

   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeTrd1(no_Block,Iron,1.5,0,0.25,5.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(1.125,34.9,-1.5, new TGeoRotation("r1",90,-90,90)));









///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Second Part of Front-Building ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   F=0;
   while (F<nF){
      N = 0; i = 0; nW = 7;

   while (i<nW){
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);      // Wall
      mBlocks = geom->MakeBox(no_Block,Iron,3.8,0.35,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (21.65-6*i)*0.13,-21.2 + (21.65-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,4.8,0.3,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (21.0-6*i)*0.13,-21.2 + (21-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3.8,0.3,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (20.4-6*i)*0.13,-21.2 + (20.4-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,4.8,0.3,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (19.7-6*i)*0.13,-21.2 + (19.7-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3.8,0.35,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (19.05-6*i)*0.13,-21.2 + (19.05-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));


      sprintf(no_Block, "B12_F%d_B%d",F, ++N);      // Windows
      mBlocks = geom->MakeBox(no_Block,Iron,3,1.4,1.8);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(23.38 + (17.4-6*i)*0.13,-21.2 + (17.4-6*i)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      i++;
      if( i >= nW){
      sprintf(no_Block, "B12_F%d_B%d",F, ++N);      // Wall.
      mBlocks = geom->MakeBox(no_Block,Iron,5.8,1,1.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,
         new TGeoCombiTrans(21.4 + (-21)*0.13,-21 + (-21)*0.99,hP+(8*F),
         new TGeoRotation("r1",-7.5,0,0)));
      }
   }

      sprintf(no_Block, "B12_F%d_B%d",++F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,5.8,22,2.2);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.4,-21,hP-12+(8*F), new TGeoRotation("r1",-7.5,0,0)));

   }
   sprintf(no_Block, "B12_F%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.8,22,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.4,-21,hP-4.2+(8*F), new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_F%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,22,14);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(18.43,-20.61,29, new TGeoRotation("r1",-7.5,0,0)));


////////////////////// RoofTop
   sprintf(no_Block, "B12_RT%d_%d", F, N=0);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,21.75,0.2);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.43,-20.75,hP-2+(8*F), new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.23,21.95,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(26.9,-21.72,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.1,22,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.1,-21.75,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));


   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.23,3.6,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(13.65,-38.03,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.02,3.8,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(13.3,-38.39,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));



   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.7,0.23,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(18.57,-42.48,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.8,0.1,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(18.54,-42.71,hP-1.7+(8*F), new TGeoRotation("r1",-7.5,0,0)));


//////////////////////// Pillars & fences
   sprintf(no_Block, "B12_PF%d", N = 0);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(28.32,-7.44,6, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(26.75,-19.33,6, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,9);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(25.19,-31.23,6, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,11);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(23.75,-42.14,4, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,11);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(13.84,-40.83,4, new TGeoRotation("r1",-7.5,0,0)));



   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,15.75,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.42,-15.48,-1, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,2,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(24.28,-39.27,-3, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.5,15.75,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(28.91,-15.68,-4, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_RT%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.8,0.5,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(18.8,-40.73,-3, new TGeoRotation("r1",-7.5,0,0)));


/////////////////////// Stair
   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3,0.5,3.25);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(28.33,-31.49,-2.75, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,6.25,1.625);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(30.56,-37.58,-4.375, new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B1_PF%d", ++N);
   mBlocks = geom->MakeTrd1(no_Block,Iron,3.25,0,0.5,6.25);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(30.56,-37.58,-2.75, new TGeoRotation("r1",-7.5,90,90)));


   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3,3,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.37,-34.89,-2.5, new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.5,3,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.74,-35.95,-3.5, new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.5,3,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.61,-36.94,-4.5, new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.5,3,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(27.48,-37.93,-5.5, new TGeoRotation("r1",-7.5,0,0)));



//////////////////////// Ground
   sprintf(no_Block, "B12_GR%d", N=0);
   mBlocks = geom->MakeBox(no_Block,Iron,4.8,21,1);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.53,-20.1,-2, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_GR%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.8,18,9);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(12.86,-16.62,6, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_GR%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.8,22,2);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.4,-21,13, new TGeoRotation("r1",-7.5,0,0)));

   sprintf(no_Block, "B12_GR%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.8,22,1.95);
   mBlocks->SetLineColor(30);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.4,-21,9.05, new TGeoRotation("r1",-7.5,0,0)));
   sprintf(no_Block, "B12_GR%d_%d", F, ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.8,22,0.05);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoCombiTrans(21.4,-21,7.05, new TGeoRotation("r1",-7.5,0,0)));







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////         Bridge-Building       ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   F=1; N = 0; nF = 4;


   sprintf(no_Block, "B2_F%d", 6);
   mBlocks = geom->MakeBox(no_Block,Iron,7,17.5,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(12,-17.5,41));

   while(F++ <= nF){
//////////////////////// Front
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,0.8,4,4);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(10,-4,-5 +(F*8)));

      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,3.5,1);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.9,-4,-2 +(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,4.5,0.2);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.9,-4,-3.2+(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,4,2.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.9,-4,-6.2+(F*8)));

      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,0.7,4,4);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(13.6,-4,-5 +(F*8)));

      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,3.5,1);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(15.4,-4,-2 +(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,4.5,0.2);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(15.4,-4,-3.2+(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.1,4,2.8);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(15.4,-4,-6.2+(F*8)));

      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,0.7,4,4);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(17.1,-4,-5 +(F*8)));


//////////////////////////// Back
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.3,13.5,1.5);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(6.8,-21.5,-2.5 +(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.3,14,0.2);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(6.8,-21.5,-4.2+(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.3,13.5,2.3);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(6.8,-21.5,-6.8+(F*8)));



      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3.7,13,1.5);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.8,-21.5,-2.5 +(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3.7,14,0.2);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.8,-21.5,-4.2+(F*8)));
      sprintf(no_Block, "B2_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3.7,13.5,2.3);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(11.8,-21.5,-6.8+(F*8)));


   }


   sprintf(no_Block, "B2_F%d_%d", 0,1);
   mBlocks = geom->MakeBox(no_Block,Iron,5,13.5,6);
   mBlocks->SetLineColor(30);
//   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(10,-21.5,-15));
   sprintf(no_Block, "B2_F%d_%d", 0,2);
   mBlocks = geom->MakeBox(no_Block,Iron,5,13.5,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(10,-21.5,-5));
   sprintf(no_Block, "B2_F%d_%d", 0,3);
   mBlocks = geom->MakeBox(no_Block,Iron,5,13.5,4);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(10,-21.5,3));



/////////////////////////// RoofTop
   sprintf(no_Block, "B2_F%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,17.4,0.1);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(12,-17.5,43.1));

   sprintf(no_Block, "B2_F%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,0.2,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(9.5,-34.7,43.5));
   sprintf(no_Block, "B2_F%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,0.05,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(9.5,-34.95,43.5));

   sprintf(no_Block, "B2_F%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.75,0.2,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(13.75,-0.3,43.5));
   sprintf(no_Block, "B2_F%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.55,0.05,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(13.55,-0.05,43.5));









///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////         Building 3            ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   F=0; N = 0; nF = 4; nW = 6;


   sprintf(no_Block, "B3_F0%d", 7);
   mBlocks = geom->MakeBox(no_Block,Iron,3,36,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,49));

   while (F++ < nF){
      i=0; N=0;

      sprintf(no_Block, "B3_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,4,36,0.2);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,14.7 +(F*8)));

      while (i++ <nW){
         sprintf(no_Block, "B3_F%d_%d",F, ++N);
         mBlocks = geom->MakeBox(no_Block,Iron,2.5,5,1.8);
         mBlocks->SetLineColor(12);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-6 -(i*12),12.8 +(F*8)));

         sprintf(no_Block, "B3_F%d_%d",F, ++N);
         mBlocks = geom->MakeBox(no_Block,Iron,2.8,1,1.8);
         mBlocks->SetLineColor(18);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-12 -(i*12),12.8 +(F*8)));

      }

      sprintf(no_Block, "B3_F%d_%d",F, ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3,36,2);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,9.2 +(F*8)));

   }

   sprintf(no_Block, "B3_F0%d", 1);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,36,2);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,13));
   sprintf(no_Block, "B3_F0%d", 2);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,36,2);
   mBlocks->SetLineColor(30);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,9));

   sprintf(no_Block, "B3_F0%d", 3);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,36,4);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,3));

   sprintf(no_Block, "B3_F0%d", 4);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,36,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,-5));

/*   sprintf(no_Block, "B3_F0%d", 5);
   mBlocks = geom->MakeBox(no_Block,Iron,2.8,36,6);
   mBlocks->SetLineColor(30);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-48,-15));
*/
   sprintf(no_Block, "B3_F0%d", 61);
   mBlocks = geom->MakeBox(no_Block,Iron,3,8,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-88,49));

   sprintf(no_Block, "B3_F0%d", 62);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,8,24);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.9,-88,23));
   sprintf(no_Block, "B3_F0%d", 63);
   mBlocks = geom->MakeBox(no_Block,Iron,2,7,24);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-88,23));
   sprintf(no_Block, "B3_F0%d", 64);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,8,24);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-0.1,-88,23));

   sprintf(no_Block, "B3_F0%d", 65);
   mBlocks = geom->MakeBox(no_Block,Iron,3,8,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-88,-5));

/////////////////////////////// Left-Side
   nF = 6;nW = 6;

   sprintf(no_Block, "B3_F2%d",7);
   mBlocks = geom->MakeBox(no_Block,Iron,7,40.5,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-43.5,49));

   for (F=0 ; F<nF ; F++){ N=0;
      for (i = 0 ; i<nW ; i++){
         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,6,2.35,2);
         mBlocks->SetLineColor(12);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-14.35-(12*i),5 + (8*F)));
         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,6.5,0.3,2);
         mBlocks->SetLineColor(18);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-17-(12*i),5 + (8*F)));
         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,6,2.35,2);
         mBlocks->SetLineColor(12);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-19.65-(12*i),5 + (8*F)));

         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,7,1,2);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-23-(12*i),5 + (8*F)));
      }

      sprintf(no_Block, "B3_F2%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,6.8,36,0.3);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-48,3.3 + (8*F)));

      sprintf(no_Block, "B3_F2%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,7,36,2);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-48,1 + (8*F)));

      for(int i=0;i<4;i++){
         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,7,0.5,1.4);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-3.5,5.6 + (8*F)));

         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,6,0.7,1.4);
         mBlocks->SetLineColor(12);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-4.7,5.6 + (8*F)));

         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,7,1.6,1.4);
         mBlocks->SetLineColor(2);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-7,5.6 + (8*F)));

         sprintf(no_Block, "B3_F2%d_%d",F,++N);
         mBlocks = geom->MakeBox(no_Block,Iron,6,0.7,1.4);
         mBlocks->SetLineColor(12);
         Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-9.3,5.6 + (8*F)));
      }

      sprintf(no_Block, "B3_F2%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,7,3.5,2.6);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-6.5,1.6 + (8*F)));
   }

   sprintf(no_Block, "B3_F2%d",71);
   mBlocks = geom->MakeBox(no_Block,Iron,7,40.5,4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-43.5,-5));

   sprintf(no_Block, "B3_F2%d",72);
   mBlocks = geom->MakeBox(no_Block,Iron,7,2,30);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-86,21));

   sprintf(no_Block, "B3_F2%d",73);
   mBlocks = geom->MakeBox(no_Block,Iron,7,1,30);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.4,-11,21));



//////////////////////////////////// Rooftop
   sprintf(no_Block, "B3_RT%d",N = 0);
   mBlocks = geom->MakeBox(no_Block,Iron,7,42.25,0.1);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.15,-45.5,51.1));
   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.75,41.75,0.1);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-54,51.1));

   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.24,41.99,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5.15,-53.99,51.5));
   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.01,42,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5.4,-54,51.5));

   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.24,3.99,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-0.35,-92,51.5));
   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.01,4,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-0.6,-92,51.5));

   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,2.99,0.24,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-95.79,51.5));
   sprintf(no_Block, "B3_RT%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(2.4,-96.04,51.5));

   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.24,42.49,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-14.14,-45.5,51.5));
   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.01,42.5,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-14.39,-45.5,51.5));


/////////////////////// Stair
   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6.99,0.24,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.15,-3.25,51.5));
   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.15,-3,51.5));

   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,0.25,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.15,-87.74,51.5));
   sprintf(no_Block, "B3_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,7,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-7.15,-87.99,51.5));



/////////////////////////////// Pillars
   N=0;
   for (i=0 ; i<6; i++) {
      sprintf(no_Block, "B3_PF%d", ++N);
      mBlocks = geom->MakeBox(no_Block,Iron,1.2,1.5,12);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.6,-12-(12*i),3));
   }
   sprintf(no_Block, "B3_PF%d", ++N);
   mBlocks = geom->MakeBox(no_Block,Iron,1.5,40,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(7,-56,-5));


////////////////////////////// Stair
   sprintf(no_Block, "B3_ST%d",N=0);
   mBlocks = geom->MakeBox(no_Block,Iron,0.5,7,5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-6.5,-88,-2));

   for(int i=0;i<5;i++){
      sprintf(no_Block, "B3_ST%d",++N);
      mBlocks = geom->MakeBox(no_Block,Iron,3,5,0.5);
      mBlocks->SetLineColor(18);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3,-86-(0.7*i),-2-(1*i)));
   }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////           Mid-Building         ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////// Left-Side

   for(int F=0;F<5;F++){ N=0;
      sprintf(no_Block, "B4_LF%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,5.5,12.5,2.3);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-7.5,9.6+(8*F)));

      sprintf(no_Block, "B4_LF%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,5.5,2,1.7);
      mBlocks->SetLineColor(2);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,3,13.6+(8*F)));

      sprintf(no_Block, "B4_LF%d_%d",F,++N);
      mBlocks = geom->MakeBox(no_Block,Iron,5,10.5,1.7);
      mBlocks->SetLineColor(12);
      Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-9.5,13.6+(8*F)));
   }

   sprintf(no_Block, "B4_LF%d_%d",9,N=0);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,12.5,6);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-7.5,53));

   sprintf(no_Block, "B4_LF%d_%d",9,++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,2,4.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,3,3));
   sprintf(no_Block, "B4_LF%d_%d",9,++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5,10.5,4.5);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-9.5,3));

   sprintf(no_Block, "B4_LF%d_%d",9,++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,12.5,5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-7.5,-3));




///////////////////////////////////// Right-Side
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.25,11,24);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.25,-9,19));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.25,4,32);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(8.75,2,27));


   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,2,1.8);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,0,44.8));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,3.5,5);
   mBlocks->SetLineColor(20);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-12.5,0,-4));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,6,2,0.3);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-4,46.3));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4,2,1.5);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-4,44.5));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,7,1.8);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-13,44.8));

   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,11,1.8);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-9,48.4));

   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.25,1.5,2);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-0,52.2));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4,2,2);
   mBlocks->SetLineColor(12);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-4,52.2));
   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,7,2);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-13,52.2));


   sprintf(no_Block, "B4_RS%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,4.5,11,2.4);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-9,56.6));

///////////////////////////////// RoofTop
   sprintf(no_Block, "B4_RT%d",N=0);
   mBlocks = geom->MakeBox(no_Block,Iron,4.25,10.9,0.2);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(4.5,-9,59));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.25,12.4,0.2);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,-7.5,59));

   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.24,12.4,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-8.79,-7.5,59.3));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.01,12.4,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-9.06,-7.5,59.3));

   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.24,13,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(8.75,-7,59.3));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,0.01,13,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(9,-7,59.3));

   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,8.75,0.24,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,-19.75,59.3));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,8.75,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(0,-20.01,59.3));

   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.25,0.24,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.5,4.55,59.3));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,5.5,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(-3.75,5.1,59.3));

   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3.5,0.24,0.5);
   mBlocks->SetLineColor(18);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5,1.55,59.3));
   sprintf(no_Block, "B4_RT%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,3.5,0.01,0.5);
   mBlocks->SetLineColor(2);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5,2.1,59.3));

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////             GROUND             ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   sprintf(no_Block, "GRD%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,40,90,2);
   mBlocks->SetLineColor(30);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5,-20,-9));

   sprintf(no_Block, "GRD%d",++N);
   mBlocks = geom->MakeBox(no_Block,Iron,30,30,2);
   mBlocks->SetLineColor(41);
   Phy_Building->AddNodeOverlap(mBlocks,1,new TGeoTranslation(5,30,-5));
   geom->CloseGeometry();




////////////////////////// Draw
   Phy_Building->SetVisibility(0);
   Phy_Building->Draw("ogl");


}
