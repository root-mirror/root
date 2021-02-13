/// \file
/// \ingroup tutorial_geom
/// Drawing the Cheongwadae building which is the Presidential Residence of the Republic of Korea, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x cheongwadae.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_cheongwadae.png width=800px
/// \macro_code
///
/// \author Hee Jun Shin (s-heejun@hanmail.net), Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void cheongwadae()
{
   TGeoManager *geom = new TGeoManager("geom","My first 3D geometry");

   //material
   TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
   TGeoMaterial *Fe = new TGeoMaterial("Fe",55.845,26,7.87);

   //creat media
   TGeoMedium *Air = new TGeoMedium("Vacuum",0,vacuum);
   TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);

   //creat volume
   TGeoVolume *top = geom->MakeBox("top",Air,300,300,300);
   geom->SetTopVolume(top);
   geom->SetTopVisible(0);
   // If you want to see the boundary, please input the number, 1 instead of 0.
   // Like this, geom->SetTopVisible(1);

char nBlocks[100];
int N = 0;
int f=0;
int di[2]; di[0] = 0; di[1] = 30;
TGeoVolume *mBlock;

   for(int k=0;k<7;k++){
      for(int i=0;i<20;i++){
         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.6,1.8,63);
         mBlock->SetLineColor(20);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-10.6-(2.6*i),-17.8+(6*k),0));

         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.7,1.8,58);
         mBlock->SetLineColor(12);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-11.9-(2.6*i),-17.8+(6*k),0));
      }
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 26,1.2,63);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-36,-14.8+(6*k),0));
   }
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 26,2,63);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-36,-21.6,0));

   for(int k=0;k<7;k++){
      for(int i=0;i<20;i++){
         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.6,1.8,63);
         mBlock->SetLineColor(20);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-10.6-(2.6*i),-17.8+(6*k),0));
         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.7,1.8,58);
         mBlock->SetLineColor(12);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-11.9-(2.6*i),-17.8+(6*k),0));

      }
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 26,1.2,63);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-36,-14.8+(6*k),0));
   }

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,22,58);
   mBlock->SetLineColor(2);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,0));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 3.5,8,0.1);
   mBlock->SetLineColor(13);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(4,-14,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 3.5,8,0.1);
   mBlock->SetLineColor(13);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-4,-14,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,20,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,17,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,14,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,11,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,8,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,5,60));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.2,0.1);
   mBlock->SetLineColor(1);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,2,60));

   for(int k=0;k<7;k++){
      for(int i=0;i<20;i++){
         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.6,1.8,63);
         mBlock->SetLineColor(20);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(10.6+(2.6*i),-17.8+(6*k),0));
         sprintf(nBlocks,"f%d_bg%d",f,N++);
         mBlock = geom->MakeBox(nBlocks, Iron, 0.7,1.8,58);
         mBlock->SetLineColor(12);
         top->AddNodeOverlap(mBlock,1,new TGeoTranslation(11.9+(2.6*i),-17.8+(6*k),0));

      }
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 26,1.2,63);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(36,-14.8+(6*k),0));
   }
   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 26,2,63);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(36,-21.6,0));


   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 82,2,82);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,24,0));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 85,0.5,85);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,26,0));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 88,2,88);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-24,0));


   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 30, 0, 180, 0, 180);
   mBlock->SetLineColor(32);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,24,0));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 0.1,30,0.1);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,40,0));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 0,30,4,360,360);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,27,0, new TGeoRotation("r1",0,90,0)));

   for(int i=0;i<8;i++){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2,22,2);
      mBlock->SetLineColor(18);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-70+(20*i),0,80));
   }

   for(int i=0;i<8;i++){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2,22,2);
      mBlock->SetLineColor(18);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-70+(20*i),0,-80));
   }

   for(int i=0;i<7;i++){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2,22,2);
      mBlock->SetLineColor(18);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-70,0,-80+(23*i)));
   }

   for(int i=0;i<7;i++){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2,22,2);
      mBlock->SetLineColor(18);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(70,0,-80+(23*i)));
   }

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 100,0.5,160);
   mBlock->SetLineColor(41);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-26,40));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 10,0.01,160);
   mBlock->SetLineColor(19);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-25,40));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-22,170));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-25,170));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-22,150));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-25,150));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-22,130));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-25,130));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-22,110));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(15,-25,110));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-22,170));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-25,170));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-22,150));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-25,150));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-22,130));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-25,130));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-22,110));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeSphere(nBlocks, Iron, 0, 5, 0, 180, 0, 180);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-25,110));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 0.1,10,0.1);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(20,-15,110));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 5,3,0.1);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(25,-8,110));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 0.1,10,0.1);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-20,-15,110));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 5,3,0.1);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-15,-8,110));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 7,1.5,5);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-24,88));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 7,1,5);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-24,92));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 7,0.5,5);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-24,96));

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 7,0.1,5);
   mBlock->SetLineColor(18);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,-24,100));

   geom->CloseGeometry();
   top->SetVisibility(0);

   top->Draw("ogl");
}
