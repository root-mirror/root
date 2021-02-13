/// \file
/// \ingroup tutorial_geom
/// Drawing a famous Korean gate, the South gate, called Namdeamoon in Korean, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x south_gate.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_south_gate.png width=800px
/// \macro_code
///
/// \author Lan Hee Yang(yangd5d5@hotmail.com), Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void south_gate()
{


  TGeoManager *geom=new TGeoManager("geom","My first 3D geometry");


  TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);//a,z,rho
  TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);

  //Create media

  TGeoMedium *Air  = new TGeoMedium("Vacuum",0,vacuum);
  TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);

  //Create volume

  TGeoVolume *top = geom->MakeBox("top",Air,1000,1000,1000);
  geom->SetTopVolume(top);
  geom->SetTopVisible(0);
  // If you want to see the boundary, please input the number, 1 instead of 0.
  // Like this, geom->SetTopVisible(1);


//base

char nBlocks[100];
int i=1;
int N = 0;
int f=0;
int di[2]; di[0] = 0; di[1] = 30;
TGeoVolume *mBlock;

while (f<11){
while (i<14){
   if (i==6 && f<8){
      i = i+3;
   }

   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 29,149,9);
   mBlock->SetLineColor(20);
   if (f<8){
   if (i<=5 && f<8){
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-120-((i-1)*60)-di[f%2],5,5+(20*f)));
   } else if (i>5 && f<8){
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(120+((i-9)*60)  +di[f%2],5,5+(20*f)));
   }
   } else {
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-420+(i*60)-di[f%2],5,5+(20*f)));
   }
   i++;
   if (i>=14 && f>=8 && f%2 == 1){
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 29,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-420+(i*60)-di[f%2],5,5+(20*f)));
   i++;
   }
   if (f%2 ==0){
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 14.5,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-405,5,5+(20*f)));
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 14.5,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(405,5,5+(20*f)));
   } else if (f<5){
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 14.5,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-105,5,5+(20*f)));
      sprintf(nBlocks,"f%d_bg%d",f,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 14.5,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(105,5,5+(20*f)));

   }
}
      sprintf(nBlocks,"f%d_bg%d",8,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 40,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-80,5,145));
      sprintf(nBlocks,"f%d_bg%d",8,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 40,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(80,5,145));

      sprintf(nBlocks,"f%d_bg%d",7,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 15,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-75,5,125));
      sprintf(nBlocks,"f%d_bg%d",7,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 15,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(75,5,125));

      sprintf(nBlocks,"f%d_bg%d",6,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 24,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-95,5,105));
      sprintf(nBlocks,"f%d_bg%d",6,N++);
      mBlock = geom->MakeBox(nBlocks, Iron, 24,149,9);
      mBlock->SetLineColor(20);
      top->AddNodeOverlap(mBlock,1,new TGeoTranslation(95,5,105));



i=1;f++;

}




//wall

f=0;
while (f<5){
i=0;
while (i<65){
   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 5.8,3,3.8);
   mBlock->SetLineColor(25);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-384+(i*12),137,218+(f*8)));
i++;

}
f++;
}



f=0;
while (f<5){
i=0;
while (i<65){
   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 5.8,3,3.8);
   mBlock->SetLineColor(25);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-384+(i*12),-137,218+(f*8)));
i++;

}
f++;
}





f=0;
while (f<7){
i=0;
while (i<22){
   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 3,5.8,3.8);
   mBlock->SetLineColor(25);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-384,-126+(i*12),218+(f*8)));
i++;

}
f++;
}



f=0;
while (f<7){
i=0;
while (i<22){
   sprintf(nBlocks,"f%d_bg%d",f,N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 3,5.8,3.8);
   mBlock->SetLineColor(25);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(384,-126+(i*12),218+(f*8)));
i++;

}
f++;
}


// arch


int k;
k=0; i=0;

while (i<5){
while(k<10){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 70,89,14, (i*36)+0.5, (i+1)*36-0.5);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,-130+(k*30),70, new TGeoRotation("r1",0,90,0)));
   k++;
}
   i++; k=0;
}

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 9,149,17);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(80,5,14));
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 9,149,18);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(80,5,51));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 9,149,17);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-80,5,14));
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 9,149,18);
   mBlock->SetLineColor(20);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(-80,5,51));





//wall's kiwa

k=0; i=0;

while (i<5){
while(k<52){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 1,3,7, 0, 180);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-382+(k*15),137,255, new TGeoRotation("r1",90,90,0)));
   k++;
}
   i++; k=0;
}





k=0; i=0;

while (i<5){
while(k<52){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 2.5,3,7, 0, 180);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-382+(k*15),-137,255, new TGeoRotation("r1",90,90,0)));
   k++;
}
   i++; k=0;
}



k=0; i=0;

while (i<5){
while(k<20){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 2.5,3,6, 0, 180);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-382,-123+(k*13),271, new TGeoRotation("r1",0,90,0)));
   k++;
}
   i++; k=0;
}




k=0; i=0;

while (i<5){
while(k<20){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 2.5,3,7, 0, 180);
   mBlock->SetLineColor(12);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(382,-123+(k*13),271, new TGeoRotation("r1",0,90,0)));
   k++;
}
   i++; k=0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// 1 floor


k=0; i=0;

while (i<5){
while(k<7){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 0,5,56, 0, 360);
   mBlock->SetLineColor(50);

   if (k<=2){

   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),80,260, new TGeoRotation("r1",0,0,0)));
   }else if (k>=4){
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),80,260, new TGeoRotation("r1",0,0,0)));
   }


k++;
}
   i++; k=0;
}



k=0; i=0;

while (i<5){
while(k<7){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeTubs(nBlocks,Iron, 0,5,56, 0, 360);
   mBlock->SetLineColor(50);

   if (k<=2){

   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),-80,260, new TGeoRotation("r1",0,0,0)));
   }else if (k>=4){
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),-80,260, new TGeoRotation("r1",0,0,0)));
   }


k++;
}
   i++; k=0;
}

// ||=====||======||=====||=====||=====||=====||


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 298,78,8);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,300));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 298,78,5);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,320));



//1
k=0; i=0;

while (i<5){
while(k<6){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron,18,10,8);
   mBlock->SetLineColor(8);
   {
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),70,300, new TGeoRotation("r1",0,0,0)));

   }
   k++;
}
   i++; k=0;
}





k=0; i=0;

while (i<5){
while(k<6){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron,18,10,8);
   mBlock->SetLineColor(8);
   {
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),-70,300, new TGeoRotation("r1",0,0,0)));

   }
   k++;
}
   i++; k=0;
}




   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,8);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-290,0,300, new TGeoRotation("r1",90,0,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,8);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(290,0,300, new TGeoRotation("r1",90,0,0)));






//2
k=0; i=0;

while (i<5){
while(k<6){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron,18,10,5);
   mBlock->SetLineColor(8);
   {
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),70,320, new TGeoRotation("r1",0,0,0)));

   }
   k++;
}
   i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,5);
      mBlock->SetLineColor(8);
      {
         top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),-70,320, new TGeoRotation("r1",0,0,0)));

      }
      k++;
   }
      i++; k=0;
}



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,5);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-290,0,320, new TGeoRotation("r1",90,0,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,5);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(290,0,320, new TGeoRotation("r1",90,0,0)));











//___||____||_____||____||____||____||____||


k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),78,345, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),-78,345, new TGeoRotation("r1",0,0,0)));
      k++;
   }
   i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,-78+(k*33),345, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,-78+(k*33),345, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}

//        ||//  ||//  ||//  ||//



k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),90,342, new TGeoRotation("r1",0,-45,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),-90,342, new TGeoRotation("r1",0,45,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-318,-78+(k*33),345, new TGeoRotation("r1",-90,45,0)));
      k++;
   }
   i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(318,-78+(k*33),345, new TGeoRotation("r1",90,45,0)));
      k++;
   }
   i++; k=0;
}


//   /// || / / / / / / / || / / / / / / / / || / / / / / / / / / / /



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 330,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,-107,362, new TGeoRotation("r1",0,-45,0)));




   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 330,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,107,362, new TGeoRotation("r1",0,45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 110,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(330,0,362, new TGeoRotation("r1",90,-45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 110,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-330,0,362, new TGeoRotation("r1",90,45,0)));




/////////////////////// add box




k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,2);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*100),-108,362, new TGeoRotation("r1",0,-45,0)));

      }
      k++;
   }
   i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,2);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*100),108,362, new TGeoRotation("r1",0,45,0)));

      }
      k++;
   }
   i++; k=0;
}


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,2);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(331,0,362, new TGeoRotation("r1",90,-45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,2);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-331,0,362, new TGeoRotation("r1",90,45,0)));





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 2nd floor


k=0; i=0;

while (i<5){
   while(k<7){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 0,5,30, 0, 360);
      mBlock->SetLineColor(50);

      if (k<=2){

      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),80,465, new TGeoRotation("r1",0,0,0)));
      }else if (k>=4){
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),80,465, new TGeoRotation("r1",0,0,0)));
      }


   k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<7){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 0,5,30, 0, 360);
      mBlock->SetLineColor(50);

      if (k<=2){

      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),-80,465, new TGeoRotation("r1",0,0,0)));
      }else if (k>=4){
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*100),-80,465, new TGeoRotation("r1",0,0,0)));
      }


   k++;
   }
   i++; k=0;
}




// ||=====||======||=====||=====||=====||=====||


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 302,80,8);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,480));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 302,80,5);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,500));


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 305,80,2.5);
   mBlock->SetLineColor(50);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,465));


///////////////////////add box






//1
k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,8);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),71,480, new TGeoRotation("r1",0,0,0)));

      }
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,8);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),-71,480, new TGeoRotation("r1",0,0,0)));

      }
      k++;
   }
      i++; k=0;
}




   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,8);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-293,0,480, new TGeoRotation("r1",90,0,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,8);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(293,0,480, new TGeoRotation("r1",90,0,0)));






//2
k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,5);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),71,500, new TGeoRotation("r1",0,0,0)));

      }
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<6){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron,18,10,5);
      mBlock->SetLineColor(8);
      {
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-250+(k*100),-71,500, new TGeoRotation("r1",0,0,0)));

      }
      k++;
   }
      i++; k=0;
}



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,5);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-293,0,500, new TGeoRotation("r1",90,0,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,5);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(293,0,500, new TGeoRotation("r1",90,0,0)));











//  1 ___||____||_____||____||____||____||____||


k=0; i=0;

while (i<5){
   while(k<25){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 1.5,5,15);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*25),78,450, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<25){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 1.5,5,15);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*25),-78,450, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<7){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,1.5,15);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,-78+(k*25),450, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while (k<7){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,1.5,15);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,-78+(k*25),450, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}




//  2 ___||____||_____||____||____||____||____||


k=0; i=0;

while (i<5){
while(k<19){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
   mBlock->SetLineColor(50);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),78,525, new TGeoRotation("r1",0,0,0)));
   k++;
}
   i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),-78,525, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,-78+(k*33),525, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,-78+(k*33),525, new TGeoRotation("r1",0,0,0)));
      k++;
   }
      i++; k=0;
}




//        ||//  ||//  ||//  ||//

//down

k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),90,522, new TGeoRotation("r1",0,-45,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<19){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300+(k*33.3),-90,522, new TGeoRotation("r1",0,45,0)));
      k++;
   }
      i++; k=0;
}


k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-318,-78+(k*33.3),525, new TGeoRotation("r1",-90,45,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<5){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 5,5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(318,-78+(k*33.3),525, new TGeoRotation("r1",90,45,0)));
      k++;
   }
      i++; k=0;
}


// up


k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),115,562, new TGeoRotation("r1",0,-115,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),-115,562, new TGeoRotation("r1",0,115,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-340,-98+(k*13),565, new TGeoRotation("r1",-90,115,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(340,-98+(k*13),565, new TGeoRotation("r1",90,115,0)));
      k++;
   }
      i++; k=0;
}


//up2



k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),115,375, new TGeoRotation("r1",0,-115,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),-115,375, new TGeoRotation("r1",0,115,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-340,-98+(k*13),375, new TGeoRotation("r1",-90,115,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(50);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(340,-98+(k*13),375, new TGeoRotation("r1",90,115,0)));
      k++;
   }
      i++; k=0;
}


//up 3

k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),115,568, new TGeoRotation("r1",0,-115,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),-115,568, new TGeoRotation("r1",0,115,0)));
      k++;
   }
      i++; k=0;
}


k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-340,-98+(k*13),568, new TGeoRotation("r1",-90,115,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(340,-98+(k*13),568, new TGeoRotation("r1",90,115,0)));
      k++;
   }
      i++; k=0;
}





//up4


k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),115,385, new TGeoRotation("r1",0,-115,0)));
      k++;
   }
      i++; k=0;
}




k=0; i=0;

while (i<5){
   while(k<50){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-320+(k*13),-115,385, new TGeoRotation("r1",0,115,0)));
      k++;
   }
      i++; k=0;
}


k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-340,-98+(k*13),385, new TGeoRotation("r1",-90,115,0)));
      k++;
   }
      i++; k=0;
}





k=0; i=0;

while (i<5){
   while(k<17){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeBox(nBlocks,Iron, 2.5,2.5,20);
      mBlock->SetLineColor(44);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(340,-98+(k*13),385, new TGeoRotation("r1",90,115,0)));
      k++;
   }
      i++; k=0;
}


// up kiwa
   //=========
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 270,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,0,620, new TGeoRotation("r1",0,0,0)));
   //===============//2
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 75,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,-50,600, new TGeoRotation("r1",0,20,-40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 75,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,50,600, new TGeoRotation("r1",0,-20,40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 75,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,50,600, new TGeoRotation("r1",0,-20,-40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 75,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,-50,600, new TGeoRotation("r1",0,20,40)));




   //===============//1
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 50,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,-80,413, new TGeoRotation("r1",0,20,-40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 50,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(300,80,413, new TGeoRotation("r1",0,-20,40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 50,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,80,413, new TGeoRotation("r1",0,-20,-40)));

   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron, 50,15,20);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-300,-80,413, new TGeoRotation("r1",0,20,40)));




// _1_

//front

k=0; i=0;
while (i<7){
   while(k<44){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-280+(k*13),70+(i*12.5),425-(i*5), new TGeoRotation("r1",0,60,0)));
      k++;
   }
      i++; k=0;
}



k=0; i=0;

while (i<7){
   while(k<44){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-280+(k*13),-70-(i*12.5),425-(i*5), new TGeoRotation("r1",0,120,0)));
      k++;
   }
      i++; k=0;
}

//_2_




k=0; i=0;
while (i<11){
   while(k<43){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);

      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*13),15+(i*12.5),620-(i*5), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}



k=0; i=0;

while (i<11){
   while(k<43){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*13),-15-(i*12.5),620-(i*5), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}




//////left
k=0; i=0;

while (i<6){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-335,81.25+(i*12.5),592.5-(i*2), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<7){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-322,69.75+(i*12.5),595-(i*2), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<8){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-309,56.25+(i*12.5),605-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<9){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-296,50+(i*12.5),610-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<10){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-283,37.5+(i*12.5),615-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}



k=0; i=0;

while (i<6){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-335,-81.25-(i*12.5),592.5-(i*2), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<7){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-322,-69.75-(i*12.5),595-(i*2), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<8){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-309,-56.25-(i*12.5),605-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<9){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-296,-50-(i*12.5),610-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<10){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-283,-37.5-(i*12.5),615-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}



//////////right



k=0; i=0;

while (i<6){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(335,81.25+(i*12.5),592.5-(i*2), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<7){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(322,69.75+(i*12.5),595-(i*2), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<8){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(309,56.25+(i*12.5),605-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<9){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(296,50+(i*12.5),610-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<10){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(283,37.5+(i*12.5),615-(i*4), new TGeoRotation("r1",0,60,0)));
      k++;
   }
   i++; k=0;
}


//





k=0; i=0;

while (i<6){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(335,-81.25-(i*12.5),592.5-(i*2), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<7){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(322,-69.75-(i*12.5),595-(i*2), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<8){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(309,-56.25-(i*12.5),605-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}

k=0; i=0;

while (i<9){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(296,-50-(i*12.5),610-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}


k=0; i=0;

while (i<10){
   while(k<11){
      sprintf(nBlocks,"ab%d",N++);
      mBlock = geom->MakeTubs(nBlocks,Iron, 3,6,6,10,170);
      mBlock->SetLineColor(13);
      top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(283,-37.5-(i*12.5),615-(i*4), new TGeoRotation("r1",0,120,0)));
      k++;
   }
   i++; k=0;
}


//   /// || / / / / / / / || / / / / / / / / || / / / / / / / / / / /


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 330,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,-110,550, new TGeoRotation("r1",0,-45,0)));




   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 330,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,110,550, new TGeoRotation("r1",0,45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 110,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(335,0,550, new TGeoRotation("r1",90,-45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 110,10,2);
   mBlock->SetLineColor(42);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-335,0,550, new TGeoRotation("r1",90,45,0)));



////////////////////////////////add box





k=0; i=0;

while (i<5){
while(k<6){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron,18,10,2);
   mBlock->SetLineColor(8);
   {
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*100),-111,550, new TGeoRotation("r1",0,-45,0)));

   }
   k++;
}
   i++; k=0;
}




k=0; i=0;

while (i<5){
while(k<6){
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks,Iron,18,10,2);
   mBlock->SetLineColor(8);
   {
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-270+(k*100),111,550, new TGeoRotation("r1",0,45,0)));

   }
   k++;
}
   i++; k=0;
}


   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,2);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(336,0,550, new TGeoRotation("r1",90,-45,0)));



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 18,10,2);
   mBlock->SetLineColor(8);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(-336,0,550, new TGeoRotation("r1",90,45,0)));




//                  |           |           |            |           |



   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 300,75,40);
   mBlock->SetLineColor(45);
   top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(0,0,450, new TGeoRotation("r1",0,0,0)));



//kiwa
   sprintf(nBlocks,"ab%d",N++);
   mBlock = geom->MakeBox(nBlocks, Iron, 305,80,2.5);
   mBlock->SetLineColor(10);
   top->AddNodeOverlap(mBlock,1,new TGeoTranslation(0,0,430));




  top->SetVisibility(0);
  geom->CloseGeometry();

  top->Draw("ogl");
}
