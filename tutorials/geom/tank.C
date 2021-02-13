/// \file
/// \ingroup tutorial_geom
/// Drawing a fine tank, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x tank.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_tank.png width=800px
/// \macro_code
///
/// \author Dong Gyu Lee (ravirus@hanmail.net), Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void tank()
{
   TGeoManager *geom = new TGeoManager("geom","My 3D Project");

//------------------Creat materials-----------------------------
   TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
   TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
   TGeoMaterial *Cu = new TGeoMaterial("Cu",63.549,29,8.92);

//------------------Creat media----------------------------------
   TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
   TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);

//------------------Create TOP volume----------------------------
   TGeoVolume *top = geom->MakeBox("top",Air,100,100,100);
   geom->SetTopVolume(top);
   geom->SetTopVisible(0);
   // If you want to see the boundary, please input the number, 1 instead of 0.
   // Like this, geom->SetTopVisible(1);


//-----------------Create Object volume--------------------------


//Now, we start real shape

//UpperBody
   TGeoVolume *pl=geom->MakeBox("pl",Iron,210,93,20);
   pl->SetLineColor(42);
   TGeoVolume *pl1=geom->MakeBox("pl1",Iron,217,50,5);
   pl1->SetLineColor(42);
   TGeoVolume *pl2=geom->MakeTrd2("pl2",Iron,219,150,50,40,10);
   pl2->SetLineColor(42);
   TGeoVolume *plu=geom->MakeTrd2("plu",Iron,210,70,100,100,5);
   plu->SetLineColor(42);
   top->AddNodeOverlap(plu,1,new TGeoTranslation(0,0,-105));
   TGeoVolume *sp=geom->MakeTubs("sp",Iron,30,40,50,10,60);//Small Plate front
   sp->SetLineColor(42);

//Top which will have the gun
   TGeoVolume *tp=geom->MakeSphere("tp",Iron,0,100,67,90,0,360);//tp is Top with gun
   tp->SetLineColor(12);
   TGeoVolume *tp1=geom->MakeSphere("tp1",Iron,90,190,0,29,0,360);//tp1 is Top with roof
   tp1->SetLineColor(12);
   TGeoVolume *mgg=geom->MakeTubs("mgg",Iron,0,25,30,42,136);//Main Gun Guard
   mgg->SetLineColor(12);
   TGeoVolume *mgg1=geom->MakeTrd2("mgg1",Iron,30.5,45,19,30,35);
   mgg1->SetLineColor(12);

   top->AddNodeOverlap(mgg1,1,new TGeoCombiTrans(-57,0,-63,new TGeoRotation("mgg",90,90,0)));
   top->AddNodeOverlap(mgg,1,new TGeoCombiTrans(-75,0,-63,new TGeoRotation("mgg",0,90,90)));

//Small Top infront Top
   TGeoVolume *stp=geom->MakeSphere("stp",Iron,0,30,67,90,0,360);//Top for driver
   stp->SetLineColor(12);
   TGeoVolume *stp1=geom->MakeSphere("stp1",Iron,115,120,0,12,0,360);//Top with roof
   stp1->SetLineColor(12);
   TGeoVolume *stpo1=geom->MakeBox("stpo1",Iron,3,1,5);
   stpo1->SetLineColor(42);//Small T P Option 1

   top->AddNodeOverlap(stpo1,1,new TGeoTranslation(-93,-32,-95));
   top->AddNodeOverlap(stpo1,1,new TGeoTranslation(-93,-38,-95));
   top->AddNodeOverlap(stp,1,new TGeoTranslation(-120,-35,-108));
   top->AddNodeOverlap(stp1,1,new TGeoCombiTrans(-185,-35,-168,new TGeoRotation("stp1",90,40,0)));

//The Main Gun1 with AddNodeOverlap
   TGeoVolume *mg1=geom->MakeCone("mg1",Iron,160,4,5,4,7);
   mg1->SetLineColor(12);
   top->AddNodeOverlap(mg1,1,new TGeoCombiTrans(-220,0,-53,new TGeoRotation("bs",90,94,0)));
   TGeoVolume *mg1o1=geom->MakeCone("mg1o1",Iron,40,4.1,8,4.1,8);
   mg1o1->SetLineColor(12);//
   top->AddNodeOverlap(mg1o1,1,new TGeoCombiTrans(-220,0,-53,new TGeoRotation("bs",90,94,0)));


//Underbody
   TGeoVolume *underbody=geom->MakeTrd2("underbody",Iron,160,210,93,93,30);
   underbody->SetLineColor(28);
   TGeoVolume *bs=geom->MakeTubs("bs",Iron,0,20,93,10,270);
   bs->SetLineColor(42);
   TGeoVolume *bsp=geom->MakeTubs("bsp",Iron,0,20,30,10,270);
   bsp->SetLineColor(42);

   TGeoVolume *Tip=geom->MakeCone("Tip",Iron,21,0,24,0,24); //Tip is wheel
   Tip->SetLineColor(12);
   TGeoVolume *Tip1=geom->MakeCone("Tip1",Iron,10,23,30,25,30);
   Tip1->SetLineColor(14);
   TGeoVolume *Tip2=geom->MakeCone("Tip2",Iron,30,0,7,0,7);
   Tip2->SetLineColor(42);

   TGeoVolume *wheel=geom->MakeCone("wheel",Iron,30,0,7,0,7);
   wheel->SetLineColor(42);
   TGeoVolume *wheel1=geom->MakeCone("wheel1",Iron,21,0,16,0,16); //innner wheel
   wheel1->SetLineColor(14);
   TGeoVolume *wheel2=geom->MakeCone("wheel2",Iron,10,15,22,15,22); //outter wheel
   wheel2->SetLineColor(12);

   TGeoVolume *Tip0=geom->MakeCone("Tip0",Iron,30,0,7,0,7);
   Tip0->SetLineColor(12);
   TGeoVolume *Tip01=geom->MakeCone("Tip01",Iron,10,7,10.5,7,10.5);
   Tip0->SetLineColor(14);

//cycle of chain with AddNodeOverlap
   char name[50];
   TGeoVolume *WH;//piece of chain
   TGeoVolume *whp;
   TGeoVolume *who;

   //consist upper chain
   for(int i=0;i<26;i++){
      sprintf(name,"wh%d",i);
      WH = geom->MakeBox(name,Iron,5.5,22,2);
      whp = geom->MakeBox(name,Iron,5,2.1,4);
      who = geom->MakeBox(name,Iron,2,6,1);
      WH->SetLineColor(12);
      whp->SetLineColor(14);
      who->SetLineColor(42);
      top->AddNodeOverlap(WH,1,new TGeoTranslation(-195+(15*i),-120,-125));
      top->AddNodeOverlap(WH,1,new TGeoTranslation(-195+(15*i),120,-125));

      top->AddNodeOverlap(whp,1,new TGeoTranslation(-195+(15*i),-120,-127));
      top->AddNodeOverlap(whp,1,new TGeoTranslation(-195+(15*i),120,-127));

      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195+(15*i),-127,-123, new TGeoRotation("who",-15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195+(15*i),-113,-123, new TGeoRotation("who",15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195+(15*i),127,-123, new TGeoRotation("who",15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195+(15*i),113,-123, new TGeoRotation("who",-15,0,0)));

   }
   //chain connector
   TGeoVolume *WHl = geom->MakeBox(name,Iron,187.5,5,1);
   WHl->SetLineColor(12);
   top->AddNodeOverlap(WHl,1,new TGeoTranslation(-7.5,-129,-125));
   top->AddNodeOverlap(WHl,1,new TGeoTranslation(-7.5,-111,-125));
   top->AddNodeOverlap(WHl,1,new TGeoTranslation(-7.5,111,-125));
   top->AddNodeOverlap(WHl,1,new TGeoTranslation(-7.5,129,-125));

   //just one side
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(34*(3.14/180))),-120,-150+(25*cos(34*(3.14/180))), new TGeoRotation("who",90,34,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(68*(3.14/180))),-120,-150+(25*cos(68*(3.14/180))), new TGeoRotation("who",90,68,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(102*(3.14/180))),-120,-150+(25*cos(102*(3.14/180))), new TGeoRotation("who",90,102,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180))),-120,-150+(25*cos(136*(3.14/180))), new TGeoRotation("who",90,136,-90)));

   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-12,-120,-150+(25*cos(136*(3.14/180)))-10, new TGeoRotation("who",90,140,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-24,-120,-150+(25*cos(136*(3.14/180)))-20, new TGeoRotation("who",90,142,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-37,-120,-150+(25*cos(136*(3.14/180)))-30, new TGeoRotation("who",90,145,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-50,-120,-150+(25*cos(136*(3.14/180)))-40, new TGeoRotation("who",90,149,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(34*(3.14/180))),-120,-150+(22.8*cos(34*(3.14/180))), new TGeoRotation("whp",90,34,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(68*(3.14/180))),-120,-150+(22.8*cos(68*(3.14/180))), new TGeoRotation("whp",90,68,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(102*(3.14/180))),-120,-150+(22.8*cos(102*(3.14/180))), new TGeoRotation("whp",90,102,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180))),-120,-150+(22.8*cos(136*(3.14/180))), new TGeoRotation("whp",90,136,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-12,-120,-150+(22.8*cos(136*(3.14/180)))-10, new TGeoRotation("whp",90,140,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-24,-120,-150+(22.8*cos(136*(3.14/180)))-20, new TGeoRotation("whp",90,142,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-37,-120,-150+(22.8*cos(136*(3.14/180)))-30, new TGeoRotation("whp",90,145,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-50,-120,-150+(22.8*cos(136*(3.14/180)))-40, new TGeoRotation("whp",90,149,-90)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(34*(3.14/180))),-127,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",97.5,34,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(68*(3.14/180))),-127,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",97.5,68,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(102*(3.14/180))),-127,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",97.5,102,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180))),-127,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",97.5,136,-97.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-12,-127,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",97.5,140,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-24,-127,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",97.5,142,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-37,-127,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",97.5,145,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-50,-127,-150+(27*cos(136*(3.14/180)))-40, new TGeoRotation("who",97.5,149,-97.5)));
//--------------------------
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(34*(3.14/180))),-113,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",82.5,34,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(68*(3.14/180))),-113,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",82.5,68,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(102*(3.14/180))),-113,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",82.5,102,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180))),-113,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",82.5,136,-82.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-12,-113,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",82.5,140,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-24,-113,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",82.5,142,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-37,-113,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",82.5,145,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-50,-113,-150+(27*cos(136*(3.14/180)))-40, new TGeoRotation("who",82.5,149,-82.5)));


   TGeoVolume *chc0=geom->MakeTubs("chc0",Iron,24.5,26.5,5,-34,0);//Small Plate front
   chc0->SetLineColor(12);
   TGeoVolume *chc1=geom->MakeTubs("chc1",Iron,24.5,26.5,5,-68,-34);//Small Plate front
   chc1->SetLineColor(12);
   TGeoVolume *chc2=geom->MakeTubs("chc2",Iron,24.5,26.5,5,-102,-68);//Small Plate front
   chc2->SetLineColor(12);
   TGeoVolume *chc3=geom->MakeTubs("chc3",Iron,24.5,26.5,5,-136,-102);//Small Plate front
   chc3->SetLineColor(12);

   top->AddNodeOverlap(chc0,1,new TGeoCombiTrans(180,-129,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1,1,new TGeoCombiTrans(180,-129,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2,1,new TGeoCombiTrans(180,-129,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3,1,new TGeoCombiTrans(180,-129,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chc0,1,new TGeoCombiTrans(180,-111,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1,1,new TGeoCombiTrans(180,-111,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2,1,new TGeoCombiTrans(180,-111,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3,1,new TGeoCombiTrans(180,-111,-150,new TGeoRotation("chc3",0,90,90)));

   TGeoVolume *chcl=geom->MakeBox("chcl",Iron,5,5,1);
   chcl->SetLineColor(12);
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-6,-111,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-18,-111,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-30,-111,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,145,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-43,-111,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,149,-90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-6,-129,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-18,-129,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-30,-129,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,145,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-43,-129,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,149,-90)));

   TGeoVolume *chc4=geom->MakeTubs("chc4",Iron,31.5,34.5,5,-175,-145);//Small Plate front
   chc4->SetLineColor(12);
   top->AddNodeOverlap(chc4,1,new TGeoCombiTrans(130,-111,-180,new TGeoRotation("chc3",0,90,90)));
   top->AddNodeOverlap(chc4,1,new TGeoCombiTrans(130,-129,-180,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(34*(3.14/180))),-120,-150+(25*cos(34*(3.14/180))), new TGeoRotation("who",90,-34,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(68*(3.14/180))),-120,-150+(25*cos(68*(3.14/180))), new TGeoRotation("who",90,-68,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(102*(3.14/180))),-120,-150+(25*cos(102*(3.14/180))), new TGeoRotation("who",90,-102,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180))),-120,-150+(25*cos(136*(3.14/180))), new TGeoRotation("who",90,-136,-90)));

   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+12,-120,-150+(25*cos(136*(3.14/180)))-10, new TGeoRotation("who",90,-140,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+24,-120,-150+(25*cos(136*(3.14/180)))-20, new TGeoRotation("who",90,-142,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+35,-120,-150+(25*cos(136*(3.14/180)))-30, new TGeoRotation("who",90,-139,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+48,-120,-150+(25*cos(136*(3.14/180)))-41, new TGeoRotation("who",90,-153,-90)));


   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(34*(3.14/180))),-120,-150+(22.8*cos(34*(3.14/180))), new TGeoRotation("whp",90,-34,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(68*(3.14/180))),-120,-150+(22.8*cos(68*(3.14/180))), new TGeoRotation("whp",90,-68,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(102*(3.14/180))),-120,-150+(22.8*cos(102*(3.14/180))), new TGeoRotation("whp",90,-102,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180))),-120,-150+(22.8*cos(136*(3.14/180))), new TGeoRotation("whp",90,-136,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+12,-120,-150+(22.8*cos(136*(3.14/180)))-10, new TGeoRotation("whp",90,-140,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+24,-120,-150+(22.8*cos(136*(3.14/180)))-20, new TGeoRotation("whp",90,-142,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+35,-120,-150+(22.8*cos(136*(3.14/180)))-30, new TGeoRotation("whp",90,-139,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+48,-120,-150+(22.8*cos(136*(3.14/180)))-41, new TGeoRotation("whp",90,-153,-90)));


   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(34*(3.14/180))),-127,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",97.5,-34,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(68*(3.14/180))),-127,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",97.5,-68,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(102*(3.14/180))),-127,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",97.5,-102,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180))),-127,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",97.5,-136,-97.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+12,-127,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",97.5,-140,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+24,-127,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",97.5,-142,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+35,-127,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",97.5,-139,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+48,-127,-150+(27*cos(136*(3.14/180)))-41, new TGeoRotation("who",97.5,-153,-97.5)));
//-------------------------
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(34*(3.14/180))),-113,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",82.5,-34,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(68*(3.14/180))),-113,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",82.5,-68,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(102*(3.14/180))),-113,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",82.5,-102,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180))),-113,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",82.5,-136,-82.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+12,-113,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",82.5,-140,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+24,-113,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",82.5,-142,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+35,-113,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",82.5,-139,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+48,-113,-150+(27*cos(136*(3.14/180)))-41, new TGeoRotation("who",82.5,-153,-82.5)));


   TGeoVolume *chc0i=geom->MakeTubs("chc0i",Iron,24.5,26.5,5,0,34);//Small Plate front
   chc0i->SetLineColor(12);
   TGeoVolume *chc1i=geom->MakeTubs("chc1i",Iron,24.5,26.5,5,34,68);//Small Plate front
   chc1i->SetLineColor(12);
   TGeoVolume *chc2i=geom->MakeTubs("chc2i",Iron,24.5,26.5,5,68,102);//Small Plate front
   chc2i->SetLineColor(12);
   TGeoVolume *chc3i=geom->MakeTubs("chc3i",Iron,24.5,26.5,5,102,136);//Small Plate front
   chc3i->SetLineColor(12);

   top->AddNodeOverlap(chc0i,1,new TGeoCombiTrans(-195,-129,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1i,1,new TGeoCombiTrans(-195,-129,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2i,1,new TGeoCombiTrans(-195,-129,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3i,1,new TGeoCombiTrans(-195,-129,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chc0i,1,new TGeoCombiTrans(-195,-111,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1i,1,new TGeoCombiTrans(-195,-111,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2i,1,new TGeoCombiTrans(-195,-111,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3i,1,new TGeoCombiTrans(-195,-111,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+06,-129,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,-140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+18,-129,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,-142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+29,-129,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,-139,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+41,-129,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,-138,-90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+06,-111,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,-140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+18,-111,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,-142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+29,-111,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,-139,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+41,-111,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,-138,-90)));

   TGeoVolume *chc4i=geom->MakeTubs("chc4i",Iron,31.5,33,5,145,175);//Small Plate front
   chc4i->SetLineColor(12);
   top->AddNodeOverlap(chc4i,1,new TGeoCombiTrans(-150,-111,-180,new TGeoRotation("chc3",0,90,90)));
   top->AddNodeOverlap(chc4i,1,new TGeoCombiTrans(-150,-129,-180,new TGeoRotation("chc3",0,90,90)));


//just other side


   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(34*(3.14/180))),120,-150+(25*cos(34*(3.14/180))), new TGeoRotation("who",90,34,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(68*(3.14/180))),120,-150+(25*cos(68*(3.14/180))), new TGeoRotation("who",90,68,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(102*(3.14/180))),120,-150+(25*cos(102*(3.14/180))), new TGeoRotation("who",90,102,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180))),120,-150+(25*cos(136*(3.14/180))), new TGeoRotation("who",90,136,-90)));

   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-12,120,-150+(25*cos(136*(3.14/180)))-10, new TGeoRotation("who",90,140,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-24,120,-150+(25*cos(136*(3.14/180)))-20, new TGeoRotation("who",90,142,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-37,120,-150+(25*cos(136*(3.14/180)))-30, new TGeoRotation("who",90,145,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-50,120,-150+(25*cos(136*(3.14/180)))-40, new TGeoRotation("who",90,149,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(34*(3.14/180))),120,-150+(22.8*cos(34*(3.14/180))), new TGeoRotation("whp",90,34,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(68*(3.14/180))),120,-150+(22.8*cos(68*(3.14/180))), new TGeoRotation("whp",90,68,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(102*(3.14/180))),120,-150+(22.8*cos(102*(3.14/180))), new TGeoRotation("whp",90,102,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180))),120,-150+(22.8*cos(136*(3.14/180))), new TGeoRotation("whp",90,136,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-12,120,-150+(22.8*cos(136*(3.14/180)))-10, new TGeoRotation("whp",90,140,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-24,120,-150+(22.8*cos(136*(3.14/180)))-20, new TGeoRotation("whp",90,142,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-37,120,-150+(22.8*cos(136*(3.14/180)))-30, new TGeoRotation("whp",90,145,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(180+(22.8*sin(136*(3.14/180)))-50,120,-150+(22.8*cos(136*(3.14/180)))-40, new TGeoRotation("whp",90,149,-90)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(34*(3.14/180))),113,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",97.5,34,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(68*(3.14/180))),113,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",97.5,68,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(102*(3.14/180))),113,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",97.5,102,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180))),113,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",97.5,136,-97.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-12,113,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",97.5,140,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-24,113,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",97.5,142,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-37,113,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",97.5,145,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-50,113,-150+(27*cos(136*(3.14/180)))-40, new TGeoRotation("who",97.5,149,-97.5)));
//--------------------------
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(34*(3.14/180))),127,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",82.5,34,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(68*(3.14/180))),127,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",82.5,68,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(102*(3.14/180))),127,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",82.5,102,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180))),127,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",82.5,136,-82.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-12,127,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",82.5,140,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-24,127,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",82.5,142,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-37,127,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",82.5,145,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(180+(27*sin(136*(3.14/180)))-50,127,-150+(27*cos(136*(3.14/180)))-40, new TGeoRotation("who",82.5,149,-82.5)));


   top->AddNodeOverlap(chc0,1,new TGeoCombiTrans(180,129,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1,1,new TGeoCombiTrans(180,129,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2,1,new TGeoCombiTrans(180,129,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3,1,new TGeoCombiTrans(180,129,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chc0,1,new TGeoCombiTrans(180,111,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1,1,new TGeoCombiTrans(180,111,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2,1,new TGeoCombiTrans(180,111,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3,1,new TGeoCombiTrans(180,111,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-6,111,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-18,111,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-30,111,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,145,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-43,111,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,149,-90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-6,129,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-18,129,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-30,129,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,145,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(180+(25*sin(136*(3.14/180)))-43,129,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,149,-90)));


   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(34*(3.14/180))),120,-150+(25*cos(34*(3.14/180))), new TGeoRotation("who",90,-34,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(68*(3.14/180))),120,-150+(25*cos(68*(3.14/180))), new TGeoRotation("who",90,-68,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(102*(3.14/180))),120,-150+(25*cos(102*(3.14/180))), new TGeoRotation("who",90,-102,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180))),120,-150+(25*cos(136*(3.14/180))), new TGeoRotation("who",90,-136,-90)));

   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+12,120,-150+(25*cos(136*(3.14/180)))-10, new TGeoRotation("who",90,-140,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+24,120,-150+(25*cos(136*(3.14/180)))-20, new TGeoRotation("who",90,-142,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+35,120,-150+(25*cos(136*(3.14/180)))-30, new TGeoRotation("who",90,-139,-90)));
   top->AddNodeOverlap(WH,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+48,120,-150+(25*cos(136*(3.14/180)))-41, new TGeoRotation("who",90,-153,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(34*(3.14/180))),120,-150+(22.8*cos(34*(3.14/180))), new TGeoRotation("whp",90,-34,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(68*(3.14/180))),120,-150+(22.8*cos(68*(3.14/180))), new TGeoRotation("whp",90,-68,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(102*(3.14/180))),120,-150+(22.8*cos(102*(3.14/180))), new TGeoRotation("whp",90,-102,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180))),120,-150+(22.8*cos(136*(3.14/180))), new TGeoRotation("whp",90,-136,-90)));

   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+12,120,-150+(22.8*cos(136*(3.14/180)))-10, new TGeoRotation("whp",90,-140,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+24,120,-150+(22.8*cos(136*(3.14/180)))-20, new TGeoRotation("whp",90,-142,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+35,120,-150+(22.8*cos(136*(3.14/180)))-30, new TGeoRotation("whp",90,-139,-90)));
   top->AddNodeOverlap(whp,1,new TGeoCombiTrans(-195-(22.8*sin(136*(3.14/180)))+48,120,-150+(22.8*cos(136*(3.14/180)))-41, new TGeoRotation("whp",90,-153,-90)));


   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(34*(3.14/180))),113,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",97.5,-34,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(68*(3.14/180))),113,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",97.5,-68,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(102*(3.14/180))),113,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",97.5,-102,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180))),113,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",97.5,-136,-97.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+12,113,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",97.5,-140,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+24,113,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",97.5,-142,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+35,113,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",97.5,-139,-97.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+48,113,-150+(27*cos(136*(3.14/180)))-41, new TGeoRotation("who",97.5,-153,-97.5)));
//-------------------------
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(34*(3.14/180))),127,-150+(27*cos(34*(3.14/180))), new TGeoRotation("who",82.5,-34,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(68*(3.14/180))),127,-150+(27*cos(68*(3.14/180))), new TGeoRotation("who",82.5,-68,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(102*(3.14/180))),127,-150+(27*cos(102*(3.14/180))), new TGeoRotation("who",82.5,-102,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180))),127,-150+(27*cos(136*(3.14/180))), new TGeoRotation("who",82.5,-136,-82.5)));

   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+12,127,-150+(27*cos(136*(3.14/180)))-10, new TGeoRotation("who",82.5,-140,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+24,127,-150+(27*cos(136*(3.14/180)))-20, new TGeoRotation("who",82.5,-142,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+35,127,-150+(27*cos(136*(3.14/180)))-30, new TGeoRotation("who",82.5,-139,-82.5)));
   top->AddNodeOverlap(who,1,new TGeoCombiTrans(-195-(27*sin(136*(3.14/180)))+48,127,-150+(27*cos(136*(3.14/180)))-41, new TGeoRotation("who",82.5,-153,-82.5)));


   top->AddNodeOverlap(chc0i,1,new TGeoCombiTrans(-195,129,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1i,1,new TGeoCombiTrans(-195,129,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2i,1,new TGeoCombiTrans(-195,129,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3i,1,new TGeoCombiTrans(-195,129,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chc0i,1,new TGeoCombiTrans(-195,111,-150,new TGeoRotation("chc0",0,90,90)));
   top->AddNodeOverlap(chc1i,1,new TGeoCombiTrans(-195,111,-150,new TGeoRotation("chc1",0,90,90)));
   top->AddNodeOverlap(chc2i,1,new TGeoCombiTrans(-195,111,-150,new TGeoRotation("chc2",0,90,90)));
   top->AddNodeOverlap(chc3i,1,new TGeoCombiTrans(-195,111,-150,new TGeoRotation("chc3",0,90,90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+06,129,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,-140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+18,129,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,-142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+29,129,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,-139,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+41,129,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,-138,-90)));

   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+06,111,-150+(25*cos(136*(3.14/180)))-5, new TGeoRotation("chcl",90,-140,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+18,111,-150+(25*cos(136*(3.14/180)))-15, new TGeoRotation("chcl",90,-142,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+29,111,-150+(25*cos(136*(3.14/180)))-25, new TGeoRotation("chcl",90,-139,-90)));
   top->AddNodeOverlap(chcl,1,new TGeoCombiTrans(-195-(25*sin(136*(3.14/180)))+41,111,-150+(25*cos(136*(3.14/180)))-35, new TGeoRotation("chcl",90,-138,-90)));

   //consist under chain
   for(int i=0;i<20;i++){
      sprintf(name,"wh%d",i);
      top->AddNodeOverlap(WH,1,new TGeoTranslation(-150+(15*i),-120,-212));
      top->AddNodeOverlap(WH,1,new TGeoTranslation(-150+(15*i),120,-212));

      top->AddNodeOverlap(whp,1,new TGeoTranslation(-150+(15*i),-120,-210));
      top->AddNodeOverlap(whp,1,new TGeoTranslation(-150+(15*i),120,-210));

      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-150+(15*i),-127,-214, new TGeoRotation("who",15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-150+(15*i),-113,-214, new TGeoRotation("who",-15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-150+(15*i),127,-214, new TGeoRotation("who",-15,0,0)));
      top->AddNodeOverlap(who,1,new TGeoCombiTrans(-150+(15*i),113,-214, new TGeoRotation("who",15,0,0)));
   }
   TGeoVolume *WHlu = geom->MakeBox(name,Iron,140,5,1);//chain connetor in under
   WHlu->SetLineColor(12);
   top->AddNodeOverlap(WHlu,1,new TGeoTranslation(-7.5,-129,-212));
   top->AddNodeOverlap(WHlu,1,new TGeoTranslation(-7.5,-111,-212));
   top->AddNodeOverlap(WHlu,1,new TGeoTranslation(-7.5,129,-212));
   top->AddNodeOverlap(WHlu,1,new TGeoTranslation(-7.5,111,-212));




//Now, we put real shape

   top->AddNodeOverlap(underbody,1,new TGeoTranslation(0,0,-160));
   top->AddNodeOverlap(pl,1,new TGeoTranslation(0,0,-130));
   top->AddNodeOverlap(tp,1,new TGeoTranslation(30,0,-83));
   top->AddNodeOverlap(tp1,1,new TGeoTranslation(30,0,-208));
   top->AddNodeOverlap(pl2,1,new TGeoTranslation(0,-120,-100));
   top->AddNodeOverlap(pl2,1,new TGeoTranslation(0,120,-100));
   top->AddNodeOverlap(pl1,1,new TGeoTranslation(0,-120,-115));
   top->AddNodeOverlap(pl1,1,new TGeoTranslation(0,120,-115));
   top->AddNodeOverlap(bs,1,new TGeoCombiTrans(180,0,-150,new TGeoRotation("bs",180,90,90)));
   top->AddNodeOverlap(bsp,1,new TGeoCombiTrans(-195,61.5,-150,new TGeoRotation("bsp",0,90,90)));
   top->AddNodeOverlap(bsp,1,new TGeoCombiTrans(-195,-61.5,-150,new TGeoRotation("bsp",0,90,90)));


   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-115,-132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-45,-132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(35,-132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(95,-132.5,-140,new TGeoRotation("Tip01",0,90,90)));

   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-115,-107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-45,-107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(35,-107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(95,-107.5,-140,new TGeoRotation("Tip01",0,90,90)));

   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(-115,-110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(-45,-110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(35,-110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(95,-110.5,-140,new TGeoRotation("Tip0",0,90,90)));

   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-150,-120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-80,-120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-10,-120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(60,-120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(130,-120,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-150,-107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-150,-132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-80,-107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-80,-132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-10,-107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-10,-132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(60,-107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(60,-132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(130,-107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(130,-132.5,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-150,-112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-80,-112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-10,-112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(60,-112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(130,-112.5,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(wheel1,1,new TGeoCombiTrans(180,-120,-150,new TGeoRotation("wheel1",0,90,90)));
   top->AddNodeOverlap(wheel1,1,new TGeoCombiTrans(-195,-120,-150,new TGeoRotation("wheel1",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(180,-107.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(180,-132.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(-195,-107.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(-195,-132.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel,1,new TGeoCombiTrans(180,-112.5,-150,new TGeoRotation("wheel",0,90,90)));
   top->AddNodeOverlap(wheel,1,new TGeoCombiTrans(-195,-112.5,-150,new TGeoRotation("wheel2",0,90,90)));

   top->AddNodeOverlap(sp,1,new TGeoCombiTrans(-209,-120,-149,new TGeoRotation("sp",0,90,90)));//sp!
   top->AddNodeOverlap(sp,1,new TGeoCombiTrans(209,-120,-149,new TGeoRotation("sp1",180,90,90)));//sp!

   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-115,132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-45,132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(35,132.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(95,132.5,-140,new TGeoRotation("Tip01",0,90,90)));

   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-115,107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(-45,107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(35,107.5,-140,new TGeoRotation("Tip01",0,90,90)));
   top->AddNodeOverlap(Tip01,1,new TGeoCombiTrans(95,107.5,-140,new TGeoRotation("Tip01",0,90,90)));

   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(-115,110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(-45,110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(35,110.5,-140,new TGeoRotation("Tip0",0,90,90)));
   top->AddNodeOverlap(Tip0,1,new TGeoCombiTrans(95,110.5,-140,new TGeoRotation("Tip0",0,90,90)));

   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-150,120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-80,120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(-10,120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(60,120,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip,1,new TGeoCombiTrans(130,120,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-150,107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-150,132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-80,107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-80,132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-10,107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(-10,132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(60,107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(60,132.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(130,107.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip1,1,new TGeoCombiTrans(130,132.5,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-150,112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-80,112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(-10,112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(60,112.5,-180,new TGeoRotation("Tip",0,90,90)));
   top->AddNodeOverlap(Tip2,1,new TGeoCombiTrans(130,112.5,-180,new TGeoRotation("Tip",0,90,90)));

   top->AddNodeOverlap(wheel,1,new TGeoCombiTrans(-195,112.5,-150,new TGeoRotation("wheel1",0,90,90)));
   top->AddNodeOverlap(wheel,1,new TGeoCombiTrans(180,112.5,-150,new TGeoRotation("wheel",0,90,90)));
   top->AddNodeOverlap(wheel1,1,new TGeoCombiTrans(180,120,-150,new TGeoRotation("wheel1",0,90,90)));
   top->AddNodeOverlap(wheel1,1,new TGeoCombiTrans(-195,120,-150,new TGeoRotation("wheel1",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(180,107.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(180,132.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(-195,107.5,-150,new TGeoRotation("wheel2",0,90,90)));
   top->AddNodeOverlap(wheel2,1,new TGeoCombiTrans(-195,132.5,-150,new TGeoRotation("wheel2",0,90,90)));

   top->AddNodeOverlap(sp,1,new TGeoCombiTrans(-209,120,-149,new TGeoRotation("sp",0,90,90)));//sp!
   top->AddNodeOverlap(sp,1,new TGeoCombiTrans(209,120,-149,new TGeoRotation("sp1",180,90,90)));//sp!
   top->SetVisibility(0);
   geom->CloseGeometry();


//------------------draw on GL viewer-------------------------------
   top->Draw("ogl");

}
