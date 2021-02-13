/// \file
/// \ingroup tutorial_geom
/// Drawing a famous Korean robot, TaekwonV, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x robot.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_robot.png width=500px
/// \macro_code
///
/// \author Jin Hui Hwang, Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void robot()
{


   TGeoManager *Robot = new TGeoManager("Robot","This is Taegwon V");

   TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
   TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);

   TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
   TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);

   // create volume

   TGeoVolume *top=Robot->MakeBox("top",Air,1000,1000,1000);
   Robot->SetTopVolume(top);
   Robot->SetTopVisible(0);
      // If you want to see the boundary, please input the number, 1 instead of 0.
      // Like this, geom->SetTopVisible(1);



   // head
   TGeoVolume *Band=Robot->MakeEltu("Band",Iron,20,20,2.5);
     Band->SetLineColor(12);
     Band->SetFillColor(12);
   TGeoVolume *Band_b=Robot->MakeSphere("Band_b",Iron,0,2,0,180,180,360);
     Band_b->SetLineColor(2);
     Band_b->SetFillColor(2);
   TGeoVolume *Head=Robot->MakeSphere("Head",Iron,0,19,0,180,180,360);
     Head->SetLineColor(17);
     Head->SetFillColor(17);
   TGeoVolume *Horn=Robot->MakeSphere("Horn",Iron,0,10,60,180,240,300);

   // drawing head
   top->AddNodeOverlap(Band,1,new TGeoTranslation(0,0,90));
   float Phi = 3.14;
   int N = 10;

   for (int i=0; i<=N;i++){
      top->AddNodeOverlap(Band_b,1,new TGeoCombiTrans(sin(2*Phi/N*i)*19,-cos(2*Phi/N*i)*19,90,
         new TGeoRotation("R1",-90+(360/N*i),-90,90)));
   }
   top->AddNodeOverlap(Head,1,new TGeoCombiTrans(0,0,87.5,new TGeoRotation("R2",0,-90,0)));

   char name[50];
   float pcs = 30;
   for (int i=1; i<pcs; i++){
      sprintf(name,"Horn%d",i);
      Horn=Robot->MakeSphere(name,Iron,
         10- 10/pcs*i ,10,180-(120/pcs)*i,180-((120/pcs) * (i-1)),240,300);
      Horn->SetLineColor(2);
      Horn->SetFillColor(2);
      top->AddNodeOverlap(Horn,1,new TGeoCombiTrans(0,8,102,new TGeoRotation("R2",0,140,0)));
      top->AddNodeOverlap(Horn,1,new TGeoCombiTrans(0,-8,102,new TGeoRotation("R2",180,140,0)));
   }

   // face
   TGeoVolume *Migan=Robot->MakeGtra("Migan",Iron,3,0,0,0,3,2,11,0,3,3,11,0);
     Migan->SetLineColor(17);
     Migan->SetFillColor(17);
   TGeoVolume *Ko=Robot->MakeGtra("Ko",Iron,7,0,0,0,3,1,5,0,3,2,5,0);
     Ko->SetLineColor(17);
     Ko->SetFillColor(17);
   TGeoVolume *Ko_m=Robot->MakeBox("Ko_m",Iron,2,8,4);
     Ko_m->SetLineColor(17);
     Ko_m->SetFillColor(17);
   TGeoVolume *Bol_1=Robot->MakeBox("Bol_1",Iron,7,5.5,7);
     Bol_1->SetLineColor(17);
     Bol_1->SetFillColor(17);
   TGeoVolume *Bol_2=Robot->MakeGtra("Bol_2",Iron,1,0,0,0,7,0,9,0,7,0,9,0);
     Bol_2->SetLineColor(17);
     Bol_2->SetFillColor(17);
   TGeoVolume *Noon=Robot->MakeBox("Noon",Iron,1,10,5);
     Noon->SetLineColor(12);
     Noon->SetFillColor(12);
   TGeoVolume *Tuck=Robot->MakeBox("Tuck",Iron,2,10,5.5);
     Tuck->SetLineColor(2);
     Tuck->SetFillColor(2);
   TGeoVolume *Tuck_1=Robot->MakeBox("Tuck_1",Iron,2,9,1);
     Tuck_1->SetLineColor(2);
     Tuck_1->SetFillColor(2);
   TGeoVolume *Tuck_2=Robot->MakeBox("Tuck_2",Iron,3,1,14);
     Tuck_2->SetLineColor(2);
     Tuck_2->SetFillColor(2);
   TGeoVolume *Tuck_j=Robot->MakeSphere("Tuck_j",Iron,0,3.5,0,180,0,360);
     Tuck_j->SetLineColor(5);
     Tuck_j->SetFillColor(5);
          TGeoVolume *Ear=Robot->MakeCons("Ear",Iron,1,0,3,0,3,0,360);
     Ear->SetLineColor(5);
     Ear->SetFillColor(5);
          TGeoVolume *Ear_2=Robot->MakeCone("Ear_2",Iron,5,0,0,0,3);
     Ear_2->SetLineColor(5);
     Ear_2->SetFillColor(5);

   // drawing face
   top->AddNodeOverlap(Migan,1,new TGeoCombiTrans(-15,0,88,new TGeoRotation("R2",-90,40,0)));
   top->AddNodeOverlap(Ko,1,new TGeoCombiTrans(-15,0,76.5,new TGeoRotation("R2",-90,-20,0)));
   top->AddNodeOverlap(Ko_m,1,new TGeoTranslation(-9,0,68));
   top->AddNodeOverlap(Bol_1,1,new TGeoCombiTrans(-7,2,76,new TGeoRotation("R2",-30,-10,0)));
   top->AddNodeOverlap(Bol_1,1,new TGeoCombiTrans(-7,-2,76,new TGeoRotation("R2",30,10,0)));
   top->AddNodeOverlap(Bol_2,1,new TGeoCombiTrans(-6.5,-10.5,76,new TGeoRotation("R2",-15,-90,-30)));
   top->AddNodeOverlap(Bol_2,1,new TGeoCombiTrans(-4,-12.5,82.5,new TGeoRotation("R2",-20,-90,-95)));
   top->AddNodeOverlap(Bol_2,1,new TGeoCombiTrans(-7.5,10.5,76,new TGeoRotation("R2",20,-90,-30)));
   top->AddNodeOverlap(Bol_2,1,new TGeoCombiTrans(-4,12.5,82.5,new TGeoRotation("R2",20,-90,-95)));
   top->AddNodeOverlap(Noon,1,new TGeoCombiTrans(-5,-7,86,new TGeoRotation("R2",60,0,0)));
   top->AddNodeOverlap(Noon,1,new TGeoCombiTrans(-5,7,86,new TGeoRotation("R2",-60,0,0)));
   top->AddNodeOverlap(Tuck,1,new TGeoTranslation(-12,0,62.5));
   for (int i=0; i<10;i++) {
      top->AddNodeOverlap(Tuck_1,1,new TGeoCombiTrans(-4.2,11,61+i,new TGeoRotation("R2",90,-20,20)));
      top->AddNodeOverlap(Tuck_1,1,new TGeoCombiTrans(-4.2,-11,61+i,new TGeoRotation("R2",90,-20,-20)));
   }
   top->AddNodeOverlap(Tuck_2,1,new TGeoTranslation(2,-15.1,76));
   top->AddNodeOverlap(Tuck_2,1,new TGeoTranslation(2,15.1,76));
   top->AddNodeOverlap(Tuck_j,1,new TGeoTranslation(-13,0,62.5));
   top->AddNodeOverlap(Ear,1,new TGeoCombiTrans(2,-16.5,80,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Ear,1,new TGeoCombiTrans(2,16.5,80,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Ear_2,1,new TGeoCombiTrans(2,-20,80,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Ear_2,1,new TGeoCombiTrans(2,20,80,new TGeoRotation("R2",0,90,0)));


   for (int i=1; i<28; i+=1) {
      float a=i*0.2;
             TGeoVolume *Hear=Robot->MakeCons("Hear",Iron,3,13+a,16+a,13+a,16+a,-60-a,60+a);
      if (i<27) {
        Hear->SetLineColor(12);
        Hear->SetFillColor(12);
      }
      else {
        Hear->SetLineColor(2);
        Hear->SetFillColor(2);
      }
      top->AddNodeOverlap(Hear,1,new TGeoTranslation(0,0,89-i));
   }
   for (int i=1; i<28; i+=1) {
      float a=i*0.2;
             TGeoVolume *Hear=Robot->MakeCons("Hear",Iron,3,13+a,16+a,13+a,16+a,-70-a,-60-a);
        Hear->SetLineColor(2);
        Hear->SetFillColor(2);
      top->AddNodeOverlap(Hear,1,new TGeoTranslation(0,0,89-i));
   }
   for (int i=1; i<28; i+=1) {
      float a=i*0.2;
             TGeoVolume *Hear=Robot->MakeCons("Hear",Iron,3,13+a,16+a,13+a,16+a,60+a,70+a);
        Hear->SetLineColor(2);
        Hear->SetFillColor(2);
      top->AddNodeOverlap(Hear,1,new TGeoTranslation(0,0,89-i));
   }

   // neck
   TGeoVolume *Mock=Robot->MakeTrd2("Mock",Iron,1,1,7,6.5,20);
     Mock->SetLineColor(17);
     Mock->SetFillColor(17);
   TGeoVolume *Mock_1=Robot->MakeTrd2("Mock_1",Iron,1,1,6,5,20);
     Mock_1->SetLineColor(17);
     Mock_1->SetFillColor(17);
   TGeoVolume *Mock_s=Robot->MakeTrd2("Mock_s",Iron,1,1,5,4.5,20);
     Mock_s->SetLineColor(17);
     Mock_s->SetFillColor(17);

   // drawing neck
          top->AddNodeOverlap(Mock,1,new TGeoCombiTrans(-5,4.7,50,new TGeoRotation("R2",-30,0,-10)));
          top->AddNodeOverlap(Mock,1,new TGeoCombiTrans(-5,-4.7,50,new TGeoRotation("R2",30,0,10)));
          top->AddNodeOverlap(Mock_1,1,new TGeoCombiTrans(11,-4,50,new TGeoRotation("R2",130,-8,10)));
          top->AddNodeOverlap(Mock_1,1,new TGeoCombiTrans(11,4,50,new TGeoRotation("R2",-130,8,-10)));
          top->AddNodeOverlap(Mock_s,1,new TGeoCombiTrans(2.5,9,50,new TGeoRotation("R2",90,0,0)));
          top->AddNodeOverlap(Mock_s,1,new TGeoCombiTrans(2.5,-9,50,new TGeoRotation("R2",90,0,0)));


   // chest
   TGeoVolume *Gasem=Robot->MakeBox("Gasem",Iron,16,50,20);
     Gasem->SetLineColor(12);
     Gasem->SetFillColor(12);
   TGeoVolume *Gasem_b1=Robot->MakeSphere("Gasem_b1",Iron,0,15,0,180,0,360);
     Gasem_b1->SetLineColor(12);
     Gasem_b1->SetFillColor(12);
   TGeoVolume *Gasem_b2=Robot->MakeSphere("Gasem_b2",Iron,0,13,0,180,0,360);
     Gasem_b2->SetLineColor(12);
     Gasem_b2->SetFillColor(12);
   TGeoVolume *Gasem_1=Robot->MakeEltu("Gasem_1",Iron,13,13,20);
     Gasem_1->SetLineColor(12);
     Gasem_1->SetFillColor(12);
   TGeoVolume *Gasem_2=Robot->MakeEltu("Gasem_2",Iron,13,13,19);
     Gasem_2->SetLineColor(12);
     Gasem_2->SetFillColor(12);
          TGeoVolume *Gasem_3=Robot->MakeCone("Gasem_3",Iron,19,0,13,0,15);
     Gasem_3->SetLineColor(12);
     Gasem_3->SetFillColor(12);
   TGeoVolume *Gasem_4=Robot->MakeEltu("Gasem_4",Iron,15,15,16);
     Gasem_4->SetLineColor(12);
     Gasem_4->SetFillColor(12);
   TGeoVolume *Gasem_5=Robot->MakeEltu("Gasem_5",Iron,13,13,16);
     Gasem_5->SetLineColor(12);
     Gasem_5->SetFillColor(12);
   TGeoVolume *Gasem_m1=Robot->MakeBox("Gasem_m1",Iron,19,19,5);
     Gasem_m1->SetLineColor(12);
     Gasem_m1->SetFillColor(12);
   TGeoVolume *Gasem_m2=Robot->MakeTrd2("Gasem_m2",Iron,13,15,2,2,19);
     Gasem_m2->SetLineColor(12);
     Gasem_m2->SetFillColor(12);
   TGeoVolume *V=Robot->MakeTrd2("V",Iron,2,2,22,30,4);
     V->SetLineColor(2);
     V->SetFillColor(2);
   TGeoVolume *V_m=Robot->MakeBox("V_m",Iron,2,7,1);
     V_m->SetLineColor(2);
     V_m->SetFillColor(2);

   // drawing chest
   top->AddNodeOverlap(Gasem,1,new TGeoTranslation(4,0,19));
   top->AddNodeOverlap(Gasem_b1,1,new TGeoTranslation(-12,50,35));
   top->AddNodeOverlap(Gasem_b1,1,new TGeoTranslation(-12,-50,35));
   top->AddNodeOverlap(Gasem_b1,1,new TGeoTranslation(20,50,35));
   top->AddNodeOverlap(Gasem_b1,1,new TGeoTranslation(20,-50,35));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,50,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,-50,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,50,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,-50,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,10,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,-10,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,10,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,-10,-5));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,10,35));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(20,-10,35));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,10,35));
   top->AddNodeOverlap(Gasem_b2,1,new TGeoTranslation(-12,-10,35));
   top->AddNodeOverlap(Gasem_1,1,new TGeoCombiTrans(20,31,-5,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_1,1,new TGeoCombiTrans(20,-31,-5,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_1,1,new TGeoCombiTrans(-12,31,-5,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_1,1,new TGeoCombiTrans(-12,-31,-5,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_2,1,new TGeoCombiTrans(20,10,13,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_2,1,new TGeoCombiTrans(20,-10,13,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_2,1,new TGeoCombiTrans(-12,10,13,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_2,1,new TGeoCombiTrans(-12,-10,13,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(-12,50,16,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(-12,-50,16,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(20,50,16,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(20,-50,16,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(-12,31,35,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(-12,-31,35,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(20,31,35,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Gasem_3,1,new TGeoCombiTrans(20,-31,35,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_4,1,new TGeoCombiTrans(4,-50,35,new TGeoRotation("R2",90,90,0)));
   top->AddNodeOverlap(Gasem_4,1,new TGeoCombiTrans(4,50,35,new TGeoRotation("R2",90,90,0)));
   top->AddNodeOverlap(Gasem_5,1,new TGeoCombiTrans(4,-50,-5,new TGeoRotation("R2",90,90,0)));
   top->AddNodeOverlap(Gasem_5,1,new TGeoCombiTrans(4,50,-5,new TGeoRotation("R2",90,90,0)));
   top->AddNodeOverlap(Gasem_m1,1,new TGeoCombiTrans(-22,30,16,new TGeoRotation("R2",90,88,0)));
   top->AddNodeOverlap(Gasem_m1,1,new TGeoCombiTrans(-22,-30,16,new TGeoRotation("R2",90,88,0)));
   top->AddNodeOverlap(Gasem_m1,1,new TGeoCombiTrans(29,30,16,new TGeoRotation("R2",90,92,0)));
   top->AddNodeOverlap(Gasem_m1,1,new TGeoCombiTrans(29,-30,16,new TGeoRotation("R2",90,92,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,-62,16,new TGeoRotation("R2",0,3,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,62,16,new TGeoRotation("R2",0,-3,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,-30,47.5,new TGeoRotation("R2",0,87,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,30,47.5,new TGeoRotation("R2",0,-87,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,-30,-16,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gasem_m2,1,new TGeoCombiTrans(2,30,-16,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(V,1,new TGeoCombiTrans(-30,18.3,16,new TGeoRotation("R2",0,-135,0)));
   top->AddNodeOverlap(V,1,new TGeoCombiTrans(-30,-18.3,16,new TGeoRotation("R2",0,135,0)));
   top->AddNodeOverlap(V_m,1,new TGeoTranslation(-30,-37,35));
   top->AddNodeOverlap(V_m,1,new TGeoTranslation(-30,37,35));

   // abdomen
   TGeoVolume *Bea=Robot->MakeEltu("Bea",Iron,20,37,25);
     Bea->SetLineColor(17);
     Bea->SetFillColor(17);
   TGeoVolume *Bea_d=Robot->MakeEltu("Bea_d",Iron,21,36,5);
     Bea_d->SetLineColor(12);
     Bea_d->SetFillColor(12);
   TGeoVolume *Beakop=Robot->MakeEltu("Beakop",Iron,15,25,5);
     Beakop->SetLineColor(10);
     Beakop->SetFillColor(10);

   // drawing abdomen
   top->AddNodeOverlap(Bea,1,new TGeoTranslation(3,0,-30));
   top->AddNodeOverlap(Bea_d,1,new TGeoTranslation(3,0,-10));
   top->AddNodeOverlap(Beakop,1,new TGeoCombiTrans(-12.1,0,-50, new TGeoRotation("R2",90,90,0)));

   // Gungdi
   TGeoVolume *Gungdi=Robot->MakeEltu("Gungdi",Iron,25,50,18);
     Gungdi->SetLineColor(12);
     Gungdi->SetFillColor(12);
   TGeoVolume *Gungdi_d=Robot->MakeEltu("Gungdi_d",Iron,5,5,5);
     Gungdi_d->SetLineColor(2);
     Gungdi_d->SetFillColor(2);

   // drawing Gungdi
   top->AddNodeOverlap(Gungdi,1,new TGeoTranslation(3,0,-70));
   for (int i=0; i<30; i++) {
      TGeoVolume *Gungdi_j=Robot->MakeEltu("Gungdi_j",Iron,24-0.2*i,49-0.5*i,1);
        Gungdi_j->SetLineColor(12);
        Gungdi_j->SetFillColor(12);
      top->AddNodeOverlap(Gungdi_j,1,new TGeoTranslation(3,0,-51+0.5*i));
   }
   for (int i=0; i<40; i++) {
      if (i<16) {
        TGeoVolume *Gungdi_h=Robot->MakeEltu("Gungdi_h",Iron,24-0.1*i,49-0.3*i,1);
          Gungdi_h->SetLineColor(12);
          Gungdi_h->SetFillColor(12);
        top->AddNodeOverlap(Gungdi_h,1,new TGeoTranslation(3,0,-88-0.5*i));
      }
      else {
        TGeoVolume *Gungdi_h=Robot->MakeEltu("Gungdi_h",Iron,27-0.3*i,52-0.5*i,1);
          Gungdi_h->SetLineColor(12);
          Gungdi_h->SetFillColor(12);
        top->AddNodeOverlap(Gungdi_h,1,new TGeoTranslation(3,0,-89-0.5*i));
      }
   }
   top->AddNodeOverlap(Gungdi_d,1,new TGeoCombiTrans(3,-45,-62,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gungdi_d,1,new TGeoCombiTrans(3,-45,-78,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gungdi_d,1,new TGeoCombiTrans(3,45,-62,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Gungdi_d,1,new TGeoCombiTrans(3,45,-78,new TGeoRotation("R2",0,90,0)));

   // feet
   TGeoVolume *Jang=Robot->MakeEltu("Jang",Iron,18,18,50);
     Jang->SetLineColor(17);
     Jang->SetFillColor(17);
   TGeoVolume *Jong=Robot->MakeEltu("Jong",Iron,22,22,50);
     Jong->SetLineColor(12);
     Jong->SetFillColor(12);
   TGeoVolume *Bal=Robot->MakeSphere("Bal",Iron,0,22,0,180,180,360);
     Bal->SetLineColor(12);
     Bal->SetFillColor(12);

   // drawing Dary
   top->AddNodeOverlap(Jang,1,new TGeoCombiTrans(3,-25,-120,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Jang,1,new TGeoCombiTrans(3,25,-120,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Jong,1,new TGeoCombiTrans(3,-25,-220,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Jong,1,new TGeoCombiTrans(3,25,-220,new TGeoRotation("R2",0,0,0)));
   for (int i=0; i<30; i++) {
             TGeoVolume *Mu=Robot->MakeCons("Mu",Iron,1,0,22.1,0,22.1,120+2*i,-120-2*i);
        Mu->SetLineColor(4);
        Mu->SetFillColor(4);
      top->AddNodeOverlap(Mu,1,new TGeoTranslation(3,-25,-171-i));
      top->AddNodeOverlap(Mu,1,new TGeoTranslation(3,25,-171-i));

   }
   top->AddNodeOverlap(Bal,1,new TGeoCombiTrans(-10,-25,-270,new TGeoRotation("R2",270,-90,0)));
   top->AddNodeOverlap(Bal,1,new TGeoCombiTrans(-10,25,-270,new TGeoRotation("R2",270,-90,0)));

   // arms
   TGeoVolume *S=Robot->MakeSphere("S",Iron,0,25,0,180,180,360);
     S->SetLineColor(17);
     S->SetFillColor(17);
   TGeoVolume *S_1=Robot->MakeSphere("S_1",Iron,0,15,0,180,0,360);
     S_1->SetLineColor(17);
     S_1->SetFillColor(17);
   TGeoVolume *Pal=Robot->MakeEltu("Pal",Iron,15,15,30);
     Pal->SetLineColor(17);
     Pal->SetFillColor(17);
   TGeoVolume *Fal=Robot->MakeEltu("Fal",Iron,17,17,30);
     Fal->SetLineColor(4);
     Fal->SetFillColor(4);
          TGeoVolume *Bbul=Robot->MakeCone("Bbul",Iron,8,0,0,0,5);
     Bbul->SetLineColor(17);
     Bbul->SetFillColor(17);

   // drawing arms
   top->AddNodeOverlap(S,1,new TGeoCombiTrans(3,73,30,new TGeoRotation("R2",0,-30,0)));
   top->AddNodeOverlap(S,1,new TGeoCombiTrans(3,-73,30,new TGeoRotation("R2",0,210,0)));
   top->AddNodeOverlap(S_1,1,new TGeoCombiTrans(3,-73,27,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(S_1,1,new TGeoCombiTrans(3,73,27,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Pal,1,new TGeoCombiTrans(3,-73,-5,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Pal,1,new TGeoCombiTrans(3,73,-5,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Fal,1,new TGeoCombiTrans(3,-73,-60,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Fal,1,new TGeoCombiTrans(3,73,-60,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Bbul,1,new TGeoCombiTrans(3,-97,-72,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Bbul,1,new TGeoCombiTrans(3,-97,-48,new TGeoRotation("R2",0,-90,0)));
   top->AddNodeOverlap(Bbul,1,new TGeoCombiTrans(3,97,-72,new TGeoRotation("R2",0,90,0)));
   top->AddNodeOverlap(Bbul,1,new TGeoCombiTrans(3,97,-48,new TGeoRotation("R2",0,90,0)));

   // hands
   TGeoVolume *Son_d=Robot->MakeBox("Son_d",Iron,15,15,7);
     Son_d->SetLineColor(4);
     Son_d->SetFillColor(4);
   TGeoVolume *Son_g=Robot->MakeBox("Son_g",Iron,4,10,4);
     Son_g->SetLineColor(4);
     Son_g->SetFillColor(4);
   TGeoVolume *Son_g1=Robot->MakeBox("Son_g1",Iron,6,6,6);
     Son_g1->SetLineColor(4);
     Son_g1->SetFillColor(4);
   TGeoVolume *Son_g2=Robot->MakeBox("Son_g2",Iron,8,3,3);
     Son_g2->SetLineColor(4);
     Son_g2->SetFillColor(4);
          TGeoVolume *Last_b=Robot->MakeCone("Last_b",Iron,10,0,0,0,4);
     Last_b->SetLineColor(17);
     Last_b->SetFillColor(17);
   TGeoVolume *Last=Robot->MakeSphere("Last",Iron,0,3,0,180,0,360);
     Last->SetLineColor(2);
     Last->SetFillColor(2);

   //drawing hands
   top->AddNodeOverlap(Son_d,1,new TGeoCombiTrans(3,-80,-105,new TGeoRotation("R2",0,90,0)));
   for (int i=0; i<4; i++) {
      top->AddNodeOverlap(Son_g,1,new TGeoCombiTrans(-6+6*i,-72,-118,new TGeoRotation("R2",0,-10,0)));
   }
   for (int i=0; i<4; i++) {
      top->AddNodeOverlap(Son_g,1,new TGeoCombiTrans(-6+6*i,-67,-113,new TGeoRotation("R2",0,110,0)));
   }
   top->AddNodeOverlap(Son_g1,1,new TGeoCombiTrans(-5,-70,-98,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Son_g2,1,new TGeoCombiTrans(-5,-65,-102,new TGeoRotation("R2",0,60,0)));
   top->AddNodeOverlap(Son_d,1,new TGeoCombiTrans(3,80,-105,new TGeoRotation("R2",0,90,0)));
   for (int i=0; i<4; i++) {
      top->AddNodeOverlap(Son_g,1,new TGeoCombiTrans(-6+6*i,72,-118,new TGeoRotation("R2",0,10,0)));
   }
   for (int i=0; i<4; i++) {
      top->AddNodeOverlap(Son_g,1,new TGeoCombiTrans(-6+6*i,67,-113,new TGeoRotation("R2",0,70,0)));
   }
   top->AddNodeOverlap(Son_g1,1,new TGeoCombiTrans(-5,70,-98,new TGeoRotation("R2",0,0,0)));
   top->AddNodeOverlap(Son_g2,1,new TGeoCombiTrans(-5,65,-102,new TGeoRotation("R2",0,60,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(3,-88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(12,-88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(-7,-88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(3,88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(12,88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last_b,1,new TGeoCombiTrans(-7,88,-103,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(3,-88,-112,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(12,-88,-112,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(-7,-88,-112,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(3,88,-112,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(12,88,-112,new TGeoRotation("R2",0,180,0)));
   top->AddNodeOverlap(Last,1,new TGeoCombiTrans(-7,88,-112,new TGeoRotation("R2",0,180,0)));

for (int i=1; i<20; i+=1) {
   if (i<7) {
     TGeoVolume *Effect=Robot->MakeCons("Effect",Iron,3,20/sin(i),21/sin(i),20/sin(i),21/sin(i),0,70);
       Effect->SetLineColor(9);
       Effect->SetFillColor(9);
     top->AddNodeOverlap(Effect,1,new TGeoTranslation(3,0,-280));
   }
   if (6<i && i<10) {
     TGeoVolume *Effect=Robot->MakeCons("Effect",Iron,5,20/sin(i),21/sin(i),20/sin(i),21/sin(i),50,120);
       Effect->SetLineColor(38);
       Effect->SetFillColor(38);
     top->AddNodeOverlap(Effect,1,new TGeoTranslation(3,0,-280));
   }
   if (9<i && i<20) {
     TGeoVolume *Effect=Robot->MakeCons("Effect",Iron,4,20/sin(i),21/sin(i),20/sin(i),21/sin(i),200,330);
       Effect->SetLineColor(33);
       Effect->SetFillColor(33);
     top->AddNodeOverlap(Effect,1,new TGeoTranslation(3,0,-280));
   }
}


   //close geometry
   top->SetVisibility(0);
   Robot->CloseGeometry();

   // in GL viewer
   top->Draw("ogl");
}
