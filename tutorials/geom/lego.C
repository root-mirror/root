/// \file
/// \ingroup tutorial_geom
/// Drawing a figure, made of lego block, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x lego.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_lego.png width=800px
/// \macro_code
///
/// \author Soon Gi Kwon(1116won@hanmail.net), Dept. of Physics, Univ. of Seoul

#include "TSystem.h"
#include "TGeoManager.h"

void lego()
{
   TGeoManager *geom = new TGeoManager("geom","My first 3D geometry");


   TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
   TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);



   TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
   TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);


   // create volume
   TGeoVolume *top=geom->MakeBox("top",Air,100,100,100);
   geom->SetTopVolume(top);
   geom->SetTopVisible(0);
   // If you want to see the boundary, please input the number, 1 instead of 0.
   // Like this, geom->SetTopVisible(1);

//----------------------------------------------------------------------

   TGeoVolume *ha1=geom->MakeSphere("ha1",Iron,0,10,80,90,0,360);
   ha1->SetLineColor(41);
   top->AddNodeOverlap(ha1,1,new TGeoCombiTrans(0,0,4,new TGeoRotation("ha1",0,0,0)));

   TGeoVolume *ha2=geom->MakeSphere("ha2",Iron,0,7,90,180,0,360);
   ha2->SetLineColor(41);
   top->AddNodeOverlap(ha2,1,new TGeoCombiTrans(0,0,4,new TGeoRotation("ha2",0,180,0)));

   TGeoVolume *ha3=geom->MakeSphere("ha3",Iron,0,7.3,80,90,0,360);
   ha3->SetLineColor(2);
   top->AddNodeOverlap(ha3,1,new TGeoCombiTrans(0,0,4.8,new TGeoRotation("ha3",0,0,0)));


   TGeoVolume *h1=geom->MakeTubs("h1",Iron,0,6,4.5,0,0);
   h1->SetLineColor(5);
   top->AddNodeOverlap(h1,1,new TGeoCombiTrans(0,0,0,new TGeoRotation("h1",0,0,0)));

   TGeoVolume *h2=geom->MakeSphere("h2",Iron,0,7.5,0,52.5,0,360);
   h2->SetLineColor(5);
   top->AddNodeOverlap(h2,1,new TGeoCombiTrans(0,0,0,new TGeoRotation("h2",0,0,0)));

   TGeoVolume *h3=geom->MakeSphere("h3",Iron,0,7.5,0,52.5,0,360);
   h3->SetLineColor(5);
   top->AddNodeOverlap(h3,1,new TGeoCombiTrans(0,0,0,new TGeoRotation("h3",180,180,0)));

   TGeoVolume *h4=geom->MakeTubs("h4",Iron,2.5,3.5,1.5,0,0);
   h4->SetLineColor(5);
   top->AddNodeOverlap(h4,1,new TGeoCombiTrans(0,0,7.5,new TGeoRotation("h4",0,0,0)));



   TGeoVolume *t1_1=geom->MakeTubs("t1_1",Iron,0,0.8,1,0,360);
   t1_1->SetLineColor(12);
   top->AddNodeOverlap(t1_1,1,new TGeoCombiTrans(-5,2,1.5,new TGeoRotation("t1_1",-90,90,0)));

   TGeoVolume *t2_1=geom->MakeTubs("t2_1",Iron,0,0.8,1,0,360);
   t2_1->SetLineColor(12);
   top->AddNodeOverlap(t2_1,1,new TGeoCombiTrans(-5,-2,1.5,new TGeoRotation("t2_1",-90,90,0)));

   TGeoVolume *fb1=geom->MakeTubs("fb1",Iron,2,2.3,1,100,260);
   fb1->SetLineColor(12);
   top->AddNodeOverlap(fb1,1,new TGeoCombiTrans(-5,0,-1,new TGeoRotation("fb1",90,90,90)));



   TGeoVolume *m1=geom->MakeBox("m1",Iron,7,8,4);
   m1->SetLineColor(2);
   top->AddNodeOverlap(m1,1,new TGeoCombiTrans(0,0,-17,new TGeoRotation("m1",90,90,0)));

   TGeoVolume *m2=geom->MakeTubs("m2",Iron,0,1,7,90,180);
   m2->SetLineColor(2);
   top->AddNodeOverlap(m2,1,new TGeoCombiTrans(-3,0,-9,new TGeoRotation("m2",0,90,0)));

   TGeoVolume *m3=geom->MakeTubs("m3",Iron,0,1,7,0,90);
   m3->SetLineColor(2);
   top->AddNodeOverlap(m3,1,new TGeoCombiTrans(3,0,-9,new TGeoRotation("m3",0,90,0)));

   TGeoVolume *m4=geom->MakeBox("m4",Iron,3,7,0.5);
   m4->SetLineColor(2);
   top->AddNodeOverlap(m4,1,new TGeoCombiTrans(0,0,-8.5,new TGeoRotation("m4",90,0,90)));

   TGeoVolume *m5=geom->MakeTubs("m5",Iron,0,1.5,1.2,0,0);
   m5->SetLineColor(5);
   top->AddNodeOverlap(m5,1,new TGeoCombiTrans(0,0,-7.8,new TGeoRotation("m5",0,0,0)));

   TGeoVolume *m6=geom->MakeTrd2("m6",Iron,4,4,0,2,8);
   m6->SetLineColor(2);
   top->AddNodeOverlap(m6,1,new TGeoCombiTrans(0,-7,-17,new TGeoRotation("m6",0,180,0)));

   TGeoVolume *m7=geom->MakeTrd2("m7",Iron,4,4,0,2,8);
   m7->SetLineColor(2);
   top->AddNodeOverlap(m7,1,new TGeoCombiTrans(0,7,-17,new TGeoRotation("m7",0,180,0)));


   TGeoVolume *md1=geom->MakeBox("md1",Iron,4,8.5,0.7);
   md1->SetLineColor(37);
   top->AddNodeOverlap(md1,1,new TGeoCombiTrans(0,0,-25.5,new TGeoRotation("md1",0,0,0)));

   TGeoVolume *md2=geom->MakeBox("md2",Iron,3,0.4,2);
   md2->SetLineColor(37);
   top->AddNodeOverlap(md2,1,new TGeoCombiTrans(0,0,-28,new TGeoRotation("md2",0,0,0)));

   TGeoVolume *d1=geom->MakeTrd2("d1",Iron,3,4,4,4,7);
   d1->SetLineColor(37);
   top->AddNodeOverlap(d1,1,new TGeoCombiTrans(-4.8,4.5,-35,new TGeoRotation("d1",90,45,-90)));

   TGeoVolume *d2=geom->MakeTrd2("d2",Iron,3,4,4,4,7);
   d2->SetLineColor(37);
   top->AddNodeOverlap(d2,1,new TGeoCombiTrans(0,-4.5,-37,new TGeoRotation("d2",0,0,0)));

   TGeoVolume *d3=geom->MakeTubs("d3",Iron,0,4,3.98,0,180);
   d3->SetLineColor(37);
   top->AddNodeOverlap(d3,1,new TGeoCombiTrans(0,4.5,-30.2,new TGeoRotation("d3",0,90,-45)));

   TGeoVolume *d4=geom->MakeTubs("d4",Iron,0,4,3.98,0,180);
   d4->SetLineColor(37);
   top->AddNodeOverlap(d4,1,new TGeoCombiTrans(0,-4.5,-30,new TGeoRotation("d4",0,90,0)));

   TGeoVolume *d5=geom->MakeBox("d5",Iron,4,4,1);
   d5->SetLineColor(37);
   top->AddNodeOverlap(d5,1,new TGeoCombiTrans(-10.2,4.5,-39,new TGeoRotation("d5",90,45,-90)));

   TGeoVolume *d6=geom->MakeBox("d6",Iron,4,4,1);
   d6->SetLineColor(37);
   top->AddNodeOverlap(d6,1,new TGeoCombiTrans(-1,-4.5,-43.4,new TGeoRotation("d6",0,0,0)));



   TGeoVolume *a1=geom->MakeTubs("a1",Iron,0,1.5,4,0,0);
   a1->SetLineColor(1);
   top->AddNodeOverlap(a1,1,new TGeoCombiTrans(0,10,-15.1,new TGeoRotation("a1",0,20,45)));

   TGeoVolume *a2=geom->MakeSphere("a2",Iron,0,1.48,0,180,0,200);
   a2->SetLineColor(1);
   top->AddNodeOverlap(a2,1,new TGeoCombiTrans(0,8.6,-11.5,new TGeoRotation("a2",120,80,20)));

   TGeoVolume *a3=geom->MakeTubs("a3",Iron,0,1.5,2.2,0,0);
   a3->SetLineColor(1);
   top->AddNodeOverlap(a3,1,new TGeoCombiTrans(0,11.3,-20.6,new TGeoRotation("a3",300,0,40)));

   TGeoVolume *a4=geom->MakeTubs("a4",Iron,0,1,1,0,0);
   a4->SetLineColor(5);
   top->AddNodeOverlap(a4,1,new TGeoCombiTrans(0,11.3,-23.8,new TGeoRotation("a4",75,0,30)));

   TGeoVolume *a5=geom->MakeTubs("a5",Iron,1.5,2.5,2,0,270);
   a5->SetLineColor(5);
   top->AddNodeOverlap(a5,1,new TGeoCombiTrans(0,11.3,-26.5,new TGeoRotation("a5",-90,90,00)));




   TGeoVolume *a1_1=geom->MakeTubs("a1_1",Iron,0,1.5,4,0,0);
   a1_1->SetLineColor(1);
   top->AddNodeOverlap(a1_1,1,new TGeoCombiTrans(0,-10,-15.1,new TGeoRotation("a1_1",0,-20,-45)));

   TGeoVolume *a2_1=geom->MakeSphere("a2_1",Iron,0,1.48,0,180,0,200);
   a2_1->SetLineColor(1);
   top->AddNodeOverlap(a2_1,1,new TGeoCombiTrans(0,-8.6,-11.5,new TGeoRotation("a2_1",120,80,-20)));

   TGeoVolume *a3_1=geom->MakeTubs("a3_1",Iron,0,1.5,2.2,0,0);
   a3_1->SetLineColor(1);
   top->AddNodeOverlap(a3_1,1,new TGeoCombiTrans(0,-11.3,-20.6,new TGeoRotation("a3_1",-300,0,-40)));

   TGeoVolume *a4_1=geom->MakeTubs("a4_1",Iron,0,1,1,0,0);
   a4_1->SetLineColor(5);
   top->AddNodeOverlap(a4_1,1,new TGeoCombiTrans(0,-11.3,-23.8,new TGeoRotation("a4_1",-75,0,-30)));

   a5=geom->MakeTubs("a5_1",Iron,1.5,2.5,2,0,270);
   a5->SetLineColor(5);
   top->AddNodeOverlap(a5,1,new TGeoCombiTrans(0,-11.3,-26.5,new TGeoRotation("a5",90,90,00)));


   //**********************************NO,2******************


   TGeoVolume *ha_1=geom->MakeSphere("ha_1",Iron,0,10,80,90,0,360);
   ha_1->SetLineColor(6);
   top->AddNodeOverlap(ha_1,1,new TGeoCombiTrans(0,36,4,new TGeoRotation("ha_1",0,0,0)));

   TGeoVolume *ha_2=geom->MakeTubs("ha_2",Iron,0,6,5,0,0);
   ha_2->SetLineColor(6);
   top->AddNodeOverlap(ha_2,1,new TGeoCombiTrans(0,36,10,new TGeoRotation("ha_2",0,180,0)));

   TGeoVolume *ha_3=geom->MakeTubs("ha_3",Iron,0,1,12,0,0);
   ha_3->SetLineColor(28);
   top->AddNodeOverlap(ha_3,1,new TGeoCombiTrans(0,36,8,new TGeoRotation("ha_3",0,90,0)));

   TGeoVolume *ha_4=geom->MakeTubs("ha_4",Iron,0,1,3,0,0);
   ha_4->SetLineColor(28);
   top->AddNodeOverlap(ha_4,1,new TGeoCombiTrans(0,22,10,new TGeoRotation("ha_4",0,0,0)));

   TGeoVolume *ha_5=geom->MakeTubs("ha_5",Iron,0,1,3,0,0);
   ha_5->SetLineColor(28);
   top->AddNodeOverlap(ha_5,1,new TGeoCombiTrans(0,46,10,new TGeoRotation("ha_5",0,0,0)));

   TGeoVolume *ha_6=geom->MakeTubs("ha_6",Iron,0,1,3,0,0);
   ha_6->SetLineColor(28);
   top->AddNodeOverlap(ha_6,1,new TGeoCombiTrans(0,24,10,new TGeoRotation("ha_6",0,0,0)));

   TGeoVolume *ha_7=geom->MakeTubs("ha_7",Iron,0,1,3,0,0);
   ha_7->SetLineColor(28);
   top->AddNodeOverlap(ha_7,1,new TGeoCombiTrans(0,48,10,new TGeoRotation("ha_7",0,0,0)));

   TGeoVolume *ha_8=geom->MakeBox("ha_8",Iron,2,0.5,2);
   ha_8->SetLineColor(19);
   top->AddNodeOverlap(ha_8,1,new TGeoCombiTrans(-4.2,36,9,new TGeoRotation("ha_8",0,45,0)));


   TGeoVolume *ha_9=geom->MakeBox("ha_9",Iron,2,0.5,2);
   ha_9->SetLineColor(19);
   top->AddNodeOverlap(ha_9,1,new TGeoCombiTrans(-4.2,36,9,new TGeoRotation("ha_9",0,135,0)));



   TGeoVolume *h_1=geom->MakeTubs("h_1",Iron,0,6,4.5,0,0);
   h_1->SetLineColor(5);
   top->AddNodeOverlap(h_1,1,new TGeoCombiTrans(0,36,0,new TGeoRotation("h_1",0,0,0)));

   TGeoVolume *h_2=geom->MakeSphere("h_2",Iron,0,7.5,0,52.5,0,360);
   h_2->SetLineColor(5);
   top->AddNodeOverlap(h_2,1,new TGeoCombiTrans(0,36,0,new TGeoRotation("h_2",0,0,0)));

   TGeoVolume *h_3=geom->MakeSphere("h_3",Iron,0,7.5,0,52.5,0,360);
   h_3->SetLineColor(5);
   top->AddNodeOverlap(h_3,1,new TGeoCombiTrans(0,36,0,new TGeoRotation("h_3",180,180,0)));

   TGeoVolume *h_4=geom->MakeTubs("h_4",Iron,2.5,3.5,1.5,0,0);
   h_4->SetLineColor(5);
   top->AddNodeOverlap(h_4,1,new TGeoCombiTrans(0,36,7.5,new TGeoRotation("h_4",0,0,0)));


   TGeoVolume *fa1=geom->MakeTubs("fa1",Iron,0,0.5,1,0,360);
   fa1->SetLineColor(12);
   top->AddNodeOverlap(fa1,1,new TGeoCombiTrans(-5,38,1.5,new TGeoRotation("fa1",-90,90,0)));

   TGeoVolume *fa2=geom->MakeTubs("fa2",Iron,0,0.5,1,0,360);
   fa2->SetLineColor(12);
   top->AddNodeOverlap(fa2,1,new TGeoCombiTrans(-5,34,1.5,new TGeoRotation("fa2",-90,90,0)));

   TGeoVolume *fa1_1=geom->MakeTubs("fa1_1",Iron,1,1.2,1,0,360);
   fa1_1->SetLineColor(12);
   top->AddNodeOverlap(fa1_1,1,new TGeoCombiTrans(-5,38,1.5,new TGeoRotation("fa1_1",-90,90,0)));

   TGeoVolume *fa2_1=geom->MakeTubs("fa2_1",Iron,1,1.2,1,0,360);
   fa2_1->SetLineColor(12);
   top->AddNodeOverlap(fa2_1,1,new TGeoCombiTrans(-5,34,1.5,new TGeoRotation("fa2_1",-90,90,0)));

   TGeoVolume *fa3=geom->MakeTubs("fa3",Iron,2,2.3,1,90,270);
   fa3->SetLineColor(12);
   top->AddNodeOverlap(fa3,1,new TGeoCombiTrans(-5,36,-1,new TGeoRotation("fa3",90,90,90)));



   TGeoVolume *m_1=geom->MakeBox("m_1",Iron,7,8,4);
   m_1->SetLineColor(25);
   top->AddNodeOverlap(m_1,1,new TGeoCombiTrans(0,36,-17,new TGeoRotation("m_1",90,90,0)));

   TGeoVolume *m_2=geom->MakeTubs("m_2",Iron,0,1,7,90,180);
   m_2->SetLineColor(25);
   top->AddNodeOverlap(m_2,1,new TGeoCombiTrans(-3,36,-9,new TGeoRotation("m_2",0,90,0)));

   TGeoVolume *m_3=geom->MakeTubs("m_3",Iron,0,1,7,0,90);
   m_3->SetLineColor(25);
   top->AddNodeOverlap(m_3,1,new TGeoCombiTrans(3,36,-9,new TGeoRotation("m_3",0,90,0)));

   TGeoVolume *m_4=geom->MakeBox("m_4",Iron,3,7,0.5);
   m_4->SetLineColor(25);
   top->AddNodeOverlap(m_4,1,new TGeoCombiTrans(0,36,-8.5,new TGeoRotation("m_4",90,0,90)));

   TGeoVolume *m_5=geom->MakeTubs("m_5",Iron,0,1.5,1.2,0,0);
   m_5->SetLineColor(5);
   top->AddNodeOverlap(m_5,1,new TGeoCombiTrans(0,36,-7.8,new TGeoRotation("m_5",0,0,0)));

   TGeoVolume *m_6=geom->MakeTrd2("m_6",Iron,4,4,0,2,8);
   m_6->SetLineColor(25);
   top->AddNodeOverlap(m_6,1,new TGeoCombiTrans(0,29,-17,new TGeoRotation("m_6",0,180,0)));

   TGeoVolume *m_7=geom->MakeTrd2("m_7",Iron,4,4,0,2,8);
   m_7->SetLineColor(25);
   top->AddNodeOverlap(m_7,1,new TGeoCombiTrans(0,43,-17,new TGeoRotation("m_7",0,180,0)));


   TGeoVolume *md_1=geom->MakeBox("md_1",Iron,4,8.5,0.7);
   md_1->SetLineColor(48);
   top->AddNodeOverlap(md_1,1,new TGeoCombiTrans(0,36,-25.5,new TGeoRotation("md_1",0,0,0)));

   TGeoVolume *md_2=geom->MakeBox("md_2",Iron,3,0.4,2);
   md_2->SetLineColor(48);
   top->AddNodeOverlap(md_2,1,new TGeoCombiTrans(0,36,-28,new TGeoRotation("md_2",0,0,0)));

   TGeoVolume *d_1=geom->MakeTrd2("d_1",Iron,3,4,4,4,7);
   d_1->SetLineColor(48);
   top->AddNodeOverlap(d_1,1,new TGeoCombiTrans(0,40.5,-37.2,new TGeoRotation("d_1",0,0,0)));

   TGeoVolume *d_2=geom->MakeTrd2("d_2",Iron,3,4,4,4,7);
   d_2->SetLineColor(48);
   top->AddNodeOverlap(d_2,1,new TGeoCombiTrans(0,31.5,-37.2,new TGeoRotation("d_2",0,0,0)));

   TGeoVolume *d_3=geom->MakeTubs("d_3",Iron,0,4,3.98,0,180);
   d_3->SetLineColor(48);
   top->AddNodeOverlap(d_3,1,new TGeoCombiTrans(0,40.5,-30.2,new TGeoRotation("d_3",0,90,0)));

   TGeoVolume *d_4=geom->MakeTubs("d_4",Iron,0,4,3.98,0,180);
   d_4->SetLineColor(48);
   top->AddNodeOverlap(d_4,1,new TGeoCombiTrans(0,31.5,-30.2,new TGeoRotation("d_4",0,90,0)));

   TGeoVolume *d_5=geom->MakeBox("d_5",Iron,4,4,1);
   d_5->SetLineColor(48);
   top->AddNodeOverlap(d_5,1,new TGeoCombiTrans(-1,40.5,-43.7,new TGeoRotation("d_5",0,0,0)));

   TGeoVolume *d_6=geom->MakeBox("d_6",Iron,4,4,1);
   d_6->SetLineColor(48);
   top->AddNodeOverlap(d_6,1,new TGeoCombiTrans(-1,31.5,-43.7,new TGeoRotation("d_6",0,0,0)));


   TGeoVolume *a_1=geom->MakeTubs("a_1",Iron,0,1.5,4,0,0);
   a_1->SetLineColor(45);
   top->AddNodeOverlap(a_1,1,new TGeoCombiTrans(0,46,-15.1,new TGeoRotation("a_1",0,20,45)));

   TGeoVolume *a_2=geom->MakeSphere("a_2",Iron,0,1.48,0,180,0,200);
   a_2->SetLineColor(45);
   top->AddNodeOverlap(a_2,1,new TGeoCombiTrans(0,44.6,-11.5,new TGeoRotation("a_2",120,80,20)));

   TGeoVolume *a_3=geom->MakeTubs("a_3",Iron,0,1.5,2.2,0,0);
   a_3->SetLineColor(45);
   top->AddNodeOverlap(a_3,1,new TGeoCombiTrans(0,47.3,-20.6,new TGeoRotation("a_3",300,0,40)));

   TGeoVolume *a_4=geom->MakeTubs("a_4",Iron,0,1,1,0,0);
   a_4->SetLineColor(12);
   top->AddNodeOverlap(a_4,1,new TGeoCombiTrans(0,47.3,-23.8,new TGeoRotation("a_4",75,0,30)));

   TGeoVolume *a_5=geom->MakeTubs("a_5",Iron,1.5,2.5,2,0,270);
   a_5->SetLineColor(12);
   top->AddNodeOverlap(a_5,1,new TGeoCombiTrans(0,47.3,-26.5,new TGeoRotation("a_5",-90,90,0)));




   TGeoVolume *Aa1=geom->MakeTubs("Aa1",Iron,0,1.5,4,0,0);
   Aa1->SetLineColor(45);
   top->AddNodeOverlap(Aa1,1,new TGeoCombiTrans(0,26,-15.1,new TGeoRotation("Aa1",0,-20,-45)));

   TGeoVolume *Aa2=geom->MakeSphere("Aa2",Iron,0,1.48,0,180,0,200);
   Aa2->SetLineColor(45);
   top->AddNodeOverlap(Aa2,1,new TGeoCombiTrans(0,27.4,-11.5,new TGeoRotation("Aa2",120,80,-20)));

   TGeoVolume *Aa3=geom->MakeTubs("Aa3",Iron,0,1.5,2.2,0,0);
   Aa3->SetLineColor(45);
   top->AddNodeOverlap(Aa3,1,new TGeoCombiTrans(0,24.7,-20.6,new TGeoRotation("Aa3",-300,0,-40)));

   TGeoVolume *Aa4=geom->MakeTubs("Aa4",Iron,0,1,1,0,0);
   Aa4->SetLineColor(12);
   top->AddNodeOverlap(Aa4,1,new TGeoCombiTrans(0,24.7,-23.8,new TGeoRotation("Aa4",-75,0,-30)));

   TGeoVolume *Aa5=geom->MakeTubs("Aa5",Iron,1.5,2.5,2,0,270);
   Aa5->SetLineColor(12);
   top->AddNodeOverlap(Aa5,1,new TGeoCombiTrans(0,24.7,-26.5,new TGeoRotation("Aa5",90,90,00)));



   TGeoVolume *bag1=geom->MakeBox("bag1",Iron,10,4,6);
   bag1->SetLineColor(19);
   top->AddNodeOverlap(bag1,1,new TGeoCombiTrans(0,48,-36,new TGeoRotation("bag1",0,0,0)));

   TGeoVolume *bag2=geom->MakeTubs("bag2",Iron,3,4,1,180,360);
   bag2->SetLineColor(19);
   top->AddNodeOverlap(bag2,1,new TGeoCombiTrans(0,48,-30,new TGeoRotation("bag2",0,270,0)));


   TGeoVolume *well=geom->MakeBox("well",Iron,5,10,3);
   well->SetLineColor(18);
   top->AddNodeOverlap(well,1,new TGeoCombiTrans(-26.5,-17,-42,new TGeoRotation("well",0,0,0)));


   TGeoVolume *K5=geom->MakeTubs("K5",Iron,0,3,3,0,0);
   K5->SetLineColor(18);
   top->AddNodeOverlap(K5,1,new TGeoCombiTrans(-27,-12.5,-39,new TGeoRotation("K5",0,0,0)));

   TGeoVolume *K4=geom->MakeTubs("K4",Iron,0,3,3,0,0);
   K4->SetLineColor(18);
   top->AddNodeOverlap(K4,1,new TGeoCombiTrans(-27,-21.5,-39,new TGeoRotation("K4",0,0,0)));



   //==============Board=========
   char nB[100];
   int Z=0,Y=0;
   TGeoVolume *bo1;

   while(Y<6){
      while(Z<10){
         sprintf(nB,"B%d_Y%d",Z,Y);
         bo1=geom->MakeTubs(nB,Iron,0,3,3,0,0);
         bo1->SetLineColor(8);
         top->AddNodeOverlap(bo1,1,new TGeoCombiTrans(-27+(Y*9),-21.5+(Z*9),-45,new TGeoRotation("bo1",0,0,0)));
         Z++;
      }
      Y++; Z=0;
   }


   TGeoVolume *bo2=geom->MakeBox("bo2",Iron,27,45,3);
   bo2->SetLineColor(8);
   top->AddNodeOverlap(bo2,1,new TGeoCombiTrans(-4.5,18,-48,new TGeoRotation("bo2",0,0,0)));


   top->SetVisibility(0);
   geom->CloseGeometry();

   top->Draw("ogl");
}
