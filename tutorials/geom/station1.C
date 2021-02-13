/// \file
/// \ingroup tutorial_geom
/// Drawing a space station, using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x station1.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_station1.png width=800px
/// \macro_code
///
/// \author Chang Yeol Lee, Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void station1()
{


 TGeoManager *geom=new TGeoManager("geom","My first 3D geometry");

 TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
 TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);
 TGeoMaterial *Cu=new TGeoMaterial("Cu",63.549,29,8.92);

 TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
 TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);
 TGeoMedium *Copper=new TGeoMedium("Copper",2,Cu);

 TGeoVolume *top=geom->MakeBox("top",Air,1000,1000,1000);
 geom->SetTopVolume(top);
 geom->SetTopVisible(0);
 // If you want to see the boundary, please input the number, 1 instead of 0.
 // Like this, geom->SetTopVisible(1);



 TGeoVolume *Cone1=geom->MakeCone("Cone1",Copper,650,0,20,0,20);
 Cone1->SetFillColor(35);
 Cone1->SetLineColor(35);
 top->AddNodeOverlap(Cone1,1,new TGeoTranslation(0,0,0));

 TGeoVolume *Cone2=geom->MakeCone("Cone2",Copper,25,0,30,0,30);
 Cone2->SetFillColor(7);
 Cone2->SetLineColor(7);
 top->AddNodeOverlap(Cone2,1,new TGeoTranslation(0,0,630));

 TGeoVolume *Cone21=geom->MakeCone("Cone21",Copper,30,0,30,0,30);
 Cone21->SetFillColor(29);
 Cone21->SetLineColor(29);
 top->AddNodeOverlap(Cone21,1,new TGeoTranslation(0,0,550));

 TGeoVolume *Cone22=geom->MakeCone("Cone22",Copper,5,0,50,0,50);
 Cone22->SetFillColor(2);
 Cone22->SetLineColor(2);
 top->AddNodeOverlap(Cone22,1,new TGeoTranslation(0,0,500));

for(int i=0;i<28;i++){
 TGeoVolume *Cone00=geom->MakeCone("Cone00",Copper,3,0,25,0,25);
 Cone00->SetFillColor(1);
 Cone00->SetLineColor(1);
 top->AddNodeOverlap(Cone00,1,new TGeoTranslation(0,0,-100+20*i));
}

 TGeoVolume *Cone3=geom->MakeCone("Cone3",Copper,60,0,70,0,0);
 Cone3->SetFillColor(13);
 Cone3->SetLineColor(13);
 top->AddNodeOverlap(Cone3,1,new TGeoTranslation(-60,0,-110));

 TGeoVolume *Cone31=geom->MakeCone("Cone31",Copper,230,0,70,0,70);
 Cone31->SetFillColor(13);
 Cone31->SetLineColor(13);
 top->AddNodeOverlap(Cone31,1,new TGeoTranslation(-60,0,-400));

for(int i=0;i<5;i++){
 Cone31=geom->MakeCone("Cone31",Copper,7,0,73,0,73);
 Cone31->SetFillColor(21);
 Cone31->SetLineColor(21);
 top->AddNodeOverlap(Cone31,1,new TGeoTranslation(-60,0,-170-(500/6*(i+1))));
 top->AddNodeOverlap(Cone31,1,new TGeoTranslation(60,0,-170-(500/6*(i+1))));
 top->AddNodeOverlap(Cone31,1,new TGeoTranslation(0,-60,-170-(500/6*(i+1))));
 top->AddNodeOverlap(Cone31,1,new TGeoTranslation(0,60,-170-(500/6*(i+1))));
}

 TGeoVolume *Cone32=geom->MakeCone("Cone32",Copper,30,60,50,0,70);
 Cone32->SetFillColor(35);
 Cone32->SetLineColor(35);
 top->AddNodeOverlap(Cone32,1,new TGeoTranslation(-60,0,-650));

 TGeoVolume *Cone321=geom->MakeCone("Cone321",Copper,5,60,50,0,50);
 Cone321->SetFillColor(2);
 Cone321->SetLineColor(2);
 top->AddNodeOverlap(Cone321,1,new TGeoTranslation(-60,0,-680));

 TGeoVolume *Cone4=geom->MakeCone("Cone4",Copper,60,0,70,0,0);
 Cone4->SetFillColor(13);
 Cone4->SetLineColor(13);
 top->AddNodeOverlap(Cone3,1,new TGeoTranslation(60,0,-110));

for(int i=1;i<=8;i++){
 TGeoVolume *Torus2=geom->MakeTorus("Torus2",Iron,120,20,40,45*i-4,8);
 Torus2->SetFillColor(18);
 Torus2->SetLineColor(18);
 top->AddNodeOverlap(Torus2,1,new TGeoTranslation(0,0,610));

 TGeoVolume *Tubs=geom->MakeTubs("Line",Iron,0,190,5,45*i-1,45*i+1);
 Tubs->SetFillColor(18);
 Tubs->SetLineColor(18);
 top->AddNodeOverlap(Tubs,1,new TGeoTranslation(0,0,610));
}

 TGeoVolume *Cone41=geom->MakeCone("Cone41",Copper,230,0,70,0,70);
 Cone41->SetFillColor(13);
 Cone41->SetLineColor(13);
 top->AddNodeOverlap(Cone41,1,new TGeoTranslation(60,0,-400));

 TGeoVolume *Cone42=geom->MakeCone("Cone42",Copper,30,60,50,0,70);
 Cone42->SetFillColor(35);
 Cone42->SetLineColor(35);
 top->AddNodeOverlap(Cone42,1,new TGeoTranslation(60,0,-650));

 TGeoVolume *Cone421=geom->MakeCone("Cone421",Copper,5,60,50,0,50);
 Cone421->SetFillColor(2);
 Cone421->SetLineColor(2);
 top->AddNodeOverlap(Cone421,1,new TGeoTranslation(60,0,-680));

 TGeoVolume *Cone5=geom->MakeCone("Cone5",Copper,60,0,70,0,0);
 Cone5->SetFillColor(13);
 Cone5->SetLineColor(13);
 top->AddNodeOverlap(Cone3,1,new TGeoTranslation(0,-60,-110));

 TGeoVolume *Cone51=geom->MakeCone("Cone51",Copper,230,0,70,0,70);
 Cone51->SetFillColor(13);
 Cone51->SetLineColor(13);
 top->AddNodeOverlap(Cone51,1,new TGeoTranslation(0,-60,-400));

 TGeoVolume *Cone52=geom->MakeCone("Cone52",Copper,30,60,50,0,70);
 Cone52->SetFillColor(35);
 Cone52->SetLineColor(35);
 top->AddNodeOverlap(Cone52,1,new TGeoTranslation(0,-60,-650));

 TGeoVolume *Cone521=geom->MakeCone("Cone521",Copper,5,60,50,0,50);
 Cone521->SetFillColor(2);
 Cone521->SetLineColor(2);
 top->AddNodeOverlap(Cone521,1,new TGeoTranslation(0,-60,-680));

 TGeoVolume *Cone6=geom->MakeCone("Cone6",Copper,60,0,70,0,0);
 Cone6->SetFillColor(13);
 Cone6->SetLineColor(13);
 top->AddNodeOverlap(Cone3,1,new TGeoTranslation(0,60,-110));

 TGeoVolume *Cone61=geom->MakeCone("Cone61",Copper,230,0,70,0,70);
 Cone61->SetFillColor(13);
 Cone61->SetLineColor(13);
 top->AddNodeOverlap(Cone61,1,new TGeoTranslation(0,60,-400));

 TGeoVolume *Cone62=geom->MakeCone("Cone62",Copper,30,60,50,0,70);
 Cone62->SetFillColor(35);
 Cone62->SetLineColor(35);
 top->AddNodeOverlap(Cone62,1,new TGeoTranslation(0,60,-650));

 TGeoVolume *Cone621=geom->MakeCone("Cone621",Copper,5,60,50,0,50);
 Cone621->SetFillColor(2);
 Cone621->SetLineColor(2);
 top->AddNodeOverlap(Cone621,1,new TGeoTranslation(0,60,-680));

 TGeoVolume *Cone7=geom->MakeCone("Cone7",Copper,50,0,40,0,5);
 Cone7->SetFillColor(13);
 Cone7->SetLineColor(13);
 top->AddNodeOverlap(Cone7,1,new TGeoCombiTrans(-90,-60,10,new TGeoRotation("Cone7",90,-90,-90)));

 TGeoVolume *Cone71=geom->MakeCone("Cone71",Copper,50,0,60,0,40);
 Cone71->SetFillColor(16);
 Cone71->SetLineColor(16);
 top->AddNodeOverlap(Cone71,1,new TGeoCombiTrans(10,-60,10,new TGeoRotation("Cone7",90,-90,-90)));

 TGeoVolume *Cone711=geom->MakeCone("Cone711",Copper,10,0,10,0,60);
 Cone711->SetFillColor(13);
 Cone711->SetLineColor(13);
 top->AddNodeOverlap(Cone711,1,new TGeoCombiTrans(70,-60,10,new TGeoRotation("Cone7",90,-90,-90)));

 TGeoVolume *Torus1=geom->MakeTorus("Torus1",Iron,120,30,20);
 Torus1->SetFillColor(33);
 Torus1->SetLineColor(33);
 top->AddNodeOverlap(Torus1,1,new TGeoTranslation(0,0,610));

 TGeoVolume *Cone8=geom->MakeCone("Cone8",Copper,50,0,40,0,5);
 Cone8->SetFillColor(13);
 Cone8->SetLineColor(13);
 top->AddNodeOverlap(Cone8,1,new TGeoCombiTrans(100,60,10,new TGeoRotation("Cone8",90,90,0)));

 TGeoVolume *Cone81=geom->MakeCone("Cone81",Copper,50,0,60,0,40);
 Cone81->SetFillColor(16);
 Cone81->SetLineColor(16);
 top->AddNodeOverlap(Cone81,1,new TGeoCombiTrans(0,60,10,new TGeoRotation("Cone8",90,90,0)));

 TGeoVolume *Cone811=geom->MakeCone("Cone811",Copper,10,0,10,0,60);
 Cone811->SetFillColor(13);
 Cone811->SetLineColor(13);
 top->AddNodeOverlap(Cone811,1,new TGeoCombiTrans(-60,60,10,new TGeoRotation("Cone8",90,90,0)));

 TGeoVolume *Box1=geom->MakeBox("Box1",Copper,10,10,3);
 Box1->SetFillColor(3);
 Box1->SetLineColor(3);
 top->AddNodeOverlap(Box1,1,new TGeoCombiTrans(-110,-50,645,new TGeoRotation("Box1",0,0,30)));

 TGeoVolume *Box2=geom->MakeBox("Box2",Copper,10,10,3);
 Box2->SetFillColor(3);
 Box2->SetLineColor(3);
 top->AddNodeOverlap(Box2,1,new TGeoCombiTrans(110,45,645,new TGeoRotation("Box2",0,0,30)));

 TGeoVolume *Box3=geom->MakeBox("Box3",Copper,10,10,3);
 Box3->SetFillColor(3);
 Box3->SetLineColor(3);
 top->AddNodeOverlap(Box3,1,new TGeoCombiTrans(-45,-110,645,new TGeoRotation("Box3",0,0,70)));

 TGeoVolume *Box4=geom->MakeBox("Box4",Copper,10,10,3);
 Box4->SetFillColor(3);
 Box4->SetLineColor(3);
 top->AddNodeOverlap(Box4,1,new TGeoCombiTrans(45,110,645,new TGeoRotation("Box4",0,0,70)));

 TGeoVolume *Box5=geom->MakeBox("Box5",Copper,10,10,3);
 Box5->SetFillColor(3);
 Box5->SetLineColor(3);
 top->AddNodeOverlap(Box5,1,new TGeoCombiTrans(45,-110,645,new TGeoRotation("Box5",0,0,30)));

 TGeoVolume *Box6=geom->MakeBox("Box6",Copper,10,10,3);
 Box6->SetFillColor(3);
 Box6->SetLineColor(3);
 top->AddNodeOverlap(Box6,1,new TGeoCombiTrans(-45,110,645,new TGeoRotation("Box6",0,0,25)));

 TGeoVolume *Box7=geom->MakeBox("Box7",Copper,10,10,3);
 Box7->SetFillColor(3);
 Box7->SetLineColor(3);
 top->AddNodeOverlap(Box7,1,new TGeoCombiTrans(110,-50,645,new TGeoRotation("Box7",0,0,60)));

 TGeoVolume *Box8=geom->MakeBox("Box8",Copper,10,10,3);
 Box8->SetFillColor(3);
 Box8->SetLineColor(3);
 top->AddNodeOverlap(Box8,1,new TGeoCombiTrans(-110,45,645,new TGeoRotation("Box8",0,0,60)));

 Torus1=geom->MakeTorus("Torus1",Iron,120,30,20);
 Torus1->SetFillColor(33);
 Torus1->SetLineColor(33);
 top->AddNodeOverlap(Torus1,1,new TGeoTranslation(0,0,610));

for(int i=1;i<=8;i++){
 TGeoVolume *Torus2=geom->MakeTorus("Torus2",Iron,120,20,40,45*i-4,8);
 Torus2->SetFillColor(18);
 Torus2->SetLineColor(18);
 top->AddNodeOverlap(Torus2,1,new TGeoTranslation(0,0,610));

 TGeoVolume *Tubs=geom->MakeTubs("Line",Iron,0,190,5,45*i-1,45*i+1);
 Tubs->SetFillColor(18);
 Tubs->SetLineColor(18);
 top->AddNodeOverlap(Tubs,1,new TGeoTranslation(0,0,610));
}

 TGeoVolume *Sphere00=geom->MakeSphere("Sphere00",Iron,0,15,0,45,0);
 Sphere00->SetFillColor(2);
 Sphere00->SetLineColor(2);
 top->AddNodeOverlap(Sphere00,1,new TGeoTranslation(-145,-145,600));

 TGeoVolume *Sphere01=geom->MakeSphere("Sphere01",Iron,0,15,0,45,0);
 Sphere01->SetFillColor(2);
 Sphere01->SetLineColor(2);
 top->AddNodeOverlap(Sphere01,1,new TGeoTranslation(0,-210,600));

 TGeoVolume *Sphere02=geom->MakeSphere("Sphere02",Iron,0,15,0,45,0);
 Sphere02->SetFillColor(2);
 Sphere02->SetLineColor(2);
 top->AddNodeOverlap(Sphere02,1,new TGeoTranslation(145,145,600));

 TGeoVolume *Sphere03=geom->MakeSphere("Sphere03",Iron,0,15,0,45,0);
 Sphere03->SetFillColor(2);
 Sphere03->SetLineColor(2);
 top->AddNodeOverlap(Sphere03,1,new TGeoTranslation(0,210,600));

 TGeoVolume *Sphere04=geom->MakeSphere("Sphere04",Iron,0,15,0,45,0);
 Sphere04->SetFillColor(2);
 Sphere04->SetLineColor(2);
 top->AddNodeOverlap(Sphere04,1,new TGeoTranslation(145,-145,600));

 TGeoVolume *Sphere05=geom->MakeSphere("Sphere05",Iron,0,15,0,45,0);
 Sphere05->SetFillColor(2);
 Sphere05->SetLineColor(2);
 top->AddNodeOverlap(Sphere05,1,new TGeoTranslation(-210,0,600));

 TGeoVolume *Sphere06=geom->MakeSphere("Sphere06",Iron,0,15,0,45,0);
 Sphere06->SetFillColor(2);
 Sphere06->SetLineColor(2);
 top->AddNodeOverlap(Sphere06,1,new TGeoTranslation(210,0,600));

 TGeoVolume *Sphere07=geom->MakeSphere("Sphere07",Iron,0,15,0,45,0);
 Sphere07->SetFillColor(2);
 Sphere07->SetLineColor(2);
 top->AddNodeOverlap(Sphere07,1,new TGeoTranslation(-145,145,600));

 TGeoVolume *Torus3=geom->MakeTorus("Torus3",Iron,190,0,10);
 Torus3->SetFillColor(18);
 Torus3->SetLineColor(18);
 top->AddNodeOverlap(Torus3,1,new TGeoTranslation(0,0,610));

 TGeoVolume *Sphere1=geom->MakeSphere("Sphere1",Iron,0,20,0,180,0,360);
 Sphere1->SetFillColor(2);
 Sphere1->SetLineColor(2);
 top->AddNodeOverlap(Sphere1,1,new TGeoTranslation(0,0,650));

 TGeoVolume *Tubs=geom->MakeTubs("Tubs",Iron,0,40,50,0,360);
 Tubs->SetFillColor(29);
 Tubs->SetLineColor(29);
 top->AddNodeOverlap(Tubs,1,new TGeoTranslation(0,0,500));

 TGeoVolume *Tubs1=geom->MakeTubs("Tubs1",Iron,50,60,230,40,150);
 Tubs1->SetFillColor(18);
 Tubs1->SetLineColor(18);
 top->AddNodeOverlap(Tubs1,1,new TGeoTranslation(-170,-30,-400));

 TGeoVolume *Tubs11=geom->MakeTubs("Tubs11",Iron,50,60,230,220,330);
 Tubs11->SetFillColor(18);
 Tubs11->SetLineColor(18);
 top->AddNodeOverlap(Tubs11,1,new TGeoTranslation(-260,35,-400));

 TGeoVolume *Sphere111=geom->MakeSphere("Sphere111",Iron,0,10,0,180,0,360);
 Sphere111->SetFillColor(2);
 Sphere111->SetLineColor(2);
 top->AddNodeOverlap(Sphere111,1,new TGeoTranslation(-310,0,-165));

 TGeoVolume *Sphere112=geom->MakeSphere("Sphere112",Iron,0,10,0,180,0,360);
 Sphere112->SetFillColor(2);
 Sphere112->SetLineColor(2);
 top->AddNodeOverlap(Sphere112,1,new TGeoTranslation(-310,0,-400));

 TGeoVolume *Sphere113=geom->MakeSphere("Sphere113",Iron,0,10,0,180,0,360);
 Sphere113->SetFillColor(2);
 Sphere113->SetLineColor(2);
 top->AddNodeOverlap(Sphere113,1,new TGeoTranslation(-310,0,-635));

 TGeoVolume *Tubs2=geom->MakeTubs("Tubs2",Iron,50,60,230,220,330);
 Tubs2->SetFillColor(18);
 Tubs2->SetLineColor(18);
 top->AddNodeOverlap(Tubs2,1,new TGeoTranslation(170,30,-400));

 TGeoVolume *Tubs21=geom->MakeTubs("Tubs21",Iron,50,60,230,400,510);
 Tubs21->SetFillColor(18);
 Tubs21->SetLineColor(18);
 top->AddNodeOverlap(Tubs21,1,new TGeoTranslation(265,-25,-400));

 TGeoVolume *Sphere211=geom->MakeSphere("Sphere211",Iron,0,10,0,180,0,360);
 Sphere211->SetFillColor(2);
 Sphere211->SetLineColor(2);
 top->AddNodeOverlap(Sphere211,1,new TGeoTranslation(310,0,-165));

 TGeoVolume *Sphere212=geom->MakeSphere("Sphere212",Iron,0,10,0,180,0,360);
 Sphere212->SetFillColor(2);
 Sphere212->SetLineColor(2);
 top->AddNodeOverlap(Sphere212,1,new TGeoTranslation(310,0,-400));

 TGeoVolume *Sphere213=geom->MakeSphere("Sphere213",Iron,0,10,0,180,0,360);
 Sphere213->SetFillColor(2);
 Sphere213->SetLineColor(2);
 top->AddNodeOverlap(Sphere213,1,new TGeoTranslation(310,0,-635));

 TGeoVolume *Tubs3=geom->MakeTubs("Tubs3",Iron,50,60,230,130,260);
 Tubs3->SetFillColor(18);
 Tubs3->SetLineColor(18);
 top->AddNodeOverlap(Tubs3,1,new TGeoTranslation(30,-170,-400));

 TGeoVolume *Tubs31=geom->MakeTubs("Tubs31",Iron,50,60,230,310,440);
 Tubs31->SetFillColor(18);
 Tubs31->SetLineColor(18);
 top->AddNodeOverlap(Tubs31,1,new TGeoTranslation(0,-275,-400));

 TGeoVolume *Sphere311=geom->MakeSphere("Sphere311",Iron,0,10,0,180,0,360);
 Sphere311->SetFillColor(2);
 Sphere311->SetLineColor(2);
 top->AddNodeOverlap(Sphere311,1,new TGeoTranslation(-35,320,-165));

 TGeoVolume *Sphere312=geom->MakeSphere("Sphere312",Iron,0,10,0,180,0,360);
 Sphere312->SetFillColor(2);
 Sphere312->SetLineColor(2);
 top->AddNodeOverlap(Sphere312,1,new TGeoTranslation(-35,320,-400));

 TGeoVolume *Sphere313=geom->MakeSphere("Sphere313",Iron,0,10,0,180,0,360);
 Sphere313->SetFillColor(2);
 Sphere313->SetLineColor(2);
 top->AddNodeOverlap(Sphere313,1,new TGeoTranslation(-35,320,-635));

 TGeoVolume *Tubs4=geom->MakeTubs("Tubs4",Iron,50,60,230,310,440);
 Tubs4->SetFillColor(18);
 Tubs4->SetLineColor(18);
 top->AddNodeOverlap(Tubs4,1,new TGeoTranslation(-30,170,-400));

 TGeoVolume *Tubs41=geom->MakeTubs("Tubs41",Iron,50,60,230,490,620);
 Tubs41->SetFillColor(18);
 Tubs41->SetLineColor(18);
 top->AddNodeOverlap(Tubs41,1,new TGeoTranslation(0,275,-400));

 TGeoVolume *Sphere411=geom->MakeSphere("Sphere411",Iron,0,10,0,180,0,360);
 Sphere411->SetFillColor(2);
 Sphere411->SetLineColor(2);
 top->AddNodeOverlap(Sphere411,1,new TGeoTranslation(30,-320,-165));

 TGeoVolume *Sphere412=geom->MakeSphere("Sphere412",Iron,0,10,0,180,0,360);
 Sphere412->SetFillColor(2);
 Sphere412->SetLineColor(2);
 top->AddNodeOverlap(Sphere412,1,new TGeoTranslation(30,-320,-400));

 TGeoVolume *Sphere413=geom->MakeSphere("Sphere413",Iron,0,10,0,180,0,360);
 Sphere413->SetFillColor(2);
 Sphere413->SetLineColor(2);
 top->AddNodeOverlap(Sphere413,1,new TGeoTranslation(30,-320,-635));

 TGeoVolume *Cone010=geom->MakeCone("Cone010",Iron,30,0,30,0,30);
 Cone010->SetFillColor(2);
 Cone010->SetLineColor(2);
 top->AddNodeOverlap(Cone010,1,new TGeoTranslation(0,0,250));

 TGeoVolume *Torus010=geom->MakeTorus("Torus010",Iron,300,50,40);
 Torus010->SetFillColor(33);
 Torus010->SetLineColor(33);
 top->AddNodeOverlap(Torus010,1,new TGeoTranslation(0,0,250));

 TGeoVolume *Torus011=geom->MakeTorus("Torus011",Iron,400,10,10);
 Torus011->SetFillColor(33);
 Torus011->SetLineColor(33);
 top->AddNodeOverlap(Torus011,1,new TGeoTranslation(0,0,250));

 TGeoVolume *Torus012=geom->MakeTorus("Torus012",Iron,200,10,10);
 Torus012->SetFillColor(33);
 Torus012->SetLineColor(33);
 top->AddNodeOverlap(Torus012,1,new TGeoTranslation(0,0,250));

 TGeoVolume *Sphere010=geom->MakeSphere("Sphere010",Iron,0,10,0,180,0,360);
 Sphere010->SetFillColor(2);
 Sphere010->SetLineColor(2);
 top->AddNodeOverlap(Sphere010,1,new TGeoTranslation(-290,-290,250));

 TGeoVolume *Sphere011=geom->MakeSphere("Sphere011",Iron,0,10,0,180,0,360);
 Sphere011->SetFillColor(2);
 Sphere011->SetLineColor(2);
 top->AddNodeOverlap(Sphere011,1,new TGeoTranslation(290,290,250));

 TGeoVolume *Sphere012=geom->MakeSphere("Sphere012",Iron,0,10,0,180,0,360);
 Sphere012->SetFillColor(2);
 Sphere012->SetLineColor(2);
 top->AddNodeOverlap(Sphere012,1,new TGeoTranslation(0,-410,250));

 TGeoVolume *Sphere013=geom->MakeSphere("Sphere013",Iron,0,10,0,180,0,360);
 Sphere013->SetFillColor(2);
 Sphere013->SetLineColor(2);
 top->AddNodeOverlap(Sphere013,1,new TGeoTranslation(0,410,250));

 TGeoVolume *Sphere014=geom->MakeSphere("Sphere014",Iron,0,10,0,180,0,360);
 Sphere014->SetFillColor(2);
 Sphere014->SetLineColor(2);
 top->AddNodeOverlap(Sphere014,1,new TGeoTranslation(290,-290,250));

 TGeoVolume *Sphere015=geom->MakeSphere("Sphere015",Iron,0,10,0,180,0,360);
 Sphere015->SetFillColor(2);
 Sphere015->SetLineColor(2);
 top->AddNodeOverlap(Sphere015,1,new TGeoTranslation(-290,290,250));

 TGeoVolume *Sphere016=geom->MakeSphere("Sphere016",Iron,0,10,0,180,0,360);
 Sphere016->SetFillColor(2);
 Sphere016->SetLineColor(2);
 top->AddNodeOverlap(Sphere016,1,new TGeoTranslation(410,0,250));

 TGeoVolume *Sphere017=geom->MakeSphere("Sphere017",Iron,0,10,0,180,0,360);
 Sphere017->SetFillColor(2);
 Sphere017->SetLineColor(2);
 top->AddNodeOverlap(Sphere017,1,new TGeoTranslation(-410,0,250));

 TGeoVolume *Box010=geom->MakeBox("Box010",Copper,10,10,3);
 Box010->SetFillColor(3);
 Box010->SetLineColor(3);
 top->AddNodeOverlap(Box1,1,new TGeoCombiTrans(-120,-280,300,new TGeoRotation("Box010",0,0,70)));

 TGeoVolume *Box011=geom->MakeBox("Box011",Copper,10,10,3);
 Box011->SetFillColor(3);
 Box011->SetLineColor(3);
 top->AddNodeOverlap(Box011,1,new TGeoCombiTrans(120,280,300,new TGeoRotation("Box011",0,0,70)));

 TGeoVolume *Box012=geom->MakeBox("Box012",Copper,10,10,3);
 Box012->SetFillColor(3);
 Box012->SetLineColor(3);
 top->AddNodeOverlap(Box012,1,new TGeoCombiTrans(120,-280,300,new TGeoRotation("Box012",0,0,30)));

 TGeoVolume *Box013=geom->MakeBox("Box013",Copper,10,10,3);
 Box013->SetFillColor(3);
 Box013->SetLineColor(3);
 top->AddNodeOverlap(Box013,1,new TGeoCombiTrans(-120,280,300,new TGeoRotation("Box013",0,0,30)));

 TGeoVolume *Box014=geom->MakeBox("Box010",Copper,10,10,3);
 Box014->SetFillColor(3);
 Box014->SetLineColor(3);
 top->AddNodeOverlap(Box014,1,new TGeoCombiTrans(270,-120,300,new TGeoRotation("Box014",0,0,70)));

 TGeoVolume *Box015=geom->MakeBox("Box015",Copper,10,10,3);
 Box015->SetFillColor(3);
 Box015->SetLineColor(3);
 top->AddNodeOverlap(Box015,1,new TGeoCombiTrans(-270,120,300,new TGeoRotation("Box015",0,0,70)));

 TGeoVolume *Box016=geom->MakeBox("Box016",Copper,10,10,3);
 Box016->SetFillColor(3);
 Box016->SetLineColor(3);
 top->AddNodeOverlap(Box016,1,new TGeoCombiTrans(270,100,300,new TGeoRotation("Box016",0,0,30)));

 TGeoVolume *Box017=geom->MakeBox("Box017",Copper,10,10,3);
 Box017->SetFillColor(3);
 Box017->SetLineColor(3);
 top->AddNodeOverlap(Box017,1,new TGeoCombiTrans(-270,-120,300,new TGeoRotation("Box017",0,0,30)));

for(int i=1;i<=8;i++){
 TGeoVolume *Torus0101=geom->MakeTorus("Torus0101",Iron,300,70,40,45*i-4,8);
 Torus0101->SetFillColor(18);
 Torus0101->SetLineColor(18);
 top->AddNodeOverlap(Torus0101,1,new TGeoTranslation(0,0,250));

 TGeoVolume *Tubs0101=geom->MakeTubs("Line",Iron,0,400,5,45*i-1,45*i+1);
 Tubs0101->SetFillColor(18);
 Tubs0101->SetLineColor(18);
 top->AddNodeOverlap(Tubs0101,1,new TGeoTranslation(0,0,250));
}

 Cone31->SetFillColor(38);
 top->SetVisibility(0);
 geom->CloseGeometry();

 top->Draw("ogl");

}
