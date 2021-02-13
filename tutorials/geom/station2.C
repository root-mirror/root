/// \file
/// \ingroup tutorial_geom
/// Drawing a space station (version 2), using ROOT geometry class.
///
/// Reviewed by Sunman Kim (sunman98@hanmail.net)
/// Supervisor: Prof. Inkyu Park (icpark@physics.uos.ac.kr)
///
/// How to run: `.x station2.C` in ROOT terminal, then use OpenGL
///
/// This macro was created for the evaluation of Computational Physics course in 2006.
/// We thank to Prof. Inkyu Park for his special lecture on ROOT and to all of ROOT team
///
/// \image html geom_station2.png width=800px
/// \macro_code
///
/// \author Dong Ryeol Lee (leedr2580@hanmail.net), Dept. of Physics, Univ. of Seoul

#include "TGeoManager.h"

void station2()
{


   TGeoManager *geom = new TGeoManager("geom","Space Station");

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


   TGeoVolume *b1=geom->MakeBox("b1",Iron,1,2,3);
   b1->SetLineColor(37);
   b1->SetFillColor(37);

   TGeoVolume *b2=geom->MakeBox("b2",Iron,1,2,3);
   b2->SetLineColor(37);
   b2->SetFillColor(37);

   TGeoVolume *b12=geom->MakeBox("b12",Iron,1,2,3);
   b12->SetLineColor(37);
   b12->SetFillColor(37);

   TGeoVolume *b22=geom->MakeBox("b22",Iron,1,2,3);
   b22->SetLineColor(37);
   b22->SetFillColor(37);

   TGeoVolume *b13=geom->MakeBox("b13",Iron,1,2,3);
   b13->SetLineColor(37);
   b13->SetFillColor(37);

   TGeoVolume *b23=geom->MakeBox("b23",Iron,1,2,3);
   b23->SetLineColor(37);
   b23->SetFillColor(37);

   TGeoVolume *b14=geom->MakeBox("b14",Iron,1,2,3);
   b14->SetLineColor(37);
   b14->SetFillColor(37);

   TGeoVolume *b24=geom->MakeBox("b24",Iron,1,2,3);
   b24->SetLineColor(37);
   b24->SetFillColor(37);

   TGeoVolume *b3=geom->MakeBox("b3",Iron,35,1,1);
   b3->SetLineColor(17);
   b3->SetFillColor(17);

   TGeoVolume *b4=geom->MakeBox("b4",Iron,35,1,1);
   b4->SetLineColor(17);
   b4->SetFillColor(17);

   TGeoVolume *b31=geom->MakeBox("b31",Iron,5,5,1);
   b31->SetLineColor(38);
   b31->SetFillColor(38);

   TGeoVolume *b41=geom->MakeBox("b41",Iron,5,5,1);
   b41->SetLineColor(38);
   b41->SetFillColor(38);

   TGeoVolume *b32=geom->MakeBox("b32",Iron,5,5,1);
   b32->SetLineColor(38);
   b32->SetFillColor(38);

   TGeoVolume *b42=geom->MakeBox("b42",Iron,5,5,1);
   b42->SetLineColor(38);
   b42->SetFillColor(38);

   TGeoVolume *b33=geom->MakeBox("b33",Iron,5,5,1);
   b33->SetLineColor(38);
   b33->SetFillColor(38);

   TGeoVolume *b43=geom->MakeBox("b43",Iron,5,5,1);
   b43->SetLineColor(38);
   b43->SetFillColor(38);

   TGeoVolume *b34=geom->MakeBox("b34",Iron,6,7,1);
   b34->SetLineColor(38);
   b34->SetFillColor(38);

   TGeoVolume *b44=geom->MakeBox("b44",Iron,6,7,1);
   b44->SetLineColor(38);
   b44->SetFillColor(38);

   TGeoVolume *b35=geom->MakeBox("b35",Iron,6,7,1);
   b35->SetLineColor(38);
   b35->SetFillColor(38);

   TGeoVolume *b45=geom->MakeBox("b45",Iron,6,7,1);
   b45->SetLineColor(38);
   b45->SetFillColor(38);

   TGeoVolume *b36=geom->MakeBox("b36",Iron,6,7,1);
   b36->SetLineColor(38);
   b36->SetFillColor(38);

   TGeoVolume *b46=geom->MakeBox("b46",Iron,6,7,1);
   b46->SetLineColor(38);
   b46->SetFillColor(38);

   TGeoVolume *b37=geom->MakeBox("b37",Iron,5,5,1);
   b37->SetLineColor(38);
   b37->SetFillColor(38);

   TGeoVolume *b47=geom->MakeBox("b47",Iron,5,5,1);
   b47->SetLineColor(38);
   b47->SetFillColor(38);

   TGeoVolume *b38=geom->MakeBox("b38",Iron,5,5,1);
   b38->SetLineColor(38);
   b38->SetFillColor(38);

   TGeoVolume *b48=geom->MakeBox("b48",Iron,5,5,1);
   b48->SetLineColor(38);
   b48->SetFillColor(38);

   TGeoVolume *b39=geom->MakeBox("b39",Iron,5,5,1);
   b39->SetLineColor(38);
   b39->SetFillColor(38);

   TGeoVolume *b49=geom->MakeBox("b49",Iron,5,5,1);
   b49->SetLineColor(38);
   b49->SetFillColor(38);

   TGeoVolume *b310=geom->MakeBox("b310",Iron,6,7,1);
   b310->SetLineColor(38);
   b310->SetFillColor(38);

   TGeoVolume *b410=geom->MakeBox("b410",Iron,6,7,1);
   b410->SetLineColor(38);
   b410->SetFillColor(38);

   TGeoVolume *b311=geom->MakeBox("b311",Iron,6,7,1);
   b311->SetLineColor(38);
   b311->SetFillColor(38);

   TGeoVolume *b411=geom->MakeBox("b411",Iron,6,7,1);
   b411->SetLineColor(38);
   b411->SetFillColor(38);

   TGeoVolume *b312=geom->MakeBox("b312",Iron,6,7,1);
   b312->SetLineColor(38);
   b312->SetFillColor(38);

   TGeoVolume *b412=geom->MakeBox("b412",Iron,6,7,1);
   b412->SetLineColor(38);
   b412->SetFillColor(38);

   TGeoVolume *b5=geom->MakeBox("b5",Iron,43,1,1);
   b5->SetLineColor(17);
   b5->SetFillColor(17);

   TGeoVolume *b6=geom->MakeBox("b6",Iron,43,1,1);
   b6->SetLineColor(17);
   b6->SetFillColor(17);

   TGeoVolume *b51=geom->MakeBox("b51",Iron,5,5,1);
   b51->SetLineColor(38);
   b51->SetFillColor(38);

   TGeoVolume *b61=geom->MakeBox("b61",Iron,5,5,1);
   b61->SetLineColor(38);
   b61->SetFillColor(38);

   TGeoVolume *b52=geom->MakeBox("b52",Iron,5,5,1);
   b52->SetLineColor(38);
   b52->SetFillColor(38);

   TGeoVolume *b62=geom->MakeBox("b62",Iron,5,5,1);
   b62->SetLineColor(38);
   b62->SetFillColor(38);

   TGeoVolume *b53=geom->MakeBox("b53",Iron,5,5,1);
   b53->SetLineColor(38);
   b53->SetFillColor(38);

   TGeoVolume *b63=geom->MakeBox("b63",Iron,5,5,1);
   b63->SetLineColor(38);
   b63->SetFillColor(38);

   TGeoVolume *b54=geom->MakeBox("b54",Iron,6,7,1);
   b54->SetLineColor(38);
   b54->SetFillColor(38);

   TGeoVolume *b64=geom->MakeBox("b64",Iron,6,7,1);
   b64->SetLineColor(38);
   b64->SetFillColor(38);

   TGeoVolume *b55=geom->MakeBox("b55",Iron,6,7,1);
   b55->SetLineColor(38);
   b55->SetFillColor(38);

   TGeoVolume *b65=geom->MakeBox("b65",Iron,6,7,1);
   b65->SetLineColor(38);
   b65->SetFillColor(38);

   TGeoVolume *b56=geom->MakeBox("b56",Iron,6,7,1);
   b56->SetLineColor(38);
   b56->SetFillColor(38);

   TGeoVolume *b66=geom->MakeBox("b66",Iron,6,7,1);
   b66->SetLineColor(38);
   b66->SetFillColor(38);

   TGeoVolume *b57=geom->MakeBox("b57",Iron,5,5,1);
   b57->SetLineColor(38);
   b57->SetFillColor(38);

   TGeoVolume *b67=geom->MakeBox("b67",Iron,5,5,1);
   b67->SetLineColor(38);
   b67->SetFillColor(38);

   TGeoVolume *b58=geom->MakeBox("b58",Iron,5,5,1);
   b58->SetLineColor(38);
   b58->SetFillColor(38);

   TGeoVolume *b68=geom->MakeBox("b68",Iron,5,5,1);
   b68->SetLineColor(38);
   b68->SetFillColor(38);

   TGeoVolume *b59=geom->MakeBox("b59",Iron,5,5,1);
   b59->SetLineColor(38);
   b59->SetFillColor(38);

   TGeoVolume *b69=geom->MakeBox("b69",Iron,5,5,1);
   b69->SetLineColor(38);
   b69->SetFillColor(38);

   TGeoVolume *b510=geom->MakeBox("b510",Iron,6,7,1);
   b510->SetLineColor(38);
   b510->SetFillColor(38);

   TGeoVolume *b610=geom->MakeBox("b610",Iron,6,7,1);
   b610->SetLineColor(38);
   b610->SetFillColor(38);

   TGeoVolume *b511=geom->MakeBox("b511",Iron,6,7,1);
   b511->SetLineColor(38);
   b511->SetFillColor(38);

   TGeoVolume *b611=geom->MakeBox("b611",Iron,6,7,1);
   b611->SetLineColor(38);
   b611->SetFillColor(38);

   TGeoVolume *b512=geom->MakeBox("b512",Iron,6,7,1);
   b512->SetLineColor(38);
   b512->SetFillColor(38);

   TGeoVolume *b612=geom->MakeBox("b612",Iron,6,7,1);
   b612->SetLineColor(38);
   b612->SetFillColor(38);

   TGeoVolume *b513=geom->MakeBox("b513",Iron,6,7,1);
   b513->SetLineColor(38);
   b513->SetFillColor(38);

   TGeoVolume *b613=geom->MakeBox("b613",Iron,6,7,1);
   b613->SetLineColor(38);
   b613->SetFillColor(38);

   TGeoVolume *b514=geom->MakeBox("b514",Iron,6,7,1);
   b514->SetLineColor(38);
   b514->SetFillColor(38);

   TGeoVolume *b614=geom->MakeBox("b614",Iron,6,7,1);
   b614->SetLineColor(38);
   b614->SetFillColor(38);

   TGeoVolume *b7=geom->MakeBox("b7",Iron,5,8,15);
   b7->SetLineColor(17);
   b7->SetFillColor(17);

   TGeoVolume *b71=geom->MakeBox("b71",Iron,1,34,1);
   b71->SetLineColor(17);
   b71->SetFillColor(17);

   TGeoVolume *b81=geom->MakeBox("b81",Iron,1,34,1);
   b81->SetLineColor(17);
   b81->SetFillColor(17);

   TGeoVolume *b72=geom->MakeBox("b72",Iron,1,6,11);
   b72->SetLineColor(18);
   b72->SetFillColor(18);

   TGeoVolume *b82=geom->MakeBox("b82",Iron,1,6,11);
   b82->SetLineColor(18);
   b82->SetFillColor(18);

   TGeoVolume *b73=geom->MakeBox("b73",Iron,1,6,11);
   b73->SetLineColor(12);
   b73->SetFillColor(12);

   TGeoVolume *b83=geom->MakeBox("b83",Iron,1,6,11);
   b83->SetLineColor(12);
   b83->SetFillColor(12);

   TGeoVolume *b74=geom->MakeBox("b74",Iron,1,6,11);
   b74->SetLineColor(18);
   b74->SetFillColor(18);

   TGeoVolume *b84=geom->MakeBox("b84",Iron,1,6,11);
   b84->SetLineColor(18);
   b84->SetFillColor(18);

   TGeoVolume *b75=geom->MakeBox("b75",Iron,1,6,11);
   b75->SetLineColor(12);
   b75->SetFillColor(12);

   TGeoVolume *b85=geom->MakeBox("b85",Iron,1,6,11);
   b85->SetLineColor(12);
   b85->SetFillColor(12);

   TGeoVolume *b76=geom->MakeBox("b76",Iron,1,6,11);
   b76->SetLineColor(18);
   b76->SetFillColor(18);

   TGeoVolume *b86=geom->MakeBox("b86",Iron,1,6,11);
   b86->SetLineColor(18);
   b86->SetFillColor(18);

   TGeoVolume *b9=geom->MakeBox("b9",Iron,2,7,5);
   b9->SetLineColor(17);
   b9->SetFillColor(17);

   TGeoVolume *b10=geom->MakeBox("b10",Iron,1,20,1);
   b10->SetLineColor(17);
   b10->SetFillColor(17);

   TGeoVolume *b111=geom->MakeBox("b111",Iron,1,20,1);
   b111->SetLineColor(17);
   b111->SetFillColor(17);

   TGeoVolume *b121=geom->MakeBox("b121",Iron,1,20,1);
   b121->SetLineColor(17);
   b121->SetFillColor(17);

   TGeoVolume *b131=geom->MakeBox("b131",Iron,1,20,1);
   b131->SetLineColor(17);
   b131->SetFillColor(17);

   TGeoVolume *n1=geom->MakeBox("n1",Iron,1,6,1);
   n1->SetLineColor(12);
   n1->SetFillColor(12);

   TGeoVolume *n2=geom->MakeBox("n2",Iron,1,6,1);
   n2->SetLineColor(12);
   n2->SetFillColor(12);

   TGeoVolume *n3=geom->MakeBox("n3",Iron,1,6,1);
   n3->SetLineColor(12);
   n3->SetFillColor(12);

   TGeoVolume *n4=geom->MakeBox("n4",Iron,1,6,1);
   n4->SetLineColor(12);
   n4->SetFillColor(12);

   TGeoVolume *n5=geom->MakeBox("n5",Iron,1,1,3);
   n5->SetLineColor(12);
   n5->SetFillColor(12);

   TGeoVolume *n6=geom->MakeBox("n6",Iron,1,6,1);
   n6->SetLineColor(12);
   n6->SetFillColor(12);

   TGeoVolume *n7=geom->MakeBox("n7",Iron,1,6,1);
   n7->SetLineColor(12);
   n7->SetFillColor(12);

   TGeoVolume *n8=geom->MakeBox("n8",Iron,1,1,3);
   n8->SetLineColor(12);
   n8->SetFillColor(12);

   TGeoVolume *n9=geom->MakeBox("n9",Iron,1,6,1);
   n9->SetLineColor(12);
   n9->SetFillColor(12);

   TGeoVolume *sp=geom->MakeSphere("sp",Iron,0,10,0,180,0,360);
   sp->SetLineColor(50);
   sp->SetFillColor(50);

   TGeoVolume *sp1=geom->MakeSphere("sp1",Iron,0,50,0,180,0,360);
   sp1->SetLineColor(9);
   sp1->SetFillColor(9);

   TGeoVolume *sp2=geom->MakeSphere("sp2",Iron,0,2,0,180,0,360);
   sp2->SetLineColor(2);
   sp2->SetFillColor(2);

   TGeoVolume *sp3=geom->MakeSphere("sp3",Iron,0,2,0,180,0,360);
   sp3->SetLineColor(4);
   sp3->SetFillColor(4);

   TGeoVolume *tbs=geom->MakeTubs("tbs",Iron,0,3,120,0,360);
   tbs->SetLineColor(10);
   tbs->SetFillColor(10);

   TGeoVolume *tbs1=geom->MakeTubs("tbs1",Iron,3,5,15,0,360);
   tbs1->SetLineColor(17);
   tbs1->SetFillColor(17);

   TGeoVolume *tbs2=geom->MakeTubs("tbs2",Iron,3,15,30,0,360);
   tbs2->SetLineColor(17);
   tbs2->SetFillColor(17);

   TGeoVolume *tbs3=geom->MakeTubs("tbs3",Iron,3,10,10,0,360);
   tbs3->SetLineColor(17);
   tbs3->SetFillColor(17);

   TGeoVolume *tbs4=geom->MakeTubs("tbs4",Iron,3,7,10,0,360);
   tbs4->SetLineColor(18);
   tbs4->SetFillColor(18);

   TGeoVolume *tbs5=geom->MakeTubs("tbs5",Iron,3,13,20,0,360);
   tbs5->SetLineColor(17);
   tbs5->SetFillColor(17);

   TGeoVolume *tbs6=geom->MakeTubs("tbs6",Iron,3,7,10,0,360);
   tbs6->SetLineColor(18);
   tbs6->SetFillColor(18);

   TGeoVolume *tbs7=geom->MakeTubs("tbs7",Iron,3,15,22,0,360);
   tbs7->SetLineColor(15);
   tbs7->SetFillColor(15);

   TGeoVolume *tbs8=geom->MakeTubs("tbs8",Iron,0,10,5,0,360);
   tbs8->SetLineColor(17);
   tbs8->SetFillColor(17);

   TGeoVolume *tbs9=geom->MakeTubs("tbs9",Iron,0,15,5,0,360);
   tbs9->SetLineColor(15);
   tbs9->SetFillColor(15);

   TGeoVolume *tbs10=geom->MakeTubs("tbs10",Iron,4,6,8,0,360);
   tbs10->SetLineColor(15);
   tbs10->SetFillColor(15);

   TGeoVolume *tbs11=geom->MakeTubs("tbs11",Iron,0,4,6,0,360);
   tbs11->SetLineColor(17);
   tbs11->SetFillColor(17);

   TGeoVolume *tbs12=geom->MakeTubs("tbs12",Iron,0,4,6,0,360);
   tbs12->SetLineColor(17);
   tbs12->SetFillColor(17);

   TGeoVolume *tbs13=geom->MakeTubs("tbs13",Iron,1.7,3.7,1,-90,180);
   tbs13->SetLineColor(12);
   tbs13->SetFillColor(12);

   TGeoVolume *tbs14=geom->MakeTubs("tbs14",Iron,1.7,3.7,1,90,0);
   tbs14->SetLineColor(12);
   tbs14->SetFillColor(12);

   TGeoVolume *tbs15=geom->MakeTubs("tbs15",Iron,14,15.2,9,0,360);
   tbs15->SetLineColor(10);
   tbs15->SetFillColor(10);

   TGeoVolume *tbs16=geom->MakeTubs("tbs16",Iron,14,15.2,12,0,360);
   tbs16->SetLineColor(10);
   tbs16->SetFillColor(10);

   TGeoVolume *tbs18=geom->MakeTubs("tbs18",Iron,14,15.2,9,80,100);
   tbs18->SetLineColor(13);
   tbs18->SetFillColor(13);

   TGeoVolume *tbs19=geom->MakeTubs("tbs19",Iron,14,15.2,9,80,100);
   tbs19->SetLineColor(13);
   tbs19->SetFillColor(13);

   TGeoVolume *tbs20=geom->MakeTubs("tbs20",Iron,12,13.2,14,80,100);
   tbs20->SetLineColor(13);
   tbs20->SetFillColor(13);

   TGeoVolume *tbs21=geom->MakeTubs("tbs21",Iron,12,13.2,14,80,100);
   tbs21->SetLineColor(13);
   tbs21->SetFillColor(13);

   TGeoVolume *tbs22=geom->MakeTubs("tbs22",Iron,14,15.2,12,80,100);
   tbs22->SetLineColor(13);
   tbs22->SetFillColor(13);

   TGeoVolume *tbs23=geom->MakeTubs("tbs23",Iron,14,15.2,12,80,100);
   tbs23->SetLineColor(13);
   tbs23->SetFillColor(13);


   TGeoVolume *Cone=geom->MakeCone("Cone",Copper,3,3,10,3,15);
   Cone->SetLineColor(17);
   Cone->SetFillColor(17);

   TGeoVolume *Cone1=geom->MakeCone("Cone1",Copper,3,3,5,3,15);
   Cone1->SetLineColor(17);
   Cone1->SetFillColor(17);

   TGeoVolume *Cone2=geom->MakeCone("Cone2",Copper,3,3,13,3,7);
   Cone2->SetLineColor(17);
   Cone2->SetFillColor(17);

   TGeoVolume *Cone3=geom->MakeCone("Cone3",Copper,3,3,10,3,7);
   Cone3->SetLineColor(17);
   Cone3->SetFillColor(17);

   TGeoVolume *Cone4=geom->MakeCone("Cone4",Copper,3,3,7,3,13);
   Cone4->SetLineColor(17);
   Cone4->SetFillColor(17);

   TGeoVolume *Cone5=geom->MakeCone("Cone5",Copper,3,3,15,3,7);
   Cone5->SetLineColor(15);
   Cone5->SetFillColor(15);

   TGeoVolume *Cone6=geom->MakeCone("Cone6",Copper,8,0,8,0,8);
   Cone6->SetLineColor(17);
   Cone6->SetFillColor(17);

   TGeoVolume *Cone7=geom->MakeCone("Cone7",Copper,1,3,5,3,6);
   Cone7->SetLineColor(18);
   Cone7->SetFillColor(18);

   TGeoVolume *Cone8=geom->MakeCone("Cone8",Copper,3,3,15,3,7);
   Cone8->SetLineColor(15);
   Cone8->SetFillColor(15);

   TGeoVolume *Cone9=geom->MakeCone("Cone9",Copper,1,3,5,3,6);
   Cone9->SetLineColor(12);
   Cone9->SetFillColor(12);

   TGeoVolume *Cone10=geom->MakeCone("Cone10",Copper,1,3,5,3,6);
   Cone10->SetLineColor(12);
   Cone10->SetFillColor(12);

   TGeoVolume *Cone11=geom->MakeCone("Cone11",Copper,1,3,5,3,6);
   Cone11->SetLineColor(14);
   Cone11->SetFillColor(14);

   TGeoVolume *Cone12=geom->MakeCone("Cone12",Copper,1,3,5,3,6);
   Cone12->SetLineColor(14);
   Cone12->SetFillColor(14);

   TGeoVolume *a1=geom->MakeBox("a1",Iron,2,1,2);
   a1->SetLineColor(10);
   a1->SetFillColor(10);
   top->AddNodeOverlap(a1,1,new TGeoCombiTrans(0,15,98, new TGeoRotation("a1",0,30,0)));

   TGeoVolume *a2=geom->MakeBox("a2",Iron,2,1,2);
   a2->SetLineColor(10);
   a2->SetFillColor(10);
   top->AddNodeOverlap(a2,1,new TGeoCombiTrans(0,15,90, new TGeoRotation("a2",0,0,0)));

   TGeoVolume *a3=geom->MakeBox("a3",Iron,2,1,2);
   a3->SetLineColor(10);
   a3->SetFillColor(10);
   top->AddNodeOverlap(a3,1,new TGeoCombiTrans(0,15,85, new TGeoRotation("a3",0,0,0)));

   TGeoVolume *a4=geom->MakeBox("a4",Iron,2,1,2);
   a4->SetLineColor(10);
   a4->SetFillColor(10);
   top->AddNodeOverlap(a4,1,new TGeoCombiTrans(3,14,76, new TGeoRotation("a4",0,0,0)));

   TGeoVolume *a5=geom->MakeBox("a5",Iron,2,1,2);
   a5->SetLineColor(10);
   a5->SetFillColor(10);
   top->AddNodeOverlap(a5,1,new TGeoCombiTrans(-7,13,75, new TGeoRotation("a5",0,0,0)));

   TGeoVolume *a6=geom->MakeBox("a6",Iron,2,1,2);
   a6->SetLineColor(10);
   a6->SetFillColor(10);
   top->AddNodeOverlap(a6,1,new TGeoCombiTrans(-7,13,71, new TGeoRotation("a6",0,0,0)));

   TGeoVolume *a7=geom->MakeBox("a7",Iron,2,1,2);
   a7->SetLineColor(10);
   a7->SetFillColor(10);
   top->AddNodeOverlap(a7,1,new TGeoCombiTrans(-6,13,66, new TGeoRotation("a7",0,40,0)));

   TGeoVolume *a8=geom->MakeBox("a8",Iron,2,1,2);
   a8->SetLineColor(10);
   a8->SetFillColor(10);
   top->AddNodeOverlap(a8,1,new TGeoCombiTrans(-7,13,60, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a9=geom->MakeBox("a9",Iron,2,1,2);
   a9->SetLineColor(10);
   a9->SetFillColor(10);
   top->AddNodeOverlap(a9,1,new TGeoCombiTrans(3,12,-1, new TGeoRotation("a9",0,0,0)));

   TGeoVolume *a10=geom->MakeBox("a12",Iron,2,1,2);
   a10->SetLineColor(10);
   a10->SetFillColor(10);
   top->AddNodeOverlap(a10,1,new TGeoCombiTrans(2,12,-6, new TGeoRotation("a5",0,0,0)));

   TGeoVolume *a11=geom->MakeBox("a11",Iron,2,1,2);
   a11->SetLineColor(10);
   a11->SetFillColor(10);
   top->AddNodeOverlap(a11,1,new TGeoCombiTrans(-3,12,-20, new TGeoRotation("a6",20,0,0)));

   TGeoVolume *a12=geom->MakeBox("a12",Iron,2,1,2);
   a12->SetLineColor(10);
   a12->SetFillColor(10);
   top->AddNodeOverlap(a12,1,new TGeoCombiTrans(-1,12,-25, new TGeoRotation("a7",0,40,0)));

   TGeoVolume *a13=geom->MakeBox("a13",Iron,2,1,2);
   a13->SetLineColor(10);
   a13->SetFillColor(10);
   top->AddNodeOverlap(a13,1,new TGeoCombiTrans(-3,12,-29, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a14=geom->MakeTubs("a14",Iron,0,1,20,0,360);
   a14->SetLineColor(36);
   a14->SetFillColor(36);
   top->AddNodeOverlap(a14,1,new TGeoCombiTrans(7.5,7.5,20, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a15=geom->MakeTubs("a15",Iron,0,1,20,0,360);
   a15->SetLineColor(36);
   a15->SetFillColor(36);
   top->AddNodeOverlap(a15,1,new TGeoCombiTrans(-7.5,7.5,20, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a16=geom->MakeTubs("a16",Iron,0,1,20,0,360);
   a16->SetLineColor(36);
   a16->SetFillColor(36);
   top->AddNodeOverlap(a16,1,new TGeoCombiTrans(7.5,-7.5,20, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a17=geom->MakeTubs("a17",Iron,0,1,20,0,360);
   a17->SetLineColor(36);
   a17->SetFillColor(36);
   top->AddNodeOverlap(a17,1,new TGeoCombiTrans(-7.5,-7.5,20, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a18=geom->MakeTubs("a18",Iron,0,1,20,0,360);
   a18->SetLineColor(36);
   a18->SetFillColor(36);
   top->AddNodeOverlap(a18,1,new TGeoCombiTrans(7.5,7.5,-50, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a19=geom->MakeTubs("a19",Iron,0,1,20,0,360);
   a19->SetLineColor(36);
   a19->SetFillColor(36);
   top->AddNodeOverlap(a19,1,new TGeoCombiTrans(-7.5,7.5,-50, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a20=geom->MakeTubs("a20",Iron,0,1,20,0,360);
   a20->SetLineColor(36);
   a20->SetFillColor(36);
   top->AddNodeOverlap(a20,1,new TGeoCombiTrans(7.5,-7.5,-50, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a21=geom->MakeTubs("a21",Iron,0,1,20,0,360);
   a21->SetLineColor(36);
   a21->SetFillColor(36);
   top->AddNodeOverlap(a21,1,new TGeoCombiTrans(-7.5,-7.5,-50, new TGeoRotation("a8",0,0,0)));

   TGeoVolume *a22=geom->MakeTubs("a22",Iron,3,4,3,0,360);
   a22->SetLineColor(10);
   a22->SetFillColor(10);
   top->AddNodeOverlap(a22,1,new TGeoCombiTrans(14,6,97, new TGeoRotation("a22",110,90,0)));

   TGeoVolume *a23=geom->MakeTubs("a23",Iron,3,4,3,0,360);
   a23->SetLineColor(14);
   a23->SetFillColor(14);
   top->AddNodeOverlap(a23,1,new TGeoCombiTrans(0,-7,14, new TGeoRotation("a22",180,90,0)));


   TGeoVolume *Cone15=geom->MakeCone("Cone15",Copper,1,3,4,4,5);
   Cone15->SetLineColor(14);
   Cone15->SetFillColor(14);
   top->AddNodeOverlap(Cone15,1,new TGeoCombiTrans(0,-11,14, new TGeoRotation("a23",0,90,0)));


   TGeoVolume *a24=geom->MakeTubs("a24",Iron,3,4,3,0,360);
   a24->SetLineColor(14);
   a24->SetFillColor(14);
   top->AddNodeOverlap(a24,1,new TGeoCombiTrans(0,-7,-46, new TGeoRotation("a23",180,90,0)));

   TGeoVolume *a25=geom->MakeTubs("a25",Iron,3,5,8,0,360);
   a25->SetLineColor(18);
   a25->SetFillColor(18);
   top->AddNodeOverlap(a25,1,new TGeoCombiTrans(0,-20,-46, new TGeoRotation("a23",180,90,0)));

   TGeoVolume *Cone13=geom->MakeCone("Cone13",Copper,1,3,4,4,5);
   Cone13->SetLineColor(14);
   Cone13->SetFillColor(14);
   top->AddNodeOverlap(Cone13,1,new TGeoCombiTrans(0,-11,-46, new TGeoRotation("a23",0,90,0)));

   TGeoVolume *Cone14=geom->MakeCone("Cone14",Copper,1,3,4,4,5);
   Cone14->SetLineColor(14);
   Cone14->SetFillColor(14);
   top->AddNodeOverlap(Cone14,1,new TGeoCombiTrans(0,-29,-46, new TGeoRotation("a23",0,270,0)));

   TGeoVolume *sp4=geom->MakeSphere("sp4",Iron,0,4,0,180,0,360);
   sp4->SetLineColor(10);
   sp4->SetFillColor(10);
   top->AddNodeOverlap(sp4,1,new TGeoCombiTrans(0,-32,-46, new TGeoRotation("a23",0,0,0)));

   TGeoVolume *Cone16=geom->MakeCone("Cone16",Copper,1,3,4,4,5);
   Cone16->SetLineColor(14);
   Cone16->SetFillColor(14);
   top->AddNodeOverlap(Cone16,1,new TGeoCombiTrans(-1,-35,-46, new TGeoRotation("a23",-30,80,0)));

   TGeoVolume *a26=geom->MakeTubs("a26",Iron,3,5,12,0,360);
   a26->SetLineColor(18);
   a26->SetFillColor(18);
   top->AddNodeOverlap(a26,1,new TGeoCombiTrans(-7.5,-46,-43.5, new TGeoRotation("a23",-30,80,0)));

   TGeoVolume *Cone17=geom->MakeCone("Cone17",Copper,1,3,4,4,5);
   Cone17->SetLineColor(14);
   Cone17->SetFillColor(14);
   top->AddNodeOverlap(Cone17,1,new TGeoCombiTrans(-13.7,-57,-41.2, new TGeoRotation("a23",-30,260,0)));

   TGeoVolume *a27=geom->MakeTubs("a27",Iron,4,6,12,0,360);
   a27->SetLineColor(18);
   a27->SetFillColor(18);
   top->AddNodeOverlap(a27,1,new TGeoCombiTrans(23.2,0,31, new TGeoRotation("a23",90,90,0)));

   TGeoVolume *Cone18=geom->MakeCone("Cone18",Copper,1,3,5,3,6);
   Cone18->SetLineColor(14);
   Cone18->SetFillColor(14);
   top->AddNodeOverlap(Cone18,1,new TGeoCombiTrans(36,0,31, new TGeoRotation("c34",270,90,0)));





   char nBlocks[50];
   int i=1;
   int N=0;
   int f=0;
   TGeoVolume *mBlock;

   f=0;
   while (f<4){
   i=0;
   while (i<30){
        sprintf(nBlocks,"f%d_bg%d",f,N++);
        mBlock = geom->MakeBox(nBlocks, Copper,2,1,3);
        mBlock->SetLineColor(46);
        top->AddNodeOverlap(mBlock,1,new TGeoCombiTrans(15+(i*5),75,-65+(f*7), new TGeoRotation("z",0,0,0)));
        i++;
     }
        f++;
     }

       TGeoVolume *mBlock1;

   f=0;
   while (f<4){
   i=0;
   while (i<30){
        sprintf(nBlocks,"f%d_bg%d",f,N++);
        mBlock1 = geom->MakeBox(nBlocks, Copper,2,1,3);
        mBlock1->SetLineColor(46);
        top->AddNodeOverlap(mBlock1,1,new TGeoCombiTrans(14+(i*5),75,-100+(f*7), new TGeoRotation("z",0,0,0)));
        i++;
     }
        f++;
     }

       TGeoVolume *mBlock2;

   f=0;
   while (f<4){
   i=0;
   while (i<30){
        sprintf(nBlocks,"f%d_bg%d",f,N++);
        mBlock2 = geom->MakeBox(nBlocks, Copper,2,1,3);
        mBlock2->SetLineColor(46);
        top->AddNodeOverlap(mBlock2,1,new TGeoCombiTrans(-160+(i*5),75,-75+(f*7), new TGeoRotation("z",0,0,0)));
        i++;
     }
        f++;
     }

       TGeoVolume *mBlock3;

   f=0;
   while (f<4){
   i=0;
   while (i<30){
        sprintf(nBlocks,"f%d_bg%d",f,N++);
        mBlock3 = geom->MakeBox(nBlocks, Copper,2,1,3);
        mBlock3->SetLineColor(46);
        top->AddNodeOverlap(mBlock3,1,new TGeoCombiTrans(-160+(i*5),75,-110+(f*7), new TGeoRotation("z",0,0,0)));
        i++;
     }
        f++;
     }


   top->AddNodeOverlap(b1,1,new TGeoCombiTrans(5,5,130, new TGeoRotation("b1",0,0,-45)));
   top->AddNodeOverlap(b2,1,new TGeoCombiTrans(-5,-5,130, new TGeoRotation("b2",0,0,-45)));
   top->AddNodeOverlap(b12,1,new TGeoCombiTrans(8.2,8.2,130, new TGeoRotation("b12",0,0,-45)));
   top->AddNodeOverlap(b22,1,new TGeoCombiTrans(-8.2,-8.2,130, new TGeoRotation("b22",0,0,-45)));
   top->AddNodeOverlap(b13,1,new TGeoCombiTrans(11.4,11.4,130, new TGeoRotation("b13",0,0,-45)));
   top->AddNodeOverlap(b23,1,new TGeoCombiTrans(-11.4,-11.4,130, new TGeoRotation("b23",0,0,-45)));
   top->AddNodeOverlap(b14,1,new TGeoCombiTrans(14.6,14.6,130, new TGeoRotation("b14",0,0,-45)));
   top->AddNodeOverlap(b24,1,new TGeoCombiTrans(-14.6,-14.6,130, new TGeoRotation("b24",0,0,-45)));
   top->AddNodeOverlap(b3,1,new TGeoCombiTrans(50,0,71, new TGeoRotation("b3",0,-30,0)));
   top->AddNodeOverlap(b4,1,new TGeoCombiTrans(-50,0,71, new TGeoRotation("b4",0,-30,0)));
   top->AddNodeOverlap(b31,1,new TGeoCombiTrans(20,5,68, new TGeoRotation("b31",0,-30,0)));
   top->AddNodeOverlap(b41,1,new TGeoCombiTrans(-20,5,68, new TGeoRotation("b41",0,-30,0)));
   top->AddNodeOverlap(b32,1,new TGeoCombiTrans(31,5,68, new TGeoRotation("b32",0,-30,0)));
   top->AddNodeOverlap(b42,1,new TGeoCombiTrans(-31,5,68, new TGeoRotation("b42",0,-30,0)));
   top->AddNodeOverlap(b33,1,new TGeoCombiTrans(42,5,68, new TGeoRotation("b33",0,-30,0)));
   top->AddNodeOverlap(b43,1,new TGeoCombiTrans(-42,5,68, new TGeoRotation("b43",0,-30,0)));
   top->AddNodeOverlap(b34,1,new TGeoCombiTrans(54,7,67, new TGeoRotation("b34",0,-30,0)));
   top->AddNodeOverlap(b44,1,new TGeoCombiTrans(-54,7,67, new TGeoRotation("b44",0,-30,0)));
   top->AddNodeOverlap(b35,1,new TGeoCombiTrans(67,7,67, new TGeoRotation("b35",0,-30,0)));
   top->AddNodeOverlap(b45,1,new TGeoCombiTrans(-67,7,67, new TGeoRotation("b45",0,-30,0)));
   top->AddNodeOverlap(b36,1,new TGeoCombiTrans(80,7,67, new TGeoRotation("b36",0,-30,0)));
   top->AddNodeOverlap(b46,1,new TGeoCombiTrans(-80,7,67, new TGeoRotation("b46",0,-30,0)));
   top->AddNodeOverlap(b37,1,new TGeoCombiTrans(20,-5,74, new TGeoRotation("b37",0,-30,0)));
   top->AddNodeOverlap(b47,1,new TGeoCombiTrans(-20,-5,74, new TGeoRotation("b47",0,-30,0)));
   top->AddNodeOverlap(b38,1,new TGeoCombiTrans(31,-5,74, new TGeoRotation("b38",0,-30,0)));
   top->AddNodeOverlap(b48,1,new TGeoCombiTrans(-31,-5,74, new TGeoRotation("b48",0,-30,0)));
   top->AddNodeOverlap(b39,1,new TGeoCombiTrans(42,-5,74, new TGeoRotation("b39",0,-30,0)));
   top->AddNodeOverlap(b49,1,new TGeoCombiTrans(-42,-5,74, new TGeoRotation("b49",0,-30,0)));
   top->AddNodeOverlap(b310,1,new TGeoCombiTrans(54,-7,75, new TGeoRotation("b310",0,-30,0)));
   top->AddNodeOverlap(b410,1,new TGeoCombiTrans(-54,-7,75, new TGeoRotation("b410",0,-30,0)));
   top->AddNodeOverlap(b311,1,new TGeoCombiTrans(67,-7,75, new TGeoRotation("b311",0,-30,0)));
   top->AddNodeOverlap(b411,1,new TGeoCombiTrans(-67,-7,75, new TGeoRotation("b411",0,-30,0)));
   top->AddNodeOverlap(b312,1,new TGeoCombiTrans(80,-7,75, new TGeoRotation("b312",0,-30,0)));
   top->AddNodeOverlap(b412,1,new TGeoCombiTrans(-80,-7,75, new TGeoRotation("b412",0,-30,0)));
   top->AddNodeOverlap(b5,1,new TGeoCombiTrans(55,0,-15, new TGeoRotation("b5",0,-30,0)));
   top->AddNodeOverlap(b6,1,new TGeoCombiTrans(-55,0,-15, new TGeoRotation("b6",0,-30,0)));
   top->AddNodeOverlap(b51,1,new TGeoCombiTrans(20,5,-18, new TGeoRotation("b51",0,-30,0)));
   top->AddNodeOverlap(b61,1,new TGeoCombiTrans(-20,5,-18, new TGeoRotation("b61",0,-30,0)));
   top->AddNodeOverlap(b52,1,new TGeoCombiTrans(31,5,-18, new TGeoRotation("b52",0,-30,0)));
   top->AddNodeOverlap(b62,1,new TGeoCombiTrans(-31,5,-18, new TGeoRotation("b62",0,-30,0)));
   top->AddNodeOverlap(b53,1,new TGeoCombiTrans(42,5,-18, new TGeoRotation("b53",0,-30,0)));
   top->AddNodeOverlap(b63,1,new TGeoCombiTrans(-42,5,-18, new TGeoRotation("b63",0,-30,0)));
   top->AddNodeOverlap(b54,1,new TGeoCombiTrans(54,7,-19, new TGeoRotation("b54",0,-30,0)));
   top->AddNodeOverlap(b64,1,new TGeoCombiTrans(-54,7,-19, new TGeoRotation("b64",0,-30,0)));
   top->AddNodeOverlap(b55,1,new TGeoCombiTrans(67,7,-19, new TGeoRotation("b55",0,-30,0)));
   top->AddNodeOverlap(b65,1,new TGeoCombiTrans(-67,7,-19, new TGeoRotation("b65",0,-30,0)));
   top->AddNodeOverlap(b56,1,new TGeoCombiTrans(80,7,-19, new TGeoRotation("b56",0,-30,0)));
   top->AddNodeOverlap(b66,1,new TGeoCombiTrans(-80,7,-19, new TGeoRotation("b66",0,-30,0)));
   top->AddNodeOverlap(b514,1,new TGeoCombiTrans(93,7,-19, new TGeoRotation("b514",0,-30,0)));
   top->AddNodeOverlap(b614,1,new TGeoCombiTrans(-93,7,-19, new TGeoRotation("b614",0,-30,0)));
   top->AddNodeOverlap(b57,1,new TGeoCombiTrans(20,-5,-12, new TGeoRotation("b57",0,-30,0)));
   top->AddNodeOverlap(b67,1,new TGeoCombiTrans(-20,-5,-12, new TGeoRotation("b67",0,-30,0)));
   top->AddNodeOverlap(b58,1,new TGeoCombiTrans(31,-5,-12, new TGeoRotation("b58",0,-30,0)));
   top->AddNodeOverlap(b68,1,new TGeoCombiTrans(-31,-5,-12, new TGeoRotation("b68",0,-30,0)));
   top->AddNodeOverlap(b59,1,new TGeoCombiTrans(42,-5,-12, new TGeoRotation("b59",0,-30,0)));
   top->AddNodeOverlap(b69,1,new TGeoCombiTrans(-42,-5,-12, new TGeoRotation("b69",0,-30,0)));
   top->AddNodeOverlap(b510,1,new TGeoCombiTrans(54,-7,-11, new TGeoRotation("b510",0,-30,0)));
   top->AddNodeOverlap(b610,1,new TGeoCombiTrans(-54,-7,-11, new TGeoRotation("b610",0,-30,0)));
   top->AddNodeOverlap(b511,1,new TGeoCombiTrans(67,-7,-11, new TGeoRotation("b511",0,-30,0)));
   top->AddNodeOverlap(b611,1,new TGeoCombiTrans(-67,-7,-11, new TGeoRotation("b611",0,-30,0)));
   top->AddNodeOverlap(b512,1,new TGeoCombiTrans(80,-7,-11, new TGeoRotation("b512",0,-30,0)));
   top->AddNodeOverlap(b612,1,new TGeoCombiTrans(-80,-7,-11, new TGeoRotation("b612",0,-30,0)));
   top->AddNodeOverlap(b513,1,new TGeoCombiTrans(93,-7,-11, new TGeoRotation("b513",0,-30,0)));
   top->AddNodeOverlap(b613,1,new TGeoCombiTrans(-93,-7,-11, new TGeoRotation("b613",0,-30,0)));
   top->AddNodeOverlap(b7,1,new TGeoCombiTrans(0,40,-80, new TGeoRotation("b7",0,90,0)));
   top->AddNodeOverlap(b71,1,new TGeoCombiTrans(0,40,-38, new TGeoRotation("b71",0,90,0)));
   top->AddNodeOverlap(b81,1,new TGeoCombiTrans(0,57,-122, new TGeoRotation("b81",0,90,0)));
   top->AddNodeOverlap(b72,1,new TGeoCombiTrans(0,40,-62, new TGeoRotation("b72",0,90,0)));
   top->AddNodeOverlap(b82,1,new TGeoCombiTrans(0,57,-98, new TGeoRotation("b82",0,90,0)));
   top->AddNodeOverlap(b73,1,new TGeoCombiTrans(0,40,-49, new TGeoRotation("b73",0,90,0)));
   top->AddNodeOverlap(b83,1,new TGeoCombiTrans(0,57,-111, new TGeoRotation("b83",0,90,0)));
   top->AddNodeOverlap(b74,1,new TGeoCombiTrans(0,40,-36, new TGeoRotation("b74",0,90,0)));
   top->AddNodeOverlap(b84,1,new TGeoCombiTrans(0,57,-124, new TGeoRotation("b84",0,90,0)));
   top->AddNodeOverlap(b75,1,new TGeoCombiTrans(0,40,-23, new TGeoRotation("b75",0,90,0)));
   top->AddNodeOverlap(b85,1,new TGeoCombiTrans(0,57,-137, new TGeoRotation("b85",0,90,0)));
   top->AddNodeOverlap(b76,1,new TGeoCombiTrans(0,40,-10, new TGeoRotation("b76",0,90,0)));
   top->AddNodeOverlap(b86,1,new TGeoCombiTrans(0,57,-150, new TGeoRotation("b86",0,90,0)));
   top->AddNodeOverlap(b9,1,new TGeoCombiTrans(0,75,-80, new TGeoRotation("b9",0,90,0)));
   top->AddNodeOverlap(b10,1,new TGeoCombiTrans(12,75,-72, new TGeoRotation("b10",0,90,0)));
   top->AddNodeOverlap(b111,1,new TGeoCombiTrans(163,75,-72, new TGeoRotation("b111",0,90,0)));
   top->AddNodeOverlap(b121,1,new TGeoCombiTrans(-12,75,-82, new TGeoRotation("b121",0,90,0)));
   top->AddNodeOverlap(b131,1,new TGeoCombiTrans(-163,75,-82, new TGeoRotation("b131",0,90,0)));
   top->AddNodeOverlap(n1,1,new TGeoCombiTrans(-15,0,-97, new TGeoRotation("n1",0,0,0)));
   top->AddNodeOverlap(n2,1,new TGeoCombiTrans(-15,0,-94, new TGeoRotation("n2",0,-25,0)));
   top->AddNodeOverlap(n3,1,new TGeoCombiTrans(-15,0,-91, new TGeoRotation("n3",0,0,0)));
   top->AddNodeOverlap(n4,1,new TGeoCombiTrans(-15,0,-85, new TGeoRotation("n4",0,15,0)));
   top->AddNodeOverlap(n5,1,new TGeoCombiTrans(-15,-2,-83, new TGeoRotation("n5",0,0,0)));
   top->AddNodeOverlap(n6,1,new TGeoCombiTrans(-15,0,-81, new TGeoRotation("n6",0,-15,0)));
   top->AddNodeOverlap(n7,1,new TGeoCombiTrans(-15,0,-65, new TGeoRotation("n7",0,15,0)));
   top->AddNodeOverlap(n8,1,new TGeoCombiTrans(-15,-2,-63, new TGeoRotation("n8",0,0,0)));
   top->AddNodeOverlap(n9,1,new TGeoCombiTrans(-15,0,-61, new TGeoRotation("n9",0,-15,0)));
   top->AddNodeOverlap(sp,1,new TGeoTranslation(100,100,150));
   top->AddNodeOverlap(sp1,1,new TGeoTranslation(-100,-100,-150));
   top->AddNodeOverlap(sp2,1,new TGeoCombiTrans(0,80,-85, new TGeoRotation("sp1",0,0,0)));
   top->AddNodeOverlap(sp3,1,new TGeoCombiTrans(0,80,-75, new TGeoRotation("sp3",0,0,0)));
   top->AddNodeOverlap(tbs,1,new TGeoCombiTrans(0,0,14, new TGeoRotation("r1",0,0,0)));
   top->AddNodeOverlap(tbs1,1,new TGeoCombiTrans(0,0,118, new TGeoRotation("r2",0,0,0)));
   top->AddNodeOverlap(tbs2,1,new TGeoCombiTrans(0,0,74, new TGeoRotation("r3",0,0,0)));
   top->AddNodeOverlap(tbs3,1,new TGeoCombiTrans(0,0,34, new TGeoRotation("r4",0,0,0)));
   top->AddNodeOverlap(tbs4,1,new TGeoCombiTrans(0,0,14, new TGeoRotation("r5",0,0,0)));
   top->AddNodeOverlap(tbs5,1,new TGeoCombiTrans(0,0,-16, new TGeoRotation("r6",0,0,0)));
   top->AddNodeOverlap(tbs6,1,new TGeoCombiTrans(0,0,-46, new TGeoRotation("r7",0,0,0)));
   top->AddNodeOverlap(tbs7,1,new TGeoCombiTrans(0,0,-78, new TGeoRotation("r8",0,0,0)));
   top->AddNodeOverlap(tbs8,1,new TGeoCombiTrans(0,20,-80, new TGeoRotation("r9",0,90,0)));
   top->AddNodeOverlap(tbs9,1,new TGeoCombiTrans(20,0,-80, new TGeoRotation("r10",90,90,0)));
   top->AddNodeOverlap(tbs10,1,new TGeoCombiTrans(30,0,-80, new TGeoRotation("r11",90,90,0)));
   top->AddNodeOverlap(tbs11,1,new TGeoCombiTrans(5,75,-80, new TGeoRotation("r12",90,90,0)));
   top->AddNodeOverlap(tbs12,1,new TGeoCombiTrans(-5,75,-80, new TGeoRotation("r13",90,90,0)));
   top->AddNodeOverlap(tbs13,1,new TGeoCombiTrans(-15,-2.6,-73, new TGeoRotation("r14",90,90,90)));
   top->AddNodeOverlap(tbs14,1,new TGeoCombiTrans(-15,2.6,-73, new TGeoRotation("r15",90,90,90)));
   top->AddNodeOverlap(tbs15,1,new TGeoCombiTrans(0,0,95, new TGeoRotation("r16",0,0,0)));
   top->AddNodeOverlap(tbs16,1,new TGeoCombiTrans(0,0,-90, new TGeoRotation("r17",0,0,20)));
   top->AddNodeOverlap(tbs18,1,new TGeoCombiTrans(0,0,57, new TGeoRotation("r19",0,0,160)));
   top->AddNodeOverlap(tbs19,1,new TGeoCombiTrans(0,0,57, new TGeoRotation("r20",0,0,200)));
   top->AddNodeOverlap(tbs20,1,new TGeoCombiTrans(0,0,-15, new TGeoRotation("r21",0,0,160)));
   top->AddNodeOverlap(tbs21,1,new TGeoCombiTrans(0,0,-15, new TGeoRotation("r22",0,0,200)));
   top->AddNodeOverlap(tbs22,1,new TGeoCombiTrans(0,0,90, new TGeoRotation("r23",0,0,160)));
   top->AddNodeOverlap(tbs23,1,new TGeoCombiTrans(0,0,90, new TGeoRotation("r24",0,0,200)));
   top->AddNodeOverlap(Cone,1,new TGeoCombiTrans(0,0,41, new TGeoRotation("c1",0,0,0)));
   top->AddNodeOverlap(Cone1,1,new TGeoCombiTrans(0,0,107, new TGeoRotation("c2",0,180,0)));
   top->AddNodeOverlap(Cone2,1,new TGeoCombiTrans(0,0,7, new TGeoRotation("c3",0,0,0)));
   top->AddNodeOverlap(Cone3,1,new TGeoCombiTrans(0,0,21, new TGeoRotation("c4",0,180,0)));
   top->AddNodeOverlap(Cone4,1,new TGeoCombiTrans(0,0,-39, new TGeoRotation("c5",0,0,0)));
   top->AddNodeOverlap(Cone5,1,new TGeoCombiTrans(0,0,-53, new TGeoRotation("c5",0,0,0)));
   top->AddNodeOverlap(Cone6,1,new TGeoCombiTrans(0,63,-80, new TGeoRotation("c6",0,90,0)));
   top->AddNodeOverlap(Cone7,1,new TGeoCombiTrans(0,0,134, new TGeoRotation("c7",0,0,0)));
   top->AddNodeOverlap(Cone8,1,new TGeoCombiTrans(0,0,-103, new TGeoRotation("c8",0,180,0)));
   top->AddNodeOverlap(Cone9,1,new TGeoCombiTrans(-10,0,31, new TGeoRotation("c9",90,-90,90)));
   top->AddNodeOverlap(Cone10,1,new TGeoCombiTrans(10,0,31, new TGeoRotation("c9",-90,-90,90)));
   top->AddNodeOverlap(Cone11,1,new TGeoCombiTrans(39,0,-80, new TGeoRotation("c10",90,-90,90)));
   top->AddNodeOverlap(Cone12,1,new TGeoCombiTrans(0,0,-107, new TGeoRotation("c11",0,0,0)));


   top->SetVisibility(0);
   geom->CloseGeometry();

   top->Draw("ogl");

}
