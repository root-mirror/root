void tprofile2poly_tprofile2d_sim1()
{
  // CREATE PLOT STURCTURES
   auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
   auto TP2D  = new TProfile2D("hprof2d","Profile of pz versus px and py",40,-4,4,40,-4,4,0,20);
   auto TP2P = new TProfile2Poly();

   TP2P->SetName("mine");
   TP2P->SetTitle("mine");


   // create the grid in such a way that it is the same as hprof2d bins
   float minx = -4; float maxx = 4;
   float miny = -4; float maxy = 4;
   float binsz = 0.2;

   for(float i=minx; i<maxx; i+=binsz){
     for(float j=miny; j<maxy; j+=binsz){
       TP2P->AddBin(i, j, i+binsz, j+binsz);
     }
   }

   // ADD EVENTS TO PLOT
   c1->Divide(2,1);
   int value = 1;

   Float_t px, py, pz;
   for ( Int_t i=0; i<30000; i++) {
      gRandom->Rannor(px,py);
      value = px*px + py*py;
      TP2D->Fill(px,py,value);
      TP2P->Fill(px,py,value);
   }
   c1->cd(1);
   TP2D->Draw("COLZ");

   c1->cd(2);
   TP2P->Draw("COLZ");
}
