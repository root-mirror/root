// Builds a polar graph in a square Canvas.

void macro3(){
    auto c = new TCanvas("myCanvas","myCanvas",600,600);
    double rmin = 0.;
    double rmax = TMath::Pi()*6.;
    const int npoints = 1000;
    double r[npoints];
    double theta[npoints];
    for (int ipt = 0; ipt < npoints; ipt++) {
        theta[ipt] = ipt*(rmax-rmin)/npoints+rmin;
        r[ipt] = TMath::Sin(theta[ipt]);
    }
    auto grP1 = new TGraphPolar(npoints,theta,r);
    grP1->SetTitle("A Fan");
    grP1->SetLineWidth(3);
    grP1->SetLineColor(2);
    grP1->Draw("L");
    gPad->Update();
    grP1->GetPolargram()->SetToRadian();
}
