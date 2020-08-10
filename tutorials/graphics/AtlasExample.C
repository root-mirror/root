/// \file
/// \ingroup tutorial_graphics
/// \notebook -js
/// Show how ATLAS Style looks like. It is based on a style file from BaBar.
///
/// \macro_image
/// \macro_code
///
/// \author  M.Sutton

const Int_t GMAX=864;

const int nren=3;
static const double mur[nren] = {1.0,0.25,4.0};
static const double muf[nren] = {1.0,0.25,4.0};
const unsigned int NUMPDF=41;

TGraphErrors* GetGraph(Int_t ir, Int_t ifs,Int_t icut, Int_t ipdf);
void AddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* MakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2);

void AtlasExample()
{
   gROOT->SetStyle("ATLAS");

   Int_t icol1=5;
   Int_t icol2=5;

   auto canvas = new TCanvas("canvas","single inclusive jets",50,50,600,600);
   canvas->SetLogy();

   Double_t ymin = 1.e-3;  Double_t ymax = 2e7;
   Double_t xmin = 60.00;  Double_t xmax = 3500.;
   auto frame = canvas->DrawFrame(xmin,ymin,xmax,ymax);
   frame->SetYTitle("d#sigma_{jet}/dE_{T,jet} [fb/GeV]");
   frame->SetXTitle("E_{T,jet}  [GeV]");
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->GetXaxis()->SetTitleOffset(1.4);

   const Int_t ncut=1;
   TGraphErrors *data[ncut];

   for (Int_t icut=0; icut<ncut; icut++) { // loop over cuts
      TGraphErrors *g1[nren][ncut];
      for (Int_t ir=0; ir<nren; ir++) { // loop over ren scale
         g1[ir][icut] = GetGraph(ir,ir,icut,0);
         if (!g1[ir][icut]) {
            cout << " g1 not  found " << endl;
            return;
         }
         g1[ir][icut]->SetLineColor(1);
         g1[ir][icut]->SetMarkerStyle(0);
      }

      char daname[100];
      sprintf(daname,"data_%d",icut);
      data[icut] = (TGraphErrors*)g1[0][icut]->Clone(daname);
      data[icut]->SetMarkerStyle(20);
      data[icut]->SetMarkerColor(1);

      // Just invent some data
      for (Int_t i=0; i< data[icut]->GetN(); i++) {
         Double_t x1,y1,e,dx1=0.;
         data[icut]->GetPoint(i,x1,y1);
         Double_t r1 = 0.4*(gRandom->Rndm(1)+2);
         Double_t r2 = 0.4*(gRandom->Rndm(1)+2);
         Double_t y;
         if (icut==0) y = r1*y1+r1*r2*r2*x1/50000.;
         else         y = r1*y1;
         e = sqrt(y*1000)/200;
         data[icut]->SetPoint(i, x1,y);
         data[icut]->SetPointError(i,dx1,e);
      }

      TGraphAsymmErrors* scale[ncut];
      TGraphAsymmErrors* scalepdf[ncut];

      scale[icut] =  MakeBand(g1[0][icut],g1[1][icut],g1[2][icut]);
      scalepdf[icut]=(TGraphAsymmErrors* ) scale[icut]->Clone("scalepdf");

      TGraphErrors *gpdf[NUMPDF][ncut];
      for (Int_t ipdf=0; ipdf<NUMPDF; ipdf++) {
         gpdf[ipdf][icut]= GetGraph(0,0,icut,ipdf);
         if (!gpdf[ipdf][icut]) {
            cout << " gpdf not  found " << endl;
            return;
         }
         gpdf[ipdf][icut]->SetLineColor(2);
         gpdf[ipdf][icut]->SetLineStyle(1);
         gpdf[ipdf][icut]->SetMarkerStyle(0);
         AddtoBand(gpdf[ipdf][icut],scalepdf[icut]);
      }

      scalepdf[icut]->SetFillColor(icol2);
      scalepdf[icut]->Draw("zE2");
      scalepdf[icut]->SetLineWidth(3);
      scale[icut]->SetFillColor(icol1);
      scale[icut]->Draw("zE2");
      g1[0][icut]->SetLineWidth(3);
      g1[0][icut]->Draw("z");
      data[icut]->Draw("P");
   }

   auto t = new TLatex; t->SetNDC();
   t->DrawLatex(0.3,  0.85, "#sqrt{s}= 14 TeV");
   t->DrawLatex(0.57, 0.85, "|#eta_{jet}|<0.5");

   auto l = new TLegend(0.45,0.65,0.8,0.8,"","NDC");
   l->SetBorderSize(0.);
   l->SetTextFont(42);
   l->AddEntry("data_0", "Data 2009", "ep");
   l->AddEntry("scalepdf", "NLO QCD", "lf");
   l->Draw();
}

TGraphErrors* GetGraph(Int_t ir, Int_t ifs,Int_t icut, Int_t ipdf)
{
   const char *cuts[5] =
      {"0.0 <= |eta| < 0.5",
       "0.5 <= |eta| < 1.0",
       "1.0 <= |eta| < 1.5",
       "1.5 <= |eta| < 2.0",
       "2.0 <= |eta| < 3.0"};

   const double mur[] = {1.0,0.25,4.0};
   const double muf[] = {1.0,0.25,4.0};

   TFile::SetCacheFileDir(".");
   TFile *file = TFile::Open("http://root.cern/files/AtlasGraphs.root", "CACHEREAD");

   char gname[100];
   char tname[100];

   if (ipdf>=0)
      sprintf(tname," E_T (mu_r=%g, mu_f=%g);%s Pdf: %d",mur[ir],muf[ifs],cuts[icut],ipdf);
   else
      sprintf(tname," E_T %s Ms= %d",cuts[icut],-ipdf);

   TGraphErrors* g1 = 0;

   for (int i=1; i<=GMAX; i++) {
      sprintf(gname,"full_%d",i);
      g1 = (TGraphErrors*) file->Get(gname);
      if (!g1) {
         cout << gname << "  not found " << endl;
         return nullptr;
      }

      const char *title = g1->GetTitle();

      if (strcmp(title,tname)==0) break;
      g1 = 0;
   }

   if (!g1) return nullptr;
   return g1;
}

TGraphAsymmErrors* MakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2) {

   TGraphAsymmErrors* g3 = new TGraphAsymmErrors();

   Double_t  x1 = 0., y1 = 0., x2 = 0., y2 = 0., y0 = 0, x3 = 0.;
   Double_t dum;
   for (Int_t i=0; i<g1->GetN(); i++) {
      g0->GetPoint(i, x1, y0);
      g1->GetPoint(i, x1, y1);
      g2->GetPoint(i, x1, y2);

      if (i==g1->GetN()-1) x2=x1;
      else                 g2->GetPoint(i+1,x2,dum);

      if (i==0)            x3=x1;
      else                 g2->GetPoint(i-1,x3,dum);

      Double_t tmp = y2;
      if (y1 < y2) {
         y2 = y1;
         y1 = tmp;
      }
      g3->SetPoint(i,x1,y0);

      Double_t binwl = (x1-x3)/2.;
      Double_t binwh = (x2-x1)/2.;
      if (binwl == 0.)  binwl = binwh;
      if (binwh == 0.)  binwh = binwl;
      g3->SetPointError(i, binwl, binwh, y0-y2, y1-y0);

   }
   return g3;
}

void AddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2) {

   Double_t x1=0., y1=0., y2=0., y0=0;

   if (g1->GetN()!=g2->GetN())
      cout << " graphs don't have the same number of elements " << endl;

   Double_t* EYhigh = g2-> GetEYhigh();
   Double_t* EYlow  = g2-> GetEYlow();

   for (Int_t i=0; i<g1->GetN(); i++) {
      g1->GetPoint(i, x1, y1);
      g2->GetPoint(i, x1, y2);

      if ( y1==0 || y2==0 )
         cerr << "check these points very carefully : AddtoBand() : point " << i << endl;

      Double_t eyh=0., eyl=0.;

      y0 = y1-y2;
      if (y0 != 0) {
         if (y0 > 0) {
            eyh = EYhigh[i];
            eyh = sqrt(eyh*eyh+y0*y0);
            g2->SetPointEYhigh(i, eyh);
         } else {
            eyl = EYlow[i];
            eyl = sqrt(eyl*eyl+y0*y0);
            g2->SetPointEYlow (i, eyl);
         }
      }
   }
}
