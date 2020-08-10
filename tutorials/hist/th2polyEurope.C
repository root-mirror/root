/// \file
/// \ingroup tutorial_hist
/// \notebook -js
/// This tutorial illustrates how to create an histogram with polygonal
/// bins (TH2Poly), fill it and draw it. The initial data are stored
/// in TMultiGraphs. They represent the european countries.
/// The histogram filling is done according to a Mercator projection,
/// therefore the bin contains should be proportional to the real surface
/// of the countries.
///
/// The initial data have been downloaded from: http://www.maproom.psu.edu/dcw/
/// This database was developed in 1991/1992 and national boundaries reflect
/// political reality as of that time.
///
/// The script is shooting npoints (script argument) randomly over the Europe area.
/// The number of points inside the countries should be proportional to the country surface
/// The estimated surface is compared to the surfaces taken from wikipedia.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Olivier Couet

void th2polyEurope(Int_t npoints=500000)
{
   Int_t i,j;
   Double_t lon1 = -25;
   Double_t lon2 =  35;
   Double_t lat1 =  34;
   Double_t lat2 =  72;
   Double_t R = (lat2-lat1)/(lon2-lon1);
   Int_t W    = 800;
   Int_t H    = (Int_t)(R*800);
   gStyle->SetStatX(0.28);
   gStyle->SetStatY(0.45);
   gStyle->SetStatW(0.15);

   // Canvas used to draw TH2Poly (the map)
   TCanvas *ce = new TCanvas("ce", "ce",0,0,W,H);
   ce->ToggleEventStatus();
   ce->SetGridx();
   ce->SetGridy();

   // Real surfaces taken from Wikipedia.
   const Int_t nx = 36;
   // see http://en.wikipedia.org/wiki/Area_and_population_of_European_countries
   const char *countries[nx] = {
      "france",     "spain",  "sweden",  "germany",       "finland",
      "norway",     "poland", "italy",   "yugoslavia",    "united_kingdom",
      "romania",    "belarus","greece",  "czechoslovakia","bulgaria",
      "iceland",    "hungary","portugal","austria",       "ireland",
      "lithuania",  "latvia", "estonia", "denmark",       "netherlands",
      "switzerland","moldova","belgium", "albania",       "cyprus",
      "luxembourg", "andorra","malta",   "liechtenstein", "san_marino", "monaco" };
   Float_t surfaces[nx] = {
      547030,        505580,   449964,      357021,        338145,
      324220,        312685,   301230,      255438,        244820,
      237500,        207600,   131940,      127711,        110910,
      103000,         93030,    89242,       83870,         70280,
      65200,         64589,    45226,       43094,         41526,
      41290,         33843,    30528,       28748,          9250,
      2586,           468,      316,         160,            61, 2};

   TH1F *h = new TH1F("h","Countries surfaces (in km^{2})",3,0,3);
   for (i=0; i<nx; i++) h->Fill(countries[i], surfaces[i]);
   h->LabelsDeflate();

   TFile::SetCacheFileDir(".");
   TFile *f;
   f = TFile::Open("http://root.cern/files/europe.root","cacheread");

   if (!f) {
      printf("Cannot access europe.root. Is internet working ?\n");
      return;
   }

   TH2Poly *p = new TH2Poly(
             "Europe",
             "Europe (bin contents are normalized to the surfaces in km^{2})",
             lon1,lon2,lat1,lat2);
   p->GetXaxis()->SetNdivisions(520);
   p->GetXaxis()->SetTitle("longitude");
   p->GetYaxis()->SetTitle("latitude");

   p->SetContour(100);

   TMultiGraph *mg;
   TKey *key;
   TIter nextkey(gDirectory->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();
      if (obj->InheritsFrom("TMultiGraph")) {
         mg = (TMultiGraph*)obj;
         p->AddBin(mg);
      }
   }

   TRandom r;
   Double_t longitude, latitude;
   Double_t x, y, pi4 = TMath::Pi()/4, alpha = TMath::Pi()/360;

   gBenchmark->Start("Partitioning");
   p->ChangePartition(100, 100);
   gBenchmark->Show("Partitioning");

   // Fill TH2Poly according to a Mercator projection.
   gBenchmark->Start("Filling");
   for (i=0; i<npoints; i++) {
      longitude = r.Uniform(lon1,lon2);
      latitude  = r.Uniform(lat1,lat2);
      x         = longitude;
      y         = 38*TMath::Log(TMath::Tan(pi4+alpha*latitude));
      p->Fill(x,y);
   }
   gBenchmark->Show("Filling");

   Int_t nbins = p->GetNumberOfBins();
   Double_t maximum = p->GetMaximum();


   // h2 contains the surfaces computed from TH2Poly.
   TH1F *h2 = (TH1F *)h->Clone("h2");
   h2->Reset();
   for (j=0; j<nx; j++) {
      for (i=0; i<nbins; i++) {
         if (strstr(countries[j],p->GetBinName(i+1))) {
            h2->Fill(countries[j],p->GetBinContent(i+1));
            h2->SetBinError(j, p->GetBinError(i+1));
         }
      }
   }

   // Normalize the TH2Poly bin contents to the real surfaces.
   Double_t scale = surfaces[0]/maximum;
   for (i=0; i<nbins; i++) p->SetBinContent(i+1, scale*p->GetBinContent(i+1));

   gStyle->SetOptStat(1111);
   p->Draw("COLZ");

   TCanvas *c1 = new TCanvas("c1", "c1",W+10,0,W-20,H);
   c1->SetRightMargin(0.047);

   scale = h->GetMaximum()/h2->GetMaximum();

   h->SetStats(0);
   h->SetLineColor(kRed-3);
   h->SetLineWidth(2);
   h->SetMarkerStyle(20);
   h->SetMarkerColor(kBlue);
   h->SetMarkerSize(0.8);
   h->Draw("LP");
   h->GetXaxis()->SetLabelFont(42);
   h->GetXaxis()->SetLabelSize(0.03);
   h->GetYaxis()->SetLabelFont(42);

   h2->Scale(scale);
   Double_t scale2=TMath::Sqrt(scale);
   for (i=0; i<nx; i++) h2->SetBinError(i+1, scale2*h2->GetBinError(i+1));
   h2->Draw("E SAME");
   h2->SetMarkerStyle(20);
   h2->SetMarkerSize(0.8);

   TLegend *leg = new TLegend(0.5,0.67,0.92,0.8,NULL,"NDC");
   leg->SetTextFont(42);
   leg->SetTextSize(0.025);
   leg->AddEntry(h,"Real countries surfaces from Wikipedia (in km^{2})","lp");
   leg->AddEntry(h2,"Countries surfaces from TH2Poly (with errors)","lp");
   leg->Draw();
   leg->Draw();

   Double_t wikiSum = h->Integral();
   Double_t polySum = h2->Integral();
   Double_t error = TMath::Abs(wikiSum-polySum)/wikiSum;
   printf("THPoly Europe surface estimation error wrt wikipedia = %f per cent when using %d points\n",100*error,npoints);
}
