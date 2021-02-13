#include "TMVA/tmvaglob.h"
#include "TMVA/Config.h"

using std::cout;
using std::endl;

// set the style
void TMVA::TMVAGlob::SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all )
{
   //signal
   // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
   // convince yourself with gROOT->GetListOfColors()->Print()
   Int_t FillColor__S = getSignalFill();
   Int_t FillStyle__S = 1001;
   Int_t LineColor__S = getSignalLine();
   Int_t LineWidth__S = 2;

   // background
   //Int_t icolor = gConfig().fVariablePlotting.fUsePaperStyle ? 2 + 100 : 2;
   Int_t FillColor__B = getBackgroundFill();
   Int_t FillStyle__B = 3554;
   Int_t LineColor__B = getBackgroundLine();
   Int_t LineWidth__B = 2;

   if (sig != NULL) {
      sig->SetLineColor( LineColor__S );
      sig->SetLineWidth( LineWidth__S );
      sig->SetFillStyle( FillStyle__S );
      sig->SetFillColor( FillColor__S );
   }

   if (bkg != NULL) {
      bkg->SetLineColor( LineColor__B );
      bkg->SetLineWidth( LineWidth__B );
      bkg->SetFillStyle( FillStyle__B );
      bkg->SetFillColor( FillColor__B );
   }

   if (all != NULL) {
      all->SetLineColor( LineColor__S );
      all->SetLineWidth( LineWidth__S );
      all->SetFillStyle( FillStyle__S );
      all->SetFillColor( FillColor__S );
   }
}

void TMVA::TMVAGlob::SetMultiClassStyle( TObjArray* hists )
{
   //signal
   // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
   // convince yourself with gROOT->GetListOfColors()->Print()
   //Int_t FillColor__S = getSignalFill();
   //Int_t FillStyle__S = 1001;
   //Int_t LineColor__S = getSignalLine();
   //Int_t LineWidth__S = 2;

   // background
   //Int_t icolor = gConfig().fVariablePlotting.fUsePaperStyle ? 2 + 100 : 2;
   //Int_t FillColor__B = getBackgroundFill();
   //Int_t FillStyle__B = 3554;
   //Int_t LineColor__B = getBackgroundLine();
   //Int_t LineWidth__B = 2;

   Int_t FillColors[10] = {38,2,3,6,7,8,9,11};
   Int_t LineColors[10] = {4,2,3,6,7,8,9,11};
   Int_t FillStyles[5] = {1001,3554,3003,3545,0};

   for(Int_t i=0; i<hists->GetEntriesFast(); ++i){
      ((TH1*)(*hists)[i])->SetFillColor(FillColors[i%10]);
      ((TH1*)(*hists)[i])->SetFillStyle(FillStyles[i%5]);
      ((TH1*)(*hists)[i])->SetLineColor(LineColors[i%10]);
      ((TH1*)(*hists)[i])->SetLineWidth(2);
   }
}

// set frame styles
void TMVA::TMVAGlob::SetFrameStyle( TH1* frame, Float_t scale )
{
   frame->SetLabelOffset( 0.012, "X" );// label offset on x axis
   frame->SetLabelOffset( 0.012, "Y" );// label offset on x axis
   frame->GetXaxis()->SetTitleOffset( 1.25 );
   frame->GetYaxis()->SetTitleOffset( 1.22 );
   frame->GetXaxis()->SetTitleSize( 0.045*scale );
   frame->GetYaxis()->SetTitleSize( 0.045*scale );
   Float_t labelSize = 0.04*scale;
   frame->GetXaxis()->SetLabelSize( labelSize );
   frame->GetYaxis()->SetLabelSize( labelSize );

   // global style settings
   gPad->SetTicks();
   gPad->SetLeftMargin  ( 0.108*scale );
   gPad->SetRightMargin ( 0.050*scale );
   gPad->SetBottomMargin( 0.120*scale  );
}

void TMVA::TMVAGlob::SetTMVAStyle() {

   TStyle *TMVAStyle = gROOT->GetStyle("TMVA");
   if(TMVAStyle!=0) {
      gROOT->SetStyle("TMVA");
      return;
   }

   TMVAStyle = new TStyle(*gROOT->GetStyle("Plain")); // our style is based on Plain
   TMVAStyle->SetName("TMVA");
   TMVAStyle->SetTitle("TMVA style based on \"Plain\" with modifications defined in tmvaglob.C");
   gROOT->GetListOfStyles()->Add(TMVAStyle);
   gROOT->SetStyle("TMVA");

   TMVAStyle->SetLineStyleString( 5, "[52 12]" );
   TMVAStyle->SetLineStyleString( 6, "[22 12]" );
   TMVAStyle->SetLineStyleString( 7, "[22 10 7 10]" );

   // the pretty color palette of old
   TMVAStyle->SetPalette((gConfig().fVariablePlotting.fUsePaperStyle ? 18 : 1),0);

   // use plain black on white colors
   TMVAStyle->SetFrameBorderMode(0);
   TMVAStyle->SetCanvasBorderMode(0);
   TMVAStyle->SetPadBorderMode(0);
   TMVAStyle->SetPadColor(0);
   TMVAStyle->SetFillStyle(0);

   TMVAStyle->SetLegendBorderSize(0);

   // title properties
   // TMVAStyle->SetTitleW(.4);
   // TMVAStyle->SetTitleH(.10);
   // MVAStyle->SetTitleX(.5);
   // TMVAStyle->SetTitleY(.9);
   TMVAStyle->SetTitleFillColor( getTitleBox() );
   TMVAStyle->SetTitleTextColor( getTitleText() );
   TMVAStyle->SetTitleBorderSize( 1 );
   TMVAStyle->SetLineColor( getTitleBorder() );
   if (!gConfig().fVariablePlotting.fUsePaperStyle) {
      TMVAStyle->SetFrameFillColor( getFrameFill() );
      TMVAStyle->SetCanvasColor( getCanvas() );
   }

   // set the paper & margin sizes
   TMVAStyle->SetPaperSize(20,26);
   TMVAStyle->SetPadTopMargin(0.10);
   TMVAStyle->SetPadRightMargin(0.05);
   TMVAStyle->SetPadBottomMargin(0.11);
   TMVAStyle->SetPadLeftMargin(0.12);

   // use bold lines and markers
   TMVAStyle->SetMarkerStyle(21);
   TMVAStyle->SetMarkerSize(0.3);
   TMVAStyle->SetHistLineWidth(2);
   TMVAStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

   // do not display any of the standard histogram decorations
   TMVAStyle->SetOptTitle(1);
   TMVAStyle->SetTitleH(0.052);

   TMVAStyle->SetOptStat(0);
   TMVAStyle->SetOptFit(0);

   // put tick marks on top and RHS of plots
   TMVAStyle->SetPadTickX(1);
   TMVAStyle->SetPadTickY(1);

}

void TMVA::TMVAGlob::DestroyCanvases()
{

   TList* loc = (TList*)gROOT->GetListOfCanvases();
   TListIter itc(loc);
   TObject *o(0);
   while ((o = itc())) delete o;
}

// set style and remove existing canvas'
void TMVA::TMVAGlob::Initialize( Bool_t useTMVAStyle )
{
   // destroy canvas'
   DestroyCanvases();

   // set style
   if (!useTMVAStyle) {
      gROOT->SetStyle("Plain");
      gStyle->SetOptStat(0);
      return;
   }

   SetTMVAStyle();
}

// checks if file with name "fin" is already open, and if not opens one
TFile* TMVA::TMVAGlob::OpenFile( const TString& fin )
{
   TFile* file = gDirectory->GetFile();
   if (file==0 || fin != file->GetName()) {
      if (file != 0) {
         gROOT->cd();
         file->Close();
      }
      cout << "--- Opening root file " << fin << " in read mode" << endl;
      file = TFile::Open( fin, "READ" );
   }
   else {
      file = gDirectory->GetFile();
   }

   file->cd();
   return file;
}

// used to create output file for canvas
void TMVA::TMVAGlob::imgconv( TCanvas* c, const TString & fname )
{
   // return;
   if (NULL == c) {
      cout << "*** Error in TMVAGlob::imgconv: canvas is NULL" << endl;
   }
   else {
      // create directory if not existing
      TString f = fname;
      TString dir = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
      // directory does not exists , try to make it
      if (gSystem->AccessPathName(dir)) {
         if (gSystem->mkdir(dir, kTRUE) != 0) {
            Error("imgconv","Error creating plot directory: %s",dir.Data());
         }
      }

      TString pngName = fname + ".png";
      TString gifName = fname + ".gif";
      TString epsName = fname + ".eps";
      TString pdfName = fname + ".pdf";
      c->cd();

      // create eps (other option: c->Print( epsName ))
      if (gConfig().fVariablePlotting.fUsePaperStyle) {
         c->Print(epsName);
      }
      else {
         if (gConfig().fVariablePlotting.fPlotFormat == Config::VariablePlotting::kGIF)
            c->Print(gifName);
         else if (gConfig().fVariablePlotting.fPlotFormat == Config::VariablePlotting::kPDF)
            c->Print(pdfName);
         else if (gConfig().fVariablePlotting.fPlotFormat == Config::VariablePlotting::kEPS)
            c->Print(epsName);
         else
            c->Print(pngName);
      }
   }
}

TImage * TMVA::TMVAGlob::findImage(const char * imageName)
{
   // looks for the image in tutorialpath
   //TString tutorialPath = "$ROOTSYS/tutorials/tmva"; // look for the image in here
   TString tutorialPath = getenv ("ROOTSYS");
   tutorialPath+="/tutorials/tmva";
   TImage *img(0);
   TString fullName = Form("%s/%s", tutorialPath.Data(), imageName);
   Bool_t fileFound = ! gSystem->AccessPathName(fullName);

   if(fileFound) {
      img = TImage::Open(fullName);
   } else {
      cout << "+++ Could not open image:  " << fullName << endl;
   }
   return img;
}

void TMVA::TMVAGlob::plot_logo( Float_t v_scale, Float_t skew )
{

   TImage *img = findImage("tmva_logo.gif");
   if (!img) {
      cout << "+++ Could not open image tmva_logo.gif" << endl;
      return;
   }

   img->SetConstRatio(kFALSE);
   UInt_t h_ = img->GetHeight();
   UInt_t w_ = img->GetWidth();

   Float_t r = w_/h_;
   gPad->Update();
   Float_t rpad = Double_t(gPad->VtoAbsPixel(0) - gPad->VtoAbsPixel(1))/(gPad->UtoAbsPixel(1) - gPad->UtoAbsPixel(0));
   r *= rpad;

   Float_t d = 0.055;
   // absolute coordinates
   Float_t x1R = 1 - gStyle->GetPadRightMargin();
   Float_t y1B = 1 - gStyle->GetPadTopMargin()+.01; // we like the logo to sit a bit above the histo

   Float_t x1L = x1R - d*r/skew;
   Float_t y1T = y1B + d*v_scale*skew;
   if (y1T>0.99) y1T = 0.99;

   TPad *p1 = new TPad("imgpad", "imgpad", x1L, y1B, x1R, y1T );
   p1->SetRightMargin(0);
   p1->SetBottomMargin(0);
   p1->SetLeftMargin(0);
   p1->SetTopMargin(0);
   p1->Draw();

   Int_t xSizeInPixel = p1->UtoAbsPixel(1) - p1->UtoAbsPixel(0);
   Int_t ySizeInPixel = p1->VtoAbsPixel(0) - p1->VtoAbsPixel(1);
   if (xSizeInPixel<=25 || ySizeInPixel<=25) {
      delete p1;
      return; // ROOT doesn't draw smaller than this
   }

   p1->cd();
   img->Draw();
}

void TMVA::TMVAGlob::NormalizeHist( TH1* h )
{
   if (h==0) return;
   if (h->GetSumw2N() == 0) h->Sumw2();
   if(h->GetSumOfWeights()!=0) {
      Float_t dx = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX();
      h->Scale( 1.0/h->GetSumOfWeights()/dx );
   }
}
void TMVA::TMVAGlob::NormalizeHists( TH1* sig, TH1* bkg )
{
   if (sig->GetSumw2N() == 0) sig->Sumw2();
   if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();

   if(sig->GetSumOfWeights()!=0) {
      Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
      sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
   }
   if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
      Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
      bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
   }
}

// the following are tools to help handling different methods and titles


void TMVA::TMVAGlob::GetMethodName( TString & name, TKey * mkey ) {
   if (mkey==0) return;
   name = mkey->GetName();
   name.ReplaceAll("Method_","");
}

void TMVA::TMVAGlob::GetMethodTitle( TString & name, TKey * ikey ) {
   if (ikey==0) return;
   name = ikey->GetName();
}

void TMVA::TMVAGlob::GetMethodName( TString & name, TDirectory * mdir ) {
   if (mdir==0) return;
   name = mdir->GetName();
   name.ReplaceAll("Method_","");
}

void TMVA::TMVAGlob::GetMethodTitle( TString & name, TDirectory * idir ) {
   if (idir==0) return;
   name = idir->GetName();
}

TKey *TMVA::TMVAGlob::NextKey( TIter & keyIter, TString className) {
   TKey *key=(TKey *)keyIter.Next();
   TKey *rkey=0;
   Bool_t loop=(key!=0);
   //
   while (loop) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (cl->InheritsFrom(className.Data())) {
         loop = kFALSE;
         rkey = key;
      } else {
         key = (TKey *)keyIter.Next();
         if (key==0) loop = kFALSE;
      }
   }
   return rkey;
}

UInt_t TMVA::TMVAGlob::GetListOfKeys( TList& keys, TString inherits, TDirectory *dir )
{
   // get a list of keys with a given inheritance
   // the list contains TKey objects
   if (dir==0) dir = gDirectory;
   TIter mnext(dir->GetListOfKeys());
   TKey *mkey;
   keys.Clear();
   keys.SetOwner(kFALSE);
   UInt_t ni=0;
   while ((mkey = (TKey*)mnext())) {
      // make sure, that we only look at TDirectory with name Method_<xxx>
      TClass *cl = gROOT->GetClass(mkey->GetClassName());
      if (cl->InheritsFrom(inherits)) {
         keys.Add(mkey);
         ni++;
      }
   }
   return ni;
}

Int_t TMVA::TMVAGlob::GetNumberOfTargets( TDirectory *dir )
{
   if (!dir) {
      cout << "tmvaglob::GetNumberOfTargets is called with *dir==NULL :( " << endl;
      return 0;
   }
   TIter next(dir->GetListOfKeys());
   TKey* key    = 0;
   Int_t noTrgts = 0;

   while ((key = (TKey*)next())) {
      if (key->GetCycle() != 1) continue;
      if (TString(key->GetName()).Contains("__Regression_target")) noTrgts++;
   }
   return noTrgts;
}

Int_t TMVA::TMVAGlob::GetNumberOfInputVariables( TDirectory *dir )
{
   TIter next(dir->GetListOfKeys());
   TKey* key    = 0;
   Int_t noVars = 0;

   while ((key = (TKey*)next())) {
      if (key->GetCycle() != 1) continue;

      // count number of variables (signal is sufficient), exclude target(s)
      if (TString(key->GetName()).Contains("__Signal") || (TString(key->GetName()).Contains("__Regression") && !(TString(key->GetName()).Contains("__Regression_target")))) noVars++;
   }

   return noVars;
}

std::vector<TString> TMVA::TMVAGlob::GetInputVariableNames(TDirectory *dir)
{
   TIter next(dir->GetListOfKeys());
   TKey* key = 0;
   std::vector<TString> names;

   while ((key = (TKey*)next())) {
      if (key->GetCycle() != 1) continue;
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TString name(key->GetName());
      Int_t pos = name.Index("__");
      if (pos != -1) {
         name.Remove(pos);
      }
      Bool_t hasname = false;
      std::vector<TString>::const_iterator iter = names.begin();
      while(iter != names.end()){
         if(name.CompareTo(*iter)==0)
            hasname=true;
         ++iter;
      }
      if(!hasname)
         names.push_back(name);
   }
   return names;
}

Int_t TMVA::TMVAGlob::GetNumberOfInputVariablesMultiClass( TDirectory *dir ){
   std::vector<TString> names(GetInputVariableNames(dir));
   return names.end() - names.begin();
}

std::vector<TString> TMVA::TMVAGlob::GetClassNames(TDirectory *dir )
{

   TIter next(dir->GetListOfKeys());
   TKey* key = 0;
   //set<std::string> varnames;
   std::vector<TString> names;

   while ((key = (TKey*)next())) {
      if (key->GetCycle() != 1) continue;
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TString name(key->GetName());
      Int_t pos = name.Index("__");
      if (pos == -1)
         continue;
      name.ReplaceAll("_Deco","");
      name.ReplaceAll("_Gauss","");
      name.ReplaceAll("_PCA","");
      name.ReplaceAll("_Id","");
      name.ReplaceAll("_vs_","");
      name.Remove(0, pos + 2);

      /*Int_t pos = name.First("__");
        name.Remove(0,pos+2);
        char c = '_';
        pos = name.Last(c);
        name.Remove(pos);
        if(name.Contains("Gauss")){
        pos = name.Last(c);
        name.Remove(pos);
        }
        pos = name.Last(c);
        if(pos!=-1)
        name.Remove(0,pos+1);
      */
      Bool_t hasname = false;
      std::vector<TString>::const_iterator iter = names.begin();
      while(iter != names.end()){
         if(name.CompareTo(*iter)==0)
            hasname=true;
         ++iter;
      }
      if(!hasname)
         names.push_back(name);
   }
   return names;
}


TKey* TMVA::TMVAGlob::FindMethod( TString name, TDirectory *dir )
{
   // find the key for a method
   if (dir==0) dir = gDirectory;
   TIter mnext(dir->GetListOfKeys());
   TKey *mkey;
   TKey *retkey=0;
   Bool_t loop=kTRUE;
   while (loop) {
      mkey = (TKey*)mnext();
      if (mkey==0) {
         loop = kFALSE;
      }
      else {
         TString clname = mkey->GetClassName();
         TClass *cl = gROOT->GetClass(clname);
         if (cl->InheritsFrom("TDirectory")) {
            TString mname = mkey->GetName(); // method name
            TString tname = "Method_"+name;  // target name
            if (mname==tname) { // target found!
               loop = kFALSE;
               retkey = mkey;
            }
         }
      }
   }
   return retkey;
}

Bool_t TMVA::TMVAGlob::ExistMethodName( TString name, TDirectory *dir )
{
   // find the key for a method
   if (dir==0) dir = gDirectory;
   TIter mnext(dir->GetListOfKeys());
   TKey *mkey;
   Bool_t loop=kTRUE;
   while (loop) {
      mkey = (TKey*)mnext();
      if (mkey==0) {
         loop = kFALSE;
      }
      else {
         TString clname  = mkey->GetClassName();
         TString keyname = mkey->GetName();
         TClass *cl = gROOT->GetClass(clname);
         if (keyname.Contains("Method") && cl->InheritsFrom("TDirectory")) {

            TDirectory* d_ = (TDirectory*)dir->Get( keyname );
            if (!d_) {
               cout << "HUUUGE TROUBLES IN TMVAGlob::ExistMethodName() --> contact authors" << endl;
               return kFALSE;
            }

            TIter mnext_(d_->GetListOfKeys());
            TKey *mkey_;
            while ((mkey_ = (TKey*)mnext_())) {
               TString clname_ = mkey_->GetClassName();
               TClass *cl_ = gROOT->GetClass(clname_);
               if (cl_->InheritsFrom("TDirectory")) {
                  TString mname = mkey_->GetName(); // method name
                  if (mname==name) { // target found!
                     return kTRUE;
                  }
               }
            }
         }
      }
   }
   return kFALSE;
}

UInt_t TMVA::TMVAGlob::GetListOfMethods( TList & methods, TDirectory *dir )
{
   // get a list of methods
   // the list contains TKey objects
   if (dir==0) dir = gDirectory;
   TIter mnext(dir->GetListOfKeys());
   TKey *mkey;
   methods.Clear();
   methods.SetOwner(kFALSE);
   UInt_t ni=0;
   while ((mkey = (TKey*)mnext())) {
      // make sure, that we only look at TDirectory with name Method_<xxx>
      TString name = mkey->GetClassName();
      TClass *cl = gROOT->GetClass(name);
      if (cl->InheritsFrom("TDirectory")) {
         if (TString(mkey->GetName()).BeginsWith("Method_")) {
            methods.Add(mkey);
            ni++;
         }
      }
   }
   cout << "--- Found " << ni << " classifier types" << endl;
   return ni;
}

UInt_t TMVA::TMVAGlob::GetListOfJobs( TFile* file, TList& jobdirs)
{
   // get a list of all jobs in all method directories
   // based on ideas by Peter and Joerg found in macro deviations.C
   TIter next(file->GetListOfKeys());
   TKey *key(0);
   while ((key = (TKey*)next())) {

      if (TString(key->GetName()).BeginsWith("Method_")) {
         if (gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) {

            TDirectory* mDir = (TDirectory*)key->ReadObj();

            TIter keyIt(mDir->GetListOfKeys());
            TKey *jobkey;
            while ((jobkey = (TKey*)keyIt())) {
               if (!gROOT->GetClass(jobkey->GetClassName())->InheritsFrom("TDirectory")) continue;

               TDirectory *jobDir = (TDirectory *)jobkey->ReadObj();
               cout << "jobdir name  " << jobDir->GetName() << endl;
               jobdirs.Add(jobDir);
            }
         }
      }
   }
   return jobdirs.GetSize();
}

UInt_t TMVA::TMVAGlob::GetListOfTitles( TDirectory *rfdir, TList & titles )
{
   // get a list of titles (i.e TDirectory) given a method dir
   UInt_t ni=0;
   if (rfdir==0) return 0;
   TList *keys = rfdir->GetListOfKeys();
   if (keys==0) {
      cout << "+++ Directory '" << rfdir->GetName() << "' contains no keys" << endl;
      return 0;
   }
   //
   TIter rfnext(rfdir->GetListOfKeys());
   TKey *rfkey;
   titles.Clear();
   titles.SetOwner(kFALSE);
   while ((rfkey = (TKey*)rfnext())) {
      // make sure, that we only look at histograms
      TClass *cl = gROOT->GetClass(rfkey->GetClassName());
      if (cl->InheritsFrom("TDirectory")) {
         titles.Add(rfkey);
         ni++;
      }
   }
   cout << "--- Found " << ni << " instance(s) of the method " << rfdir->GetName() << endl;
   return ni;
}

UInt_t TMVA::TMVAGlob::GetListOfTitles( TString & methodName, TList & titles, TDirectory *dir )
{
   // get the list of all titles for a given method
   // if the input dir is 0, gDirectory is used
   // returns a list of keys
   UInt_t ni=0;
   if (dir==0) dir = gDirectory;
   TDirectory* rfdir = (TDirectory*)dir->Get( methodName );
   if (rfdir==0) {
      cout << "+++ Could not locate directory '" << methodName << endl;
      return 0;
   }

   return GetListOfTitles( rfdir, titles );

   TList *keys = rfdir->GetListOfKeys();
   if (keys==0) {
      cout << "+++ Directory '" << methodName << "' contains no keys" << endl;
      return 0;
   }
   //
   TIter rfnext(rfdir->GetListOfKeys());
   TKey *rfkey;
   titles.Clear();
   titles.SetOwner(kFALSE);
   while ((rfkey = (TKey*)rfnext())) {
      // make sure, that we only look at histograms
      TClass *cl = gROOT->GetClass(rfkey->GetClassName());
      if (cl->InheritsFrom("TDirectory")) {
         titles.Add(rfkey);
         ni++;
      }
   }
   cout << "--- Found " << ni << " instance(s) of the method " << methodName << endl;
   return ni;
}

TDirectory *TMVA::TMVAGlob::GetInputVariablesDir( TMVAGlob::TypeOfPlot type, TDirectory *dir )
{
   // get the InputVariables directory
   const TString directories[TMVAGlob::kNumOfMethods] = { "InputVariables_Id",
                                                          "InputVariables_Deco",
                                                          "InputVariables_PCA",
                                                          "InputVariables_Gauss_Deco" };
   if (dir==0) dir = gDirectory;

   // get top dir containing all hists of the variables
   dir = (TDirectory*)gDirectory->Get( directories[type] );
   if (dir==0) {
      cout << "+++ Could not locate input variable directory '" << directories[type] << endl;
      return 0;
   }
   return dir;
}

TDirectory *TMVA::TMVAGlob::GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir )
{
   // get the CorrelationPlots directory
   if (dir==0) dir = GetInputVariablesDir( type, 0 );
   if (dir==0) return 0;
   //
   TDirectory* corrdir = (TDirectory*)dir->Get( "CorrelationPlots" );
   if (corrdir==0) {
      cout << "+++ Could not find CorrelationPlots directory 'CorrelationPlots'" << endl;
      return 0;
   }
   return corrdir;
}
