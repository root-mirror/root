void GetWebHistogram()
{
   // example of script called from an Action on Demand when a TRef object
   // is dereferenced. See Event.h, member fWebHistogram

   const char *URL = "http://root.cern/files/pippa.root";
   printf("GetWebHistogram from URL: %s\n",URL);
   TFile *f= TFile::Open(URL);
   f->cd("DM/CJ");
   TH1 *h6 = (TH1*)gDirectory->Get("h6");
   h6->SetDirectory(0);
   delete f;
   TRef::SetStaticObject(h6);
}
