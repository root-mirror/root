void CreateExampleFile()
{
   gROOT->ProcessLine(".! prepareHistFactory .");
   gROOT->ProcessLine(".! hist2workspace config/example.xml");
}
