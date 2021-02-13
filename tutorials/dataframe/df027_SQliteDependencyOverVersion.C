/// \file
/// \ingroup tutorial_dataframe
/// \notebook -jss
/// Plot the ROOT downloads based on the version reading a remote sqlite3 file.
///
/// This tutorial uses the Reduce method which allows to extract the minimum time
/// stored in the SQlite3 database.
/// The next step is to create a TH1F Histogram, which will be filled with the values stored in
/// two different columns from the database. This procedure is simplified with a lambda
/// expression that takes as parameters the values stored in the "Time" and "Version" columns.
///
/// \macro_code
/// \macro_image
///
/// \date August 2018
/// \author Alexandra-Maria Dobrescu

void df027_SQliteDependencyOverVersion () {

   auto rdfb = ROOT::RDF::MakeSqliteDataFrame("http://root.cern/files/root_download_stats.sqlite", "SELECT * FROM accesslog;" );

   auto minTimeStr = *rdfb.Reduce([](std::string a, std::string b) {return std::min(a, b);}, "Time", std::string("Z"));

   std::cout << "Minimum time is '" << minTimeStr << "'" << std::endl;

   double minTime = TDatime(minTimeStr.c_str()).Convert();
   double now = TDatime().Convert();

   auto rdf = rdfb.Define("datime", [](const std::string &time){return TDatime(time.c_str()).Convert();}, {"Time"});

   auto h614 = rdf.Filter([](const std::string &v){ return 0 == v.find("6.14");}, {"Version"})
                  .Histo1D({"h614", "Download time for version 6.14", 16, minTime, now}, {"datime"});

   auto h616 = rdf.Filter([](const std::string &v){ return 0 == v.find("6.16");}, {"Version"})
                  .Histo1D({"h616", "Download time for version 6.16", 16, minTime, now}, {"datime"});

   auto h618 = rdf.Filter([](const std::string &v){ return 0 == v.find("6.18");}, {"Version"})
                  .Histo1D({"h618", "Download time for version 6.18", 16, minTime, now}, {"datime"});

   // Add here a newer version!

   auto histoList = {h614, h616, h618};
   auto canvases = new std::vector<TCanvas*>(histoList.size());

   gStyle->SetTimeOffset(0);
   gStyle->SetOptStat(0);
   auto histoIdx = 0U;
   for (auto histo : histoList) {
      canvases->at(histoIdx) = new TCanvas();
      histo->GetXaxis()->LabelsOption("v");
      histo->GetXaxis()->SetTimeDisplay(1);
      histo->GetXaxis()->SetLabelSize(0.02);
      histo->GetXaxis()->SetNdivisions(512, kFALSE);
      histo->GetXaxis()->SetTimeFormat("%Y-%m-%d");
      histo->DrawClone();
      histoIdx++;
   }

}
