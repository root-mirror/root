/// \file
/// \ingroup tutorial_dataframe
/// \notebook -js
/// Read an sqlite3 databases with RDataFrame and plot statistics on ROOT downloads.
///
/// Plot the downloads of different ROOT versions reading a remote sqlite3 file with RSqliteDS.
/// Then a TH1F histogram is created and filled
/// using a lambda expression which receives the recorded
/// values in the "version" column of the sqlite3 database.
/// The histogram shows the usage of the ROOT development version.
///
/// \macro_code
/// \macro_image
///
/// \date August 2018
/// \author Alexandra-Maria Dobrescu

void df030_SQliteVersionsOfROOT() {

   auto rdf = ROOT::RDF::MakeSqliteDataFrame("http://root.cern/files/root_download_stats.sqlite", "SELECT Version FROM accesslog;");

   TH1F hVersionOfRoot("hVersionOfRoot", "Development Versions of ROOT", 8, 0, -1);

   auto fillVersionHisto = [&hVersionOfRoot] (const std::string &version) {
      TString copyVersion = version;
      TString shortVersion(copyVersion(0,4));
      hVersionOfRoot.Fill(shortVersion, 1);
   };

   rdf.Foreach( fillVersionHisto, { "Version" } );

   auto VersionOfRootHistogram = new TCanvas();

   gStyle->SetOptStat(0);
   hVersionOfRoot.GetXaxis()->LabelsOption("a");
   hVersionOfRoot.LabelsDeflate("X");
   hVersionOfRoot.DrawClone("");
}
