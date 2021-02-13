/// \file
/// \ingroup tutorial_dataframe
/// \notebook -js
/// Use RDataFrame to display data about ROOT downloads.
///
/// In order to display the Platform Distribution of ROOT, we choose to create two TH1F
/// histograms: one that includes all types of platforms, other filtering and classifying them.
/// This procedure is using a lambda expression taking as parameter the values
/// stored in the "Platform" column from the database. At the end, the histograms are filled
/// with their specific demand regarding the platform's type.
///
/// \macro_code
/// \macro_image
///
/// \date August 2018
/// \author Alexandra-Maria Dobrescu

void df029_SQlitePlatformDistribution() {

   auto rdf = ROOT::RDF::MakeSqliteDataFrame( "http://root.cern/files/root_download_stats.sqlite", "SELECT * FROM accesslog;" );

   TH1F hRootPlatform("hrootPlatform", "Platform Distribution", 7, 0, -1);
   TH1F hShortRootPlatform("hShortRootPlatform", "Short Platform Distribution", 7, 0, -1);

   auto fillRootPlatform = [&hRootPlatform, &hShortRootPlatform] ( const std::string &platform ) {
      TString Platform = platform;
      TString Platform_0(Platform(0,5));
      TString Platform_1(Platform(0,6));
      TString Platform_2(Platform(0,8));

      if ( Platform.Contains("win32") ){
        hShortRootPlatform.Fill(Platform_0,1);
      } else if ( Platform.Contains("Linux") ){
        hShortRootPlatform.Fill(Platform_0,1);
      } else if ( Platform.Contains("source") ){
        hShortRootPlatform.Fill(Platform_1,1);
      } else if ( Platform.Contains("macosx64") ){
        hShortRootPlatform.Fill(Platform_2,1);
      } else if ( Platform.Contains("IRIX64") ){
        hShortRootPlatform.Fill(Platform_1,1);
      }

      hRootPlatform.Fill(Platform,1);
   };

   rdf.Foreach( fillRootPlatform, { "Platform" } );

   auto PlatformDistributionHistogram = new TCanvas();

   hRootPlatform.GetXaxis()->LabelsOption("a");
   hRootPlatform.LabelsDeflate("X");
   hRootPlatform.DrawClone();

   auto shortPlatformDistributionHistogram = new TCanvas();

   hShortRootPlatform.GetXaxis()->LabelsOption("a");
   hShortRootPlatform.LabelsDeflate("X");
   hShortRootPlatform.DrawClone();
}
