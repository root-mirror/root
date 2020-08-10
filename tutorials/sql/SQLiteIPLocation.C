/// \file
/// \ingroup tutorial_sql
/// \notebook -js
///
/// This tutorial demonstrates how TSQLServer can be used to create a
/// connection with a SQlite3 database. It accesses the Sqlite data base.
/// Download from https://root.cern/download/root_download_stats.sqlite
/// The world map is hold by a TH2Poly histogram which, after filling, will show
/// the world wide dispersion of ROOT's users.
/// To histogram filling, is done having as input parameters
/// the two columns of the database: "IPLongitude' - for the longitude, and the
/// "IPLatitude" - for the latitude.
/// The data related to the latitude and the longitude has been provided from the
/// log files storing the users IP.
/// This product includes GeoLite2 data created by MaxMind, available from
/// <a href="http://www.maxmind.com">http://www.maxmind.com</a>.
///
/// \macro_code
///
/// \author Alexandra-Maria Dobrescu 08/2018

#include <TSQLiteServer.h>
#include <TSQLiteResult.h>
#include <TSQLRow.h>
#include <TString.h>
#include <TH2F.h>

void SQLiteIPLocation() {

   TSQLServer *db = TSQLServer::Connect("sqlite://root_download_stats.sqlite", "", "");

   TFile *F = TFile::Open("http://root.cern/files/WM.root");
   TH2Poly *WM;
   WM = (TH2Poly*) F->Get("WM");
   const char *location = "SELECT IPLatitude, IPLongitude FROM accesslog;";
   TSQLResult *locationRes = db->Query(location);

   while (TSQLRow *row = locationRes->Next()) {
      if (!row->GetField(0)[0])
         continue;
      std::string sLatitude(row->GetField(0));
      std::string sLongitude(row->GetField(1));
      float latitude = std::stof(sLatitude);
      float longitude = std::stof(sLongitude);
      WM->Fill(longitude, latitude);

     delete row;
   }

   TCanvas *locationHistogram = new TCanvas();

   locationHistogram->SetLogz(1);
   locationHistogram->ToggleEventStatus();
   WM->Draw("colz");
}
