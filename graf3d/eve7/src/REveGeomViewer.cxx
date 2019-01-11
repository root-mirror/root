// @(#)root/eve7:$Id$
// Author: Sergey Linev, 13.12.2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveGeomViewer.hxx>

#include <ROOT/RWebWindowsManager.hxx>
#include <ROOT/TLogger.hxx>

#include "TSystem.h"
#include "TROOT.h"
#include "THttpServer.h"
#include "TBufferJSON.h"
#include "TGeoManager.h"

ROOT::Experimental::REveGeomViewer::REveGeomViewer(TGeoManager *mgr)
{

   TString evedir = TString::Format("%s/eve7", TROOT::GetEtcDir().Data());

   if (gSystem->ExpandPathName(evedir)) {
      R__WARNING_HERE("webeve") << "Problems resolve path " << evedir << " for HTML sources";
      evedir = ".";
   }

   fWebWindow = ROOT::Experimental::RWebWindowsManager::Instance()->CreateWindow();

   fWebWindow->GetServer()->AddLocation("/evedir/",  evedir.Data());
   fWebWindow->SetDefaultPage(Form("file:%s/geom.html", evedir.Data()));

   // this is call-back, invoked when message received via websocket
   fWebWindow->SetDataCallBack([this](unsigned connid, const std::string &arg) { this->WebWindowCallback(connid, arg); });
   fWebWindow->SetGeometry(900, 700); // configure predefined window geometry
   fWebWindow->SetConnLimit(1); // the only connection is allowed
   fWebWindow->SetMaxQueueLength(30); // number of allowed entries in the window queue

   if (mgr) SetGeometry(mgr);
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// destructor

ROOT::Experimental::REveGeomViewer::~REveGeomViewer()
{

}

//////////////////////////////////////////////////////////////////////////////////////////////
/// assign new geometry to the viewer

void ROOT::Experimental::REveGeomViewer::SetGeometry(TGeoManager *mgr)
{
   fGeoManager = mgr;

   fDesc.Build(fGeoManager);

   if (!fGeoManager) return;

   // take maximal setting
   auto maxnodes = fGeoManager->GetMaxVisNodes();
   if (maxnodes > 5000)
      maxnodes = 5000;
   else if (maxnodes < 1000)
      maxnodes = 1000;

   fDesc.SetMaxVisNodes(maxnodes);
   fDesc.SetMaxVisFaces(maxnodes * 100);
   fDesc.SetNSegments(fGeoManager->GetNsegments());
}


/////////////////////////////////////////////////////////////////////////////////
/// Select visible top volume, all other volumes will be disabled

void ROOT::Experimental::REveGeomViewer::SelectVolume(const std::string &volname)
{
   if (!fGeoManager || volname.empty()) {
      fDesc.SelectVolume(nullptr);
   } else {
      fDesc.SelectVolume(fGeoManager->GetVolume(volname.c_str()));
   }
}

/////////////////////////////////////////////////////////////////////////////////
/// Show or update geometry in web window
/// If web browser already started - just refresh drawing like "reload" button does
/// If no web window exists or \param always_start_new_browser configured, starts new window

void ROOT::Experimental::REveGeomViewer::Show(const RWebDisplayArgs &args, bool always_start_new_browser)
{
   auto number = fWebWindow->NumConnections();

   if ((number == 0) || always_start_new_browser) {
      fWebWindow->Show(args);
   } else {
      for (int n=0;n<number;++n)
         WebWindowCallback(fWebWindow->GetConnectionId(n),"RELOAD");
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// convert JSON into stack array

std::vector<int> ROOT::Experimental::REveGeomViewer::GetStackFromJson(const std::string &json)
{
   std::vector<int> *stack{nullptr}, res;

   if (TBufferJSON::FromJSON(stack, json.c_str())) {
      res = *stack;
      delete stack;
   } else {
      R__ERROR_HERE("webeve") << "Fail convert " << json << " into vector<int>";
   }

   return res;
}


//////////////////////////////////////////////////////////////////////////////////////////////
/// receive data from client

void ROOT::Experimental::REveGeomViewer::WebWindowCallback(unsigned connid, const std::string &arg)
{
   printf("Recv %s\n", arg.c_str());

   if ((arg == "CONN_READY") || (arg == "RELOAD")) {

      if (arg == "RELOAD")
         fDesc.Build(fGeoManager);

      std::string sbuf = "DESCR:";
      sbuf.append(TBufferJSON::ToJSON(&fDesc,103).Data());
      printf("Send description %d\n", (int) sbuf.length());
      fWebWindow->Send(connid, sbuf);

      if (!fDesc.HasDrawData())
         fDesc.CollectVisibles();

      auto &json = fDesc.GetDrawJson();
      auto &binary = fDesc.GetDrawBinary();

      printf("Produce JSON %d binary %d\n", (int) json.length(), (int) binary.size());

      fWebWindow->Send(connid, json);

      fWebWindow->SendBinary(connid, &binary[0], binary.size());
   } else if (arg == "QUIT_ROOT") {

      RWebWindowsManager::Instance()->Terminate();

   } else if (arg.compare(0, 7, "SEARCH:") == 0) {
      std::string query = arg.substr(7);

      std::string json;
      std::vector<char> binary;

      auto nmatches = fDesc.SearchVisibles(query, json, binary);

      printf("Searches %s found %d json %d binary %d\n", query.c_str(), nmatches, (int) json.length(), (int) binary.size());

      // send reply with appropriate header - NOFOUND, FOUND0:, FOUND1:
      fWebWindow->Send(connid, json);

      if (binary.size() > 0)
         fWebWindow->SendBinary(connid, &binary[0], binary.size());
   } else if (arg.compare(0,4,"GET:") == 0) {
      // provide exact shape

      auto stack = GetStackFromJson(arg.substr(4));

      auto nodeid = fDesc.FindNodeId(stack);

      std::string json{"SHAPE:"};
      std::vector<char> binary;

      fDesc.ProduceDrawingFor(nodeid, json, binary);

      printf("Produce shape for stack json %d binary %d\n", (int) json.length(), (int) binary.size());

      fWebWindow->Send(connid, json);

      if (binary.size() > 0)
         fWebWindow->SendBinary(connid, &binary[0], binary.size());

   } else if ((arg.compare(0,7,"SETVI0:") == 0) || (arg.compare(0,7,"SETVI1:") == 0)) {
      // change visibility for specified nodeid

      auto nodeid = std::stoi(arg.substr(7));

      bool selected = (arg[5] == '1');

      if (fDesc.ChangeNodeVisibility(nodeid, selected)) {

         // send modified entry only for specified node, when disabled client will automatically remove node from drawing
         std::string json0 = "MODIF:";
         json0.append(TBufferJSON::ToJSON(&fDesc.GetGeomNode(nodeid),103).Data());
         fWebWindow->Send(connid, json0);

         if (selected && fDesc.IsPrincipalNode(nodeid)) {
            // we need to send changes in drawing nodes

            std::string json{"APPND:"};
            std::vector<char> binary;

            fDesc.ProduceDrawingFor(nodeid, json, binary);

            if (binary.size() > 0) {
               printf("Send appending JSON %d binary %d\n", (int) json.length(), (int) binary.size());

               fWebWindow->Send(connid, json);
               fWebWindow->SendBinary(connid, &binary[0], binary.size());
            }
         }
      }
   }
}
