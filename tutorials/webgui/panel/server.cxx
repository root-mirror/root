/// \file
/// \ingroup tutorial_webgui
///  This macro demonstrates simple openui5 panel, shown with RWebWindow
/// \macro_code
///
/// \author Sergey Linev


#include <ROOT/RWebWindow.hxx>
#include "TBufferJSON.h"
#include <vector>
#include <string>

/** Simple structure for ComboBox item */
struct ComboBoxItem {
   std::string fId;
   std::string fName;
   ComboBoxItem() = default;
   ComboBoxItem(const std::string &id, const std::string &name) : fId(id), fName(name) {}
};

/** Full model used to configure openui5 widget */
struct TestPanelModel {
   std::string fSampleText;
   std::vector<ComboBoxItem> fComboItems;
   std::string fSelectId;
   std::string fButtonText;
};

std::shared_ptr<ROOT::Experimental::RWebWindow> window;
std::unique_ptr<TestPanelModel> model;
int sendcnt = 0;

void ProcessConnection(unsigned connid)
{
   printf("connection established %u\n", connid);
   TString json = TBufferJSON::ToJSON(model.get());
   window->Send(connid, std::string("MODEL:") + json.Data());
}

void ProcessCloseConnection(unsigned connid)
{
   printf("connection closed %u\n", connid);
}

void ProcessData(unsigned connid, const std::string &arg)
{
   if (arg == "REFRESH") {
      // send model to client again
      printf("Send model again\n");
      model->fButtonText = Form("Custom button %d", ++sendcnt);
      TString json = TBufferJSON::ToJSON(model.get());
      window->Send(connid, std::string("MODEL:") + json.Data());
   } else if (arg.find("MODEL:") == 0) {
      // analyze new model send from client
      auto m = TBufferJSON::FromJSON<TestPanelModel>(arg.substr(6));
      if (m) {
         printf("New model, selected: %s\n", m->fSelectId.c_str());
         std::swap(model, m);
      } else {
         printf("Fail to decode model: %s\n", arg.c_str());
      }
   }
}

void server()
{
   // prepare model
   model = std::make_unique<TestPanelModel>();
   model->fSampleText = "This is openui5 widget";
   model->fComboItems = {{"item1", "Text 1"}, {"item2", "Text 2"}, {"item3", "Text 3"}, {"item4", "Text 4"}};
   model->fSelectId = "item2";
   model->fButtonText = "Custom button";

   // create window
   window = ROOT::Experimental::RWebWindow::Create();

   // Important - defines name of openui5 widget
   // "localapp" prefix will be point on current directory, where script executed
   // "localapp.view.TestPanel" means file ./view/TestPanel.view.xml will be loaded
   window->SetPanelName("localapp.view.TestPanel");

   // Provide window client version to control browser cache
   // When value changed, URL for JSROOT, UI5 and local files will differ
   // Therefore web browser automatically reload all these files
   // window->SetClientVersion("1.2");

   // these are different callbacks
   window->SetCallBacks(ProcessConnection, ProcessData, ProcessCloseConnection);

   window->SetGeometry(400, 500); // configure window geometry

   window->Show();
}
