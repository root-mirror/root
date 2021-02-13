## OPENUI5 Panel example

This is simplest way to use openui5 widget with RWebWindow.

It is normal xml::View, but controller should be derived from rootui5/panel/Controller.
This class provides methods, which simplify handling of communication between server and client

First of all, when creating RWebWindow, one should configure panel name. Like:

     auto win = ROOT::Experimental::RWebWindow::Create();

     win->SetPanelName("localapp.view.TestPanel");

Namespace "localapp" in this case corresponds to openui5 files, which will be loaded from current directory.
Therefore `"localapp.view.TestPanel"` means view, which will be loaded from `./view/TestPanel.view.xml` file.

Controller is configured in the XML file and called `"localapp.controller.TestPanel"`.
Means it will be loaded from `./controller/TestPanel.controller.js` file.

In the controller one use `onPanelInit` and `onPanelExit` methods to handle initialization and close of widget.
Method `panelSend` should be used to send string data to the server.

