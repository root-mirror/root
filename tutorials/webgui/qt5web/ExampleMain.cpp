// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <QTimer>
#include <QApplication>

#include "TApplication.h"
#include "TSystem.h"

#include "ExampleWidget.h"

int main(int argc, char **argv)
{
   argc = 1; // hide all additional parameters from ROOT and Qt
   TApplication app("uno", &argc, argv); // ROOT application

   char* argv2[3];
   argv2[0] = argv[0];
   argv2[1] = 0;

   QApplication myapp(argc, argv2); // Qt application

   // let run ROOT event loop
   QTimer timer;
   QObject::connect( &timer, &QTimer::timeout, []() { gSystem->ProcessEvents();  });
   timer.setSingleShot(false);
   timer.start(20);

   // top widget
   ExampleWidget* widget = new ExampleWidget();

   QObject::connect(&myapp, &QApplication::lastWindowClosed, &myapp, &QApplication::quit);

   widget->show();

   int res = myapp.exec();
   return res;
}
