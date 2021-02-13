# qt5web example

Demonstration how different ROOT web-based widgets can be embed into qt application
To compile example, just call `make` in the directory.
ROOT should be compiled with configured `-Dqt5web=ON`.
As a result, `qt5web` application should be created.
It includes `QTabWidget` with three tabs - standard Qt5 components for demo,
TCanvas and RCanvas. Both show different histograms drawing.

## How to include RCanvas into other Qt-based project

Most easy way - just include `RCanvasWidget.h` and `RCanvasWidget.cpp` files
in the project and let compile, linking with ROOT basic libraries `root-config --libs` plus `-lROOTWebDisplay -lROOTGpadv7`. `RCanvasWidget` is just `QWidget` which internally embed `RCanvas` drawing. See `ExampleWidget.ui` file how to embed such custom widget in normal qt ui file.

To let ROOT work inside Qt5 event loop, one should instantiate `TApplication` object and
regularly call `gSystem->ProcessEvents()` - see how it is done in `ExampleMain.cpp`.

Author: Sergey Linev
