QT += network widgets webengine webenginewidgets

TEMPLATE	= app
LANGUAGE	= C++
MOC_DIR     =.moc
OBJECTS_DIR =.obj

CONFIG	+= qt warn_off thread

INCLUDEPATH += $$(ROOTSYS)/include
DEPENDPATH  += $$(ROOTSYS)/include

TARGET      = qt5web
DESTDIR     = .
PROJECTNAME = RootQt5Web

win32:QMAKE_LFLAGS  += /nodefaultlib:msvcrt /verbose:lib msvcrt.lib

win32:QMAKE_CXXFLAGS  += -MD

# this is necessary to solve error with non-initialized gSystem
unix:QMAKE_CXXFLAGS += -fPIC

FORMS += ExampleWidget.ui

HEADERS += ExampleWidget.h \
           TCanvasWidget.h \
           RCanvasWidget.h \
           RGeomViewerWidget.h

SOURCES += ExampleWidget.cpp \
           TCanvasWidget.cpp \
           RCanvasWidget.cpp \
           RGeomViewerWidget.cpp \
           ExampleMain.cpp
