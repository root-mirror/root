ROOTQL
======

This is a Quick Look plugin that allows quick inspection of the content
of ROOT files.

Quick Look is available on MacOS X since version 10.5 (Leopard). To use QL
select a file icon in the Finder and hit the space bar. For all file types
supported by QL you will get a window showing the file content, for file types
not supported you will get a generic window showing some basic file info.

The idea of QL is that file content can be shown without the heavy application
startup process. Generating a QL view of a ROOT file depends on the size of the
file and number of keys, but generally it is a quick operation.

Get the binary for the ROOTQL plugin from:

   ftp://root.cern/root/ROOTQL.tgz

To install the plugin, after untarring the above file, just drag the
ROOTQL.qlgenerator icon to either /Library/QuickLook or to ~/Library/QuickLook.
You may have to create that folder if it does not exist. Once installed
you may have to refresh the QL plugin cache by executing:
   /usr/bin/qlmanage -r

To build from source, get it from git using:

   git clone http://root.cern/git/root.git root
   cd root/misc/rootql

Open the ROOTQL project in Xcode and click on "Build" (make sure the Active
Build Configuration is set to "Release"). A command line short cut to open
the Xcode project is to type "open ROOTQL.xcodeproj" in the Terminal app.
Move the resulting plugin from the build/Release directory to either
the /Library/QuickLook or ~/Library/QuickLook directory.

Cheers, Fons.
