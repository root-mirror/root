ROOTSL
======

This is a Spotlight plugin that allows ROOT files to be indexed by SL.
Once indexed SL can find ROOT files based on the names and titles of the
objects in the files.

Spotlight is available on MacOS X since version 10.4 (Tiger). To use SL
select the SL icon on the top right of the menubar and type in a search string.

Get the binary for the ROOTSL plugin from:

   ftp://root.cern/root/ROOTSL.tgz

To install the plugin, after untarring the above file, just double click the
ROOTSL.mdimporter icon. If you have no admin rights you will be asked for
the admin password. The plugin will be installed in /Library/Spotlight.
You can also install the plugin in your private area by dragging the
plugin to ~/Library/Spotlight. You may have to create that folder if
it does not exist. Once installed you have to tell SL to import existing
files by executing:
   /usr/bin/mdimport -r [~]/Library/Spotlight/ROOTSL.mdimporter
Spotlight will then, in the background, index all *.root files.

To build from source, get it from svn using:

   git clone http://root.cern/git/root.git root
   cd root/misc/rootsl

Open the ROOTSL project in Xcode and click on "Build" (make sure the Active
Build Configuration is set to "Release"). A command line short cut to open
the Xcode project is to type "open ROOTSL.xcodeproj" in the Terminal app.
Copy the resulting plugin from build/Release to the /Library/Spotlight
directory by double clicking the icon, or by typing in the shell
"open ROOTSL.mdimporter".

Cheers, Fons.
