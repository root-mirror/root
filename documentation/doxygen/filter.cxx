// @(#)root/test:$name:  $:$id: filter.cxx,v 1.0 exp $
// Author: O.Couet

/// The ROOT doxygen filter implements ROOT's specific directives used to generate
/// the ROOT reference guide.
///
/// ## In the ROOT classes
///
/// ### `Begin_Macro` and `End_Macro`
/// The two tags where used the THtml version to generate images from ROOT code.
/// The genererated picture is inlined exactly at the place where the macro is
/// defined. The Macro can be defined in two way:
///  - by direct in-lining of the the C++ code
///  - by a reference to a C++ file
/// The tag `Begin_Macro` can have the parameter `(source)`. The directive becomes:
/// `Begin_Macro(source)`. This parameter allows to show the macro's code in addition.
///
/// ## In the ROOT tutorials
///
/// ROOT tutorials are also included in the ROOT documentation. The tutorials'
/// macros headers should look like:
///
/// ~~~ {.cpp}
/// \file
/// \ingroup tutorial_hist
/// Getting Contours From TH2D.
///
/// #### Image produced by `.x ContourList.C`
/// The contours values are drawn next to each contour.
/// \macro_image
///
/// #### Output produced by `.x ContourList.C`
/// It shows that 6 contours and 12 graphs were found.
/// \macro_output
///
/// #### `ContourList.C`
/// \macro_code
///
/// \authors  Josh de Bever, Olivier Couet
/// ~~~
///
/// This example shows that three new directives have been implemented:
///
///  1. `\macro_image`
///  The images produced by this macro are shown. A caption can be added to document
///  the pictures: `\macro_image This is a picture`
///
///  2. `\macro_code`
///  The macro code is shown.  A caption can be added: `\macro_code This is code`
///
///  3. `\macro_output`
///  The output produced by this macro is shown. A caption can be added:
///  `\macro_output This the macro output`
///
/// Note that the doxygen directive `\authors` or `\author` must be the last one
/// of the macro header.

#include <unistd.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <stdarg.h>
#include <memory>

using namespace std;

// Auxiliary functions
void   FilterClass();
void   FilterTutorial();
void   GetClassName();
int    NumberOfImages();
string ImagesList(string&);
void   ExecuteMacro();
void   ExecuteCommand(string);
void   ReplaceAll(string&, const string&, const string&);
string StringFormat(const string fmt_str, ...);
bool   EndsWith(string const &, string const &);
bool   BeginsWith(const string&, const string&);

// Global variables.
FILE  *f;              // Pointer to the file being parsed.
char   gLine[255];     // Current line in the current input file
string gFileName;      // Input file name
string gLineString;    // Current line (as a string) in the current input file
string gClassName;     // Current class name
string gImageName;     // Current image name
string gMacroName;     // Current macro name
string gCwd;           // Current working directory
string gOutDir;        // Output directory
string gSourceDir;     // Source directory
string gOutputName;    // File containing a macro std::out
bool   gHeader;        // True if the input file is a header
bool   gSource;        // True if the input file is a source file
bool   gPython;        // True if the input file is a Python script.
bool   gImageSource;   // True the source of the current macro should be shown
int    gInMacro;       // >0 if parsing a macro in a class documentation.
int    gImageID;       // Image Identifier.
int    gMacroID;       // Macro identifier in class documentation.
int    gShowTutSource; // >0 if the tutorial source code should be shown


////////////////////////////////////////////////////////////////////////////////
/// Filter ROOT files for Doxygen.

int main(int argc, char *argv[])
{
   // Initialisation
   gFileName      = argv[1];
   gHeader        = false;
   gSource        = false;
   gPython        = false;
   gImageSource   = false;
   gInMacro       = 0;
   gImageID       = 0;
   gMacroID       = 0;
   gOutputName    = "stdout.dat";
   gShowTutSource = 0;
   if (EndsWith(gFileName,".cxx")) gSource = true;
   if (EndsWith(gFileName,".h"))   gHeader = true;
   if (EndsWith(gFileName,".py"))  gPython = true;
   GetClassName();

   // Retrieve the current working directory
   int last = gFileName.rfind("/");
   gCwd     = gFileName.substr(0,last);

   // Retrieve the output directory
   gOutDir = getenv("DOXYGEN_OUTPUT_DIRECTORY");
   ReplaceAll(gOutDir,"\"","");

   // Retrieve the source directory
   gSourceDir = getenv("DOXYGEN_SOURCE_DIRECTORY");
   ReplaceAll(gSourceDir,"\"","");

   // Open the input file name.
   f = fopen(gFileName.c_str(),"r");

   if (gFileName.find("tutorials") != string::npos) FilterTutorial();
   else                                             FilterClass();
}

////////////////////////////////////////////////////////////////////////////////
/// Filter ROOT class for Doxygen.

void FilterClass()
{
   // File for inline macros.
   FILE *m = 0;

   // File header.
   if (gHeader) {
      while (fgets(gLine,255,f)) {
         gLineString = gLine;
         printf("%s",gLineString.c_str());
      }
      fclose(f);
      return;
   }

   // Source file.
   if (gSource) {
      size_t spos = 0;
      while (fgets(gLine,255,f)) {
         gLineString = gLine;

         if (gInMacro && gLineString.find("End_Macro") != string::npos) {
            ReplaceAll(gLineString,"End_Macro","");
            gImageSource = false;
            gInMacro = 0;
            spos = 0;
            if (m) {
               fclose(m);
               m = 0;
               ExecuteCommand(StringFormat("root -l -b -q \"makeimage.C(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",true,false)\""
                                              , StringFormat("%s_%3.3d.C", gClassName.c_str(), gMacroID).c_str()
                                              , StringFormat("%s_%3.3d.png", gClassName.c_str(), gImageID).c_str()
                                              , gOutDir.c_str()));
               ExecuteCommand(StringFormat("rm %s_%3.3d.C", gClassName.c_str(), gMacroID));
            }
         }

         if (gInMacro) {
            if (spos) gLineString = gLineString.substr(spos);
            if (gInMacro == 1) {
               if (EndsWith(gLineString,".C\n") || (gLineString.find(".C(") != string::npos)) {
                  ExecuteMacro();
                  gInMacro++;
               } else {
                  gMacroID++;
                  m = fopen(StringFormat("%s_%3.3d.C", gClassName.c_str(), gMacroID).c_str(), "w");
                  if (m) fprintf(m,"%s",gLineString.c_str());
                  if (BeginsWith(gLineString,"{")) {
                     if (gImageSource) {
                        ReplaceAll(gLineString,"{"
                                                , StringFormat("\\include %s_%3.3d.C"
                                                , gClassName.c_str()
                                                , gMacroID));
                     } else {
                        gLineString = "\n";
                     }
                  }
                  gInMacro++;
               }
            } else {
               if (m) fprintf(m,"%s",gLineString.c_str());
               if (BeginsWith(gLineString,"}")) {
                  ReplaceAll(gLineString,"}", StringFormat("\\image html pict1_%s_%3.3d.png", gClassName.c_str(), gImageID));
               } else {
                  gLineString = "\n";
               }
               gInMacro++;
            }
         }

         if (gLineString.find("Begin_Macro") != string::npos &&
             gLineString.find("End_Macro") == string::npos) {
            if (BeginsWith(gLineString, "///")) {
               spos = gLineString.find_first_not_of(' ', 3);
            }
            if (gLineString.find("source") != string::npos) gImageSource = true;
            gImageID++;
            gInMacro++;
            gLineString = "\n";
         }

         size_t l = gLineString.length();
         size_t b = 0;
         do {
            size_t e = gLineString.find('\n', b);
            if (e != string::npos) e++;
            if (spos) printf("%-*s%s", (int)spos, "///",
                              gLineString.substr(b, e - b).c_str());
            else printf("%s", gLineString.substr(b, e - b).c_str());
            b = e;
         } while (b < l);
      }
      fclose(f);
      return;
   }

   // Output anything not header nor source
   while (fgets(gLine,255,f)) {
      gLineString = gLine;
      printf("%s",gLineString.c_str());
   }
   fclose(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Filter ROOT tutorials for Doxygen.

void FilterTutorial()
{
   // File for inline macros.
   FILE *m = 0;

   // Extract the macro name
   int i1      = gFileName.rfind('/')+1;
   int i2      = gFileName.rfind('C');
   gMacroName  = gFileName.substr(i1,i2-i1+1);
   gImageName  = StringFormat("%s.png", gMacroName.c_str()); // Image name
   gOutputName = StringFormat("%s.out", gMacroName.c_str()); // output name

   // Parse the source and generate the image if needed
   while (fgets(gLine,255,f)) {
      gLineString = gLine;

      // \macro_image found
      if (gLineString.find("\\macro_image") != string::npos) {
         bool nobatch = (gLineString.find("(nobatch)") != string::npos);
         ReplaceAll(gLineString,"(nobatch)","");
         if (gPython) {
            if (nobatch) {
               ExecuteCommand(StringFormat("./makeimage.py %s %s %s 0 1 0",
                                          gFileName.c_str(), gImageName.c_str(), gOutDir.c_str()));
            } else {
               ExecuteCommand(StringFormat("./makeimage.py %s %s %s 0 1 1",
                                          gFileName.c_str(), gImageName.c_str(), gOutDir.c_str()));
            }
         } else {
            if (nobatch) {
               ExecuteCommand(StringFormat("root -l -q \"makeimage.C(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",false,false)\"",
                                            gFileName.c_str(), gImageName.c_str(), gOutDir.c_str()));
            } else {
               ExecuteCommand(StringFormat("root -l -b -q \"makeimage.C(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",false,false)\"",
                                            gFileName.c_str(), gImageName.c_str(), gOutDir.c_str()));
            }
         }
         ReplaceAll(gLineString, "\\macro_image", ImagesList(gImageName));
         remove(gOutputName.c_str());
      }

      // \macro_code found
      if (gLineString.find("\\macro_code") != string::npos) {
         gShowTutSource = 1;
         m = fopen(StringFormat("%s/macros/%s",gOutDir.c_str(),gMacroName.c_str()).c_str(), "w");
         ReplaceAll(gLineString, "\\macro_code", StringFormat("\\include %s",gMacroName.c_str()));
      }

      // notebook found
      if (gLineString.find("\\notebook") != string::npos) {
         ExecuteCommand(StringFormat("python converttonotebook.py %s %s/notebooks/",
                                          gFileName.c_str(), gOutDir.c_str()));
         ReplaceAll(gLineString, "\\notebook -js", StringFormat( "\\htmlonly <a href=\"http://nbviewer.jupyter.org/url/root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src= notebook.gif alt=\"View in nbviewer\" style=\"height:1em\" ></a> <a href=\"https://cern.ch/swanserver/cgi-bin/go?projurl=https://root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src=\"http://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png\"  alt=\"Open in SWAN\" style=\"height:1em\" ></a> \\endhtmlonly", gMacroName.c_str() , gMacroName.c_str()) );         
         
         ReplaceAll(gLineString, "\\notebook -nodraw", StringFormat( "\\htmlonly <a href=\"http://nbviewer.jupyter.org/url/root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src= notebook.gif alt=\"View in nbviewer\" style=\"height:1em\" ></a> <a href=\"https://cern.ch/swanserver/cgi-bin/go?projurl=https://root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src=\"http://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png\"  alt=\"Open in SWAN\" style=\"height:1em\" ></a> \\endhtmlonly", gMacroName.c_str() , gMacroName.c_str()) );
         
         ReplaceAll(gLineString, "\\notebook", StringFormat( "\\htmlonly <a href=\"http://nbviewer.jupyter.org/url/root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src= notebook.gif alt=\"View in nbviewer\" style=\"height:1em\" ></a> <a href=\"https://cern.ch/swanserver/cgi-bin/go?projurl=https://root.cern.ch/doc/master/notebooks/%s.nbconvert.ipynb\" target=\"_blank\"><img src=\"http://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png\"  alt=\"Open in SWAN\" style=\"height:1em\" ></a> \\endhtmlonly", gMacroName.c_str() , gMacroName.c_str()) );

         
      }
      // \macro_output found
      if (gLineString.find("\\macro_output") != string::npos) {
         if (!gPython) ExecuteCommand(StringFormat("root -l -b -q %s", gFileName.c_str()).c_str());
         else          ExecuteCommand(StringFormat("python %s", gFileName.c_str()).c_str());
         rename(gOutputName.c_str(), StringFormat("%s/macros/%s",gOutDir.c_str(), gOutputName.c_str()).c_str());
         ReplaceAll(gLineString, "\\macro_output", StringFormat("\\include %s",gOutputName.c_str()));
      }

      // \author is the last comment line.
      if (gLineString.find("\\author")  != string::npos) {
         if (gPython) printf("%s",StringFormat("%s \n## \\cond \n",gLineString.c_str()).c_str());
         else         printf("%s",StringFormat("%s \n/// \\cond \n",gLineString.c_str()).c_str());
         if (gShowTutSource == 1) gShowTutSource = 2;
      } else {
         printf("%s",gLineString.c_str());
         if (m && gShowTutSource == 2) fprintf(m,"%s",gLineString.c_str());
      }
   }

   if (m) {
      if (gPython) printf("## \\endcond \n");
      else         printf("/// \\endcond \n");
      fclose(m);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the class name.

void GetClassName()
{
   int i1 = 0;
   int i2 = 0;

   FILE *f = fopen(gFileName.c_str(),"r");

   // File header.
   if (gHeader) {
      while (fgets(gLine,255,f)) {
         gLineString = gLine;
         if (gLineString.find("ClassDef") != string::npos) {
            i1         = gLineString.find("(")+1;
            i2         = gLineString.find(",")-1;
            gClassName = gLineString.substr(i1,i2-i1+1);
            fclose(f);
            return;
         }
      }
   }

   // Source file.
   if (gSource) {
      while (fgets(gLine,255,f)) {
         gLineString = gLine;
         if (gLineString.find("ClassImp") != string::npos ||
             gLineString.find("NamespaceImp") != string::npos) {
            i1         = gLineString.find("(")+1;
            i2         = gLineString.find(")")-1;
            gClassName = gLineString.substr(i1,i2-i1+1);
            fclose(f);
            return;
         }
      }
   }

   fclose(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the macro in gLineString and produce the corresponding picture.

void ExecuteMacro()
{
   // Name of the next Image to be generated
   gImageName = StringFormat("%s_%3.3d.png", gClassName.c_str(), gImageID);

   // Retrieve the macro to be executed.
   if (gLineString.find("../../..") != string::npos) {
      ReplaceAll(gLineString,"../../..", gSourceDir.c_str());
   } else {
      gLineString.insert(0, StringFormat("%s/../doc/macros/",gCwd.c_str()));
   }
   int i1     = gLineString.rfind('/')+1;
   int i2     = gLineString.rfind('C');
   gMacroName = gLineString.substr(i1,i2-i1+1);

   // Build the ROOT command to be executed.
   gLineString.insert(0, StringFormat("root -l -b -q \"makeimage.C(\\\""));
   size_t l = gLineString.length();
   gLineString.replace(l-1,1,StringFormat("\\\",\\\"%s\\\",\\\"%s\\\",true,false)\"", gImageName.c_str(), gOutDir.c_str()));

   // Execute the macro
   ExecuteCommand(gLineString);

   // Inline the directives to show the picture and/or the code
   if (gImageSource) {
      gLineString = StringFormat("\\include %s\n\\image html pict1_%s\n", gMacroName.c_str(), gImageName.c_str());
   } else {
      gLineString = StringFormat("\n\\image html pict1_%s\n", gImageName.c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a command making sure stdout will not go in the doxygen file.

void ExecuteCommand(string command)
{
   int o = dup(fileno(stdout));
   freopen(gOutputName.c_str(),"a",stdout);
   system(command.c_str());
   dup2(o,fileno(stdout));
   close(o);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the number of images in NumberOfImages.dat after makeimage.C is executed.

int NumberOfImages()
{
   int ImageNum;
   FILE *f = fopen("NumberOfImages.dat", "r");
   fscanf(f, "%d", &ImageNum);
   fclose(f);
   remove("NumberOfImages.dat");
   return ImageNum;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace all instances of a string with another string.

void ReplaceAll(string& str, const string& from, const string& to) {
   if (from.empty()) return;
   string wsRet;
   wsRet.reserve(str.length());
   size_t start_pos = 0, pos;
   while ((pos = str.find(from, start_pos)) != string::npos) {
      wsRet += str.substr(start_pos, pos - start_pos);
      wsRet += to;
      pos += from.length();
      start_pos = pos;
   }
   wsRet += str.substr(start_pos);
   str.swap(wsRet);
}

////////////////////////////////////////////////////////////////////////////////
/// std::string formatting like sprintf.

string StringFormat(const string fmt_str, ...) {
   int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
   string str;
   unique_ptr<char[]> formatted;
   va_list ap;
   while (1) {
      formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
      strcpy(&formatted[0], fmt_str.c_str());
      va_start(ap, fmt_str);
      final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
      va_end(ap);
      if (final_n < 0 || final_n >= n) n += abs(final_n - n + 1);
      else break;
   }
   return string(formatted.get());
}

////////////////////////////////////////////////////////////////////////////////
/// Return the image list after a tutorial macro execution.

string ImagesList(string& name) {

   int N = NumberOfImages();

   char val[300];
   int len = 0;
   for (int i = 1; i <= N; i++){
      if (i>1) sprintf(&val[len]," \n/// \\image html pict%d_%s",i,name.c_str());
      else     sprintf(&val[len],"\\image html pict%d_%s",i,name.c_str());
      len = (int)strlen(val);
   }
   return (string)val;
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a string ends with another string.

bool EndsWith(string const &fullString, string const &ending) {
   if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
   } else {
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a string begins with another string.

bool BeginsWith(const string& haystack, const string& needle) {
   return needle.length() <= haystack.length() && equal(needle.begin(), needle.end(), haystack.begin());
}
