// This script generates the html pages for the ROOT tutorials hierarchy.
// It creates $ROOTSYS/htmldoc if not already there
// It creates $ROOTSYS/htmldoc/tutorials if not already there
// It creates $ROOTSYS/htmldoc/tutorials/index.html (index to all directory tutorials)
// It creates $ROOTSYS/htmldoc/tutorials/dir/index.html with index of tutorials in dir
// It creates $ROOTSYS/htmldoc/tutorials/dir/*.C.html with one html file for each tutorial
// Author: Rene Brun

#include "THtml.h"
#include "TDocOutput.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TList.h"
#include <iostream>
#include <fstream>

using namespace std;

void scandir(THtml& html, const char *dir, const char *title, TObjLink* toplnk);


void AppendLink(TString& links, int id, const TNamed* n)
{
   // Used to construct the context links (prev/up/next) at the top of each html page.
   // "links" contains the html code that AppendLink() adds to.
   // "id"    contains the pass (id == 0: prev, 1: up, 2: next)
   // "n"     points to a TNamed that has the URL (without trailing ".html") as the name
   //         and the title as the title.

   static const char* tag[] = {"&lt;", "^", "&gt;"}; // "<" "^" ">"
   static const char* name[] = {"prev", "up", "next"}; // for CSS
   static const TNamed emptyName;

   // for CSS
   TString arrowid("contextheadarrow_"); arrowid += name[id];
   TString entryid("contextheadentry_"); entryid += name[id];
   TString ardivid("divcontextheadar_"); ardivid += name[id];
   if (!n) n = &emptyName;

   TString title;
   // format for prev, next: script: title (if not empty)
   // format for up: title
   if (id != 1) title = n->GetName();
   if (title.Length()) title += ": ";
   title += n->GetTitle();
   const char* mytag = tag[id];
   // if we don't have a link we don't write <^>
   if (!n->GetName()[0]) mytag = "";
   // td for "<" "^" ">"
   TString arrow = TString::Format("<td id=\"%s\"><div id=\"%s\">"
                                   "<a class=\"contextheadarrow\" href=\"%s.html\">%s</a></div></td>",
                                   arrowid.Data(), ardivid.Data(), n->GetName(), mytag);
   // td for the text
   TString entry = TString::Format("<td class=\"contextheadentry\" id=\"%s\">"
                                   "<a class=\"contextheadentry\" href=\"%s.html\">%s</a></td>",
                                   entryid.Data(), n->GetName(), title.Data());
   if (id == 2)
      // "next" has the text first
      links += entry + arrow;
   else 
      links += arrow + entry;
}


void MakeTopLinks(TString &links, const char* name, const char* title, const char* upLink, const char* upTitle,
                  TObjLink *lnk, const char* dir)
{
// Create the html code for the navigation box at the top of each page,
// showing a link to the previous and next tutorial and to the upper level.
// "links"   will hold the html code after the function returns,
// "title"   is the title for the page, displayed below the naviation links
// "upLink"  is the URL for the up link,
// "upTitle" is the text to display for the link.
// "lnk"     has TNamed as prev and next, with name as URL and title as
//           title for the previous and next links.

   links = "<div id=\"toplinks\"><div class=\"descrhead\">"
           "<table class=\"descrtitle\" id=\"contexttitle\"><tr class=\"descrtitle\">";
   TObjLink *prevlnk = lnk ? lnk->Prev() : 0;
   TObjLink *nextlnk = lnk ? lnk->Next() : 0;

   TNamed* prevname = prevlnk ? (TNamed*)prevlnk->GetObject() : 0;
   AppendLink(links, 0, prevname);

   TNamed upname;
   if (upLink && upLink[0])
      upname.SetNameTitle(upLink, upTitle);
   AppendLink(links, 1, &upname);

   TNamed* nextname = nextlnk ? (TNamed*)nextlnk->GetObject() : 0;
   AppendLink(links, 2, nextname);

   links += TString("</tr></table></div><h1 class=\"convert\">") + title + "</h1></div>\n";
   TString suburl = dir;
   TString subtitle = dir;
   if (name) {
      if (!subtitle.EndsWith("/")) {
         subtitle += '/';
      }
      subtitle += TString(name);
      suburl = subtitle + "?view=markup";
   }
   links += TString::Format("<div class=\"location\"><h2>From <a href=\"http://root.cern/viewvc/trunk/tutorials/%s\">$ROOTSYS/tutorials/%s</a></h2></div>",
                            suburl.Data(), subtitle.Data());
}

void writeHeader(THtml& html, ostream& out, const char *title, const char* relPath="../") {
   // Write the html file header
   // "html"    THtml object to use
   // "out"     where to write to
   // "title"   title to display in the browser's bar
   // "relPath" relative path to the root of the documentation output,
   // i.e. where to find ROOT.css generated below by THtml::CreateAuxiliaryFiles()

   TDocOutput docout(html);
   docout.WriteHtmlHeader(out, title, relPath);
}

void writeTrailer(THtml& html, ostream& out) {
   // Write the html file trailer
   // "html"    THtml object to use
   // "out"     where to write to

   TDocOutput docout(html);
   docout.WriteHtmlFooter(out);
}

void writeItem(ostream& out, Int_t numb, const char *ref, const char *name, const char *title, Bool_t isnew) {
   // Write a list entry in the directory index.
   // "out"     where to write to
   // "numb"    number of the current line
   // "ref"     URL of the line's page
   // "name"    name to display for the link
   //
   const char *imagenew = "";
   cout << "writingItem: " << numb << ", ref=" << ref << ", name=" << name << ", title=" << title << endl;
   if (isnew) imagenew = " <img src=\"http://root.cern/root/images/new01.gif\" alt=\"new\" align=\"top\" />";
   out << "<li class=\"idxl" << numb%2 << "\">";
   out << "<a href=\"" << ref << "\"><span class=\"typename\">" << numb << ". " << name << "</span></a> "
       << title << imagenew << "</li>" << endl;
}

void writeItemDir(THtml& html, ostream& out, TObjLink* lnk) {
   // Process a tutorial directory: add a list entry for it in the
   // topmost index and add all tutorials in this directory.
   // "html"   THtml object to use
   // "out"    where to write the html code to
   // "lnk"    a TObjLink with Prev() and Next() pointing to TNamed that
   //          hold the URL (name) and the title (title) for the previous
   //          and next tutorial directory.

   static int n=0;
   const char *dir = lnk->GetObject()->GetName();
   const char *title = lnk->GetObject()->GetTitle();
   out << "<li class=\"idxl" << (n++)%2 << "\"><a href=\"" << dir << "/index.html\">"
       << "<span class=\"typename\">" << dir << "</span></a>" << title << "</li>" << endl;

   scandir(html, dir, title, lnk);
}

void writeTutorials(THtml& html) {
   // Process all tutorials by looking over the directories in
   // $ROOTSYS/tutorials, generating index pages showing them
   // and index pages for each directory showing its tutorials,
   // and by converting all tutorials to html code, including
   // the graphics output where possible. The latter is done
   // using THtml::Convert().

   // tutorials and their titles; ordered by "significance"
   const char* tutorials[][2] = {
      {"hist",     "Histograms"},
      {"graphics", "Basic Graphics"}, 
      {"graphs",   "TGraph, TGraphErrors, etc"}, 
      {"gui",      "Graphics User Interface"}, 
      {"fit",      "Fitting tutorials"}, 
      {"fitsio",   "CFITSIO interface"}, 
      {"io",       "Input/Output"}, 
      {"tree",     "Trees I/O, Queries, Graphics"}, 
      {"math",     "Math tutorials"}, 
      {"matrix",   "Matrix packages tutorials"}, 
      {"geom",     "Geometry package"}, 
      {"gl",       "OpenGL examples"}, 
      {"eve",      "Event Display"}, 
      {"fft",      "Fast Fourier Transforms"}, 
      {"foam",     "TFoam example"}, 
      {"image",    "Image Processing"}, 
      {"mlp",      "Neural Networks"}, 
      {"net",      "Network, Client/server"}, 
      {"physics",  "Physics misc"}, 
      {"proof",    "PROOF tutorials"}, 
      {"pyroot",   "Python-ROOT"}, 
      {"pythia",   "Pythia event generator"}, 
      {"quadp",    "Quadratic Programming package"}, 
      {"roofit",   "RooFit tutorials"}, 
      {"roostats", "Roostats tutorials"}, 
      {"spectrum", "Peak Finder, Deconvolutions"}, 
      {"splot",    "TSPlot example"}, 
      {"sql",      "SQL Data Bases interfaces"}, 
      {"thread",   "Multi-Threading examples"}, 
      {"unuran",   "The Unuran package"}, 
      {"xml",      "XML tools"},
      {0, 0}
   };

   // the output file for the directory index
   ofstream fptop("htmldoc/tutorials/index.html");
   writeHeader(html, fptop,"ROOT Tutorials");
   TString topLinks;
   MakeTopLinks(topLinks, 0, "ROOT Tutorials", "../index", "ROOT", 0, "");
   fptop << topLinks << endl;
   fptop << "<ul id=\"indx\">" << endl;

   // Iterate over all tutorial directories.
   // We need prev and next, so keep prev, curr, and next
   // in a TList containing three TNamed, and sweep through
   // the char array tutorials.
   TList contextList;
   TNamed prev;
   TNamed curr(tutorials[0][0], tutorials[0][1]);
   TNamed next(tutorials[1][0], tutorials[1][1]);
   contextList.AddLast(&prev);
   contextList.AddLast(&curr);
   contextList.AddLast(&next);
   TObjLink* lnk = contextList.FirstLink();
   lnk = lnk->Next(); // "curr" is the second link
   const char** iTut = tutorials[2];
   while (iTut[0]) {
      writeItemDir(html, fptop, lnk);
      prev = curr;
      curr = next;
      next.SetNameTitle(iTut[0], iTut[1]);
      ++iTut; // skip name
      ++iTut; // skip title
   }

   fptop << "</ul>" << endl;
   fptop << "<p><a href=\"http://root.cern/drupal/content/downloading-root\">Download ROOT</a> and run the tutorials in $ROOTSYS/tutorials yourself!</p>" << endl;
   writeTrailer(html, fptop);
}

void GetMacroTitle(const char *fullpath, TString &comment, Bool_t &compile) {
   // Find the best line with a title by scanning the first 50 lines of a macro.
   // "fullpath" location of the macro
   // "comment"  is set to the comment (i.e. title) found in the macro
   // "compile"  is set to true if the macro should be compiled, i.e. the
   //            title line starts with "//+ " (note the space)
   compile = kFALSE;
   FILE *fp = fopen(fullpath,"r");
   char line[250];
   int nlines = 0;
   while (fgets(line,240,fp)) {
      nlines++;
      char *com = strstr(line,"//");
      if (com) {
         if (strstr(line,"Author")) continue;
         if (strstr(line,"@(#)")) continue;
         if (strstr(line,"****")) continue;
         if (strstr(line,"////")) continue;
         if (strstr(line,"====")) continue;
         if (strstr(line,"....")) continue;
         if (strstr(line,"----")) continue;
         if (strstr(line,"____")) continue;
         if (strlen(com+1)  < 5)  continue;
         if (!strncmp(com,"//+ ", 4)) {
            compile = kTRUE;
            com += 2; // skip "+ ", too.
         }
         comment = com+2;
         break;
      }
      if (nlines > 50) break;
   }
   fclose(fp);
}

Bool_t IsNew(const char *filename) {
   // Check if filename in SVN is newer than 6 months
   gSystem->Exec(Form("svn info %s > MakeTutorials-tmp.log",filename));
   FILE *fpdate = fopen("MakeTutorials-tmp.log","r");
   char line[250];
   Bool_t isnew = kFALSE;
   TDatime today;
   Int_t now = 365*(today.GetYear()-1)+12*(today.GetMonth()-1) + today.GetDay();
   Int_t year,month,day;
   while (fgets(line,240,fpdate)) {
      const char *com = strstr(line,"Last Changed Date: ");
      if (com) {
         sscanf(&com[19],"%d-%d-%d",&year,&month,&day);
         Int_t filedate = 365*(year-1) + 12*(month-1) + day; //see TDatime::GetDate
         if (now-filedate< 6*30) isnew = kTRUE;
         break;
      } 
   }
   fclose(fpdate);
   gSystem->Unlink("MakeTutorials-tmp.log");
   return isnew;
}

Bool_t CreateOutput_Dir(const char* dir) {
   // Whether THtml::Convert() should run the tutorials in the
   // directory "dir" and store their output

   if (strstr(dir,"net")) return kFALSE;
   if (strstr(dir,"xml")) return kFALSE;
   if (strstr(dir,"sql")) return kFALSE;
   if (strstr(dir,"proof")) return kFALSE;
   if (strstr(dir,"foam"))  return kFALSE;
   if (strstr(dir,"unuran"))  return kFALSE;
   if (strstr(dir,"roofit"))  return kFALSE;
   if (strstr(dir,"thread"))  return kFALSE;
   return kTRUE;
}
Bool_t CreateOutput_Tutorial(const char* tut) {
   // Whether THtml::Convert() should run the tutorial "tut"
   // and store its output

   static const char* vetoed[] = {
      "geodemo",
      "peaks2",
      "testUnfold",
      "readCode",
      "importCode",
      "hadd",
      "line3Dfit",
      "gtime",
      "games",
      "guiWithCINT",
      "Qt",
      "rs401d_FeldmanCousins",
      "graph_edit_playback",
      "fitpanel_playback",
      "guitest_playback",
      "geom_cms_playback",
      "gviz3d.C",
      0
   };

   for (const char** iVetoed = vetoed; *iVetoed; ++iVetoed)
      if (strstr(tut, *iVetoed))
         return kFALSE;

   return kTRUE;
}

void scandir(THtml& html, const char *dir, const char *title, TObjLink* toplnk) {
   // Process a directory containing tutorials by converting all tutorials to
   // html and creating an index of all tutorials in the directory.

   TString fullpath("htmldoc/tutorials/");
   fullpath += dir;
   if (!gSystem->OpenDirectory(fullpath)) gSystem->MakeDirectory(fullpath);
   fullpath += "/index.html";
   // The index for the current directory
   ofstream fpind(fullpath);
   writeHeader(html, fpind, title, "../../");

   TString topLinks;
   // Creates links to prev: "hist.html", up: ".html", next: "graph.html".
   MakeTopLinks(topLinks, 0, title, ".", "ROOT Tutorials", toplnk, dir);
   // But we need links to prev: "../hist/index.html", up: "../index.html", next: "graph/index.html",
   // so the following works:
   topLinks.ReplaceAll("href=\"", "href=\"../");
   topLinks.ReplaceAll("href=\"../http://", "href=\"http://");
   topLinks.ReplaceAll("href=\"../https://", "href=\"https://");
   topLinks.ReplaceAll(".html\"", "/index.html\"");
   // Also prepend "ROOT Tutorials" to the current title:
   topLinks.ReplaceAll("<h1 class=\"convert\">", "<h1 class=\"convert\">ROOT Tutorials: ");
   fpind << topLinks << endl;
   fpind << "<ul id=\"indx\">" << endl;

   TString outpath("htmldoc/tutorials/");
   outpath += dir;
   TString inpath("$ROOTSYS/tutorials/");
   inpath += dir;
   inpath += "/";
   gSystem->ExpandPathName(inpath);
   void *thedir = gSystem->OpenDirectory(inpath);
   if (!thedir) {
      printf("MakeTutorials.C: error opening directory %s", inpath.Data());
      return;
   }
   const char *direntry;
   THashList h;
   while ((direntry = gSystem->GetDirEntry(thedir))) {
      if(*direntry =='.') continue;
      const char *CC = strstr(direntry,".C");
      // must end on ".C"
      if (!CC || *(CC+2)) continue;
      // do not even document these; they are part of another tutorial:
      if(strstr(direntry,"h1anal")) continue;
      if(strstr(direntry,"hsimpleProxy")) continue;
      if(strstr(direntry,"tv3")) continue;
      if(strstr(direntry,"tvdemo")) continue;
      if(strstr(direntry,"na49")) continue;
      if(strstr(direntry,"fit1_C")) continue;
      if(strstr(direntry,"c1.C")) continue;
      if(strstr(direntry,"MDF.C")) continue;
      if(strstr(direntry,"cms_calo_detail")) continue;
      TString atut(inpath + direntry);
      TString comment;
      Bool_t compile;
      GetMacroTitle(atut,comment, compile);
      TNamed *named = new TNamed(direntry,comment.Data());
      if (compile) named->SetBit(BIT(14));
      h.Add(named);
   }
   h.Sort();
   int numb = 0;
   TObjLink *lnk = h.FirstLink();
   while (lnk) {
      TNamed* named = (TNamed*)lnk->GetObject();
      Bool_t compile = named->TestBit(BIT(14));
      direntry = named->GetName();
      TString atut(inpath + direntry);
      numb++;
      TString iname(direntry);
      iname += ".html";
      writeItem(fpind, numb, iname, direntry, named->GetTitle(), IsNew(atut));
      Int_t includeOutput = THtml::kNoOutput;
      if (!gROOT->IsBatch()) {
         if (compile)
            includeOutput = THtml::kCompiledOutput;
         else
            includeOutput = THtml::kInterpretedOutput;
         includeOutput |= THtml::kSeparateProcessOutput;
      }
      if (!CreateOutput_Dir(dir) || !CreateOutput_Tutorial(direntry))
         includeOutput = THtml::kNoOutput;

      TString links;
      TString tutTitle(named->GetName());
      tutTitle += ": ";
      tutTitle += named->GetTitle();
      MakeTopLinks(links,named->GetName(),tutTitle,"index",title,lnk, dir);
      html.Convert(atut,named->GetTitle(),outpath,"../../",includeOutput,links);
      gROOT->GetListOfCanvases()->Delete();
      gROOT->CloseFiles();
      gROOT->GetListOfFunctions()->Delete();
      gROOT->GetListOfBrowsers()->Delete();
      gROOT->GetListOfGeometries()->Delete();
      //gROOT->GetListOfSpecials()->Delete();
      // Create some styles
      gStyle = 0;
      TStyle::BuildStyles();
      gROOT->SetStyle("Default");
      lnk = lnk->Next();
   }
   fpind << "</ul>" << endl;
   writeTrailer(html, fpind);
}

void MakeTutorials() {
   // Bring the ROOT tutorials on the web, see http://root.cern/root/html/tutorials/.
   // Demonstrates the use of THtml:Convert() in a realistic context.

   if (!gSystem->OpenDirectory("htmldoc")) gSystem->MakeDirectory("htmldoc");
   if (!gSystem->OpenDirectory("htmldoc/tutorials")) gSystem->MakeDirectory("htmldoc/tutorials");
   gEnv->SetValue("Unix.*.Root.Html.SourceDir", "$(ROOTSYS)");
   gEnv->SetValue("Root.Html.ViewCVS","http://root.cern/viewcvs/trunk/%f?view=log");
   gEnv->SetValue("Root.Html.Search", "http://www.google.com/search?q=%s+site%3A%u");
   THtml html;
   html.LoadAllLibs();
   //gROOT->ProcessLine(".x htmlLoadlibs.C");
   html.CreateAuxiliaryFiles();
   writeTutorials(html);
}
