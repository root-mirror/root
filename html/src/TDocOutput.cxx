// @(#)root/html:$Id: 7ff9b72609794c66acf6b369c4eeddfbfc63cf55 $
// Author: Axel Naumann 2007-01-09

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TDocOutput.h"

#include "Riostream.h"
#include "TClassDocOutput.h"
#include "TClassEdit.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TDocInfo.h"
#include "TDocParser.h"
#include "THtml.h"
#include "TInterpreter.h"
#include "TMethod.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TUrl.h"
#include "TVirtualMutex.h"
#include "TVirtualPad.h"
#include "TVirtualViewer3D.h"
#include <vector>
#include <list>
#include <set>
#include <sstream>
#include <stdlib.h>

namespace {

   typedef std::vector<std::string> Words_t;
   typedef Words_t::const_iterator SectionStart_t;

   class TSectionInfo {
   public:
      TSectionInfo(SectionStart_t start, size_t chars, size_t size):
         fStart(start), fChars(chars), fSize(size) {};

         SectionStart_t fStart;
         size_t fChars;
         size_t fSize;
   };
   typedef std::list<TSectionInfo> SectionStarts_t;

   static void Sections_BuildIndex(SectionStarts_t& sectionStarts,
      SectionStart_t begin, SectionStart_t end,
      size_t maxPerSection)
   {
      // for each assumed section border, check that previous entry's
      // char[selectionChar] differs, else move section start forward

      SectionStart_t cursor = begin;
      if (sectionStarts.empty() || sectionStarts.back().fStart != cursor)
         sectionStarts.push_back(TSectionInfo(cursor, 1, 0));

      SectionStarts_t::iterator prevSection = sectionStarts.end();
      --prevSection;

      while (cursor != end) {
         size_t numLeft = end - cursor;
         size_t assumedNumSections = (numLeft + maxPerSection - 1 ) / maxPerSection;
         size_t step = ((numLeft + assumedNumSections - 1) / assumedNumSections);
         if (!step || step >= numLeft) return;
         cursor += step;
         if (cursor == end) break;

         SectionStart_t addWhichOne = prevSection->fStart;

         size_t selectionChar=1;
         for (; selectionChar <= cursor->length() && addWhichOne == prevSection->fStart;
            ++selectionChar) {
            SectionStart_t checkPrev = cursor;
            while (--checkPrev != prevSection->fStart
               && !strncasecmp(checkPrev->c_str(), cursor->c_str(), selectionChar)) { }

            SectionStart_t checkNext = cursor;
            while (++checkNext != end
               && !strncasecmp(checkNext->c_str(), cursor->c_str(), selectionChar)) { }

            // if the previous matching one is closer but not previous section start, take it!
            if (checkPrev != prevSection->fStart) {
               if ((cursor - checkPrev) <= (checkNext - cursor))
                  addWhichOne = ++checkPrev;
               else if (checkNext != end
                  && (size_t)(checkNext - cursor) < maxPerSection) {
                  addWhichOne = checkNext;
               }
            }
         }
         if (addWhichOne == prevSection->fStart)
            addWhichOne = cursor;

         selectionChar = 1;
         while (selectionChar <= prevSection->fStart->length()
            && selectionChar <= addWhichOne->length()
            && !strncasecmp(prevSection->fStart->c_str(), addWhichOne->c_str(), selectionChar))
            ++selectionChar;

         sectionStarts.push_back(TSectionInfo(addWhichOne, selectionChar, 0));
         cursor = addWhichOne;
         ++prevSection;
      } // while cursor != end
   }

   static void Sections_SetSize(SectionStarts_t& sectionStarts, const Words_t &words)
   {
      // Update the length of the sections
      for (SectionStarts_t::iterator iSectionStart = sectionStarts.begin();
         iSectionStart != sectionStarts.end(); ++iSectionStart) {
         SectionStarts_t::iterator next = iSectionStart;
         ++next;
         if (next == sectionStarts.end()) {
            iSectionStart->fSize = (words.end() - iSectionStart->fStart);
            break;
         }
         iSectionStart->fSize = (next->fStart - iSectionStart->fStart);
      }
   }

   static void Sections_PostMerge(SectionStarts_t& sectionStarts, const size_t maxPerSection)
   {
      // Merge sections that ended up being too small, up to maxPerSection entries
      for (SectionStarts_t::iterator iSectionStart = sectionStarts.begin();
         iSectionStart != sectionStarts.end();) {
         SectionStarts_t::iterator iNextSectionStart = iSectionStart;
         ++iNextSectionStart;
         if (iNextSectionStart == sectionStarts.end()) break;
         if (iNextSectionStart->fSize + iSectionStart->fSize < maxPerSection) {
            iSectionStart->fSize += iNextSectionStart->fSize;
            sectionStarts.erase(iNextSectionStart);
         } else ++iSectionStart;
      }
   }

   static void GetIndexChars(const Words_t& words, UInt_t numSectionsIn,
      std::vector<std::string> &sectionMarkersOut)
   {
      // Given a list of words (class names, in this case), this function builds an
      // optimal set of about numSectionIn sections (even if almost all words start
      // with a "T"...), and returns the significant characters for each section start
      // in sectionMarkersOut.

      const size_t maxPerSection = (words.size() + numSectionsIn - 1)/ numSectionsIn;
      SectionStarts_t sectionStarts;
      Sections_BuildIndex(sectionStarts, words.begin(), words.end(), maxPerSection);
      Sections_SetSize(sectionStarts, words);
      Sections_PostMerge(sectionStarts, maxPerSection);

      // convert to index markers
      sectionMarkersOut.clear();
      sectionMarkersOut.resize(sectionStarts.size());
      size_t idx = 0;
      for (SectionStarts_t::iterator iSectionStart = sectionStarts.begin();
         iSectionStart != sectionStarts.end(); ++iSectionStart)
         sectionMarkersOut[idx++] =
            iSectionStart->fStart->substr(0, iSectionStart->fChars);
   }

   static void GetIndexChars(const std::list<std::string>& wordsIn, UInt_t numSectionsIn,
      std::vector<std::string> &sectionMarkersOut)
   {
      // initialize word vector
      Words_t words(wordsIn.size());
      size_t idx = 0;
      for (std::list<std::string>::const_iterator iWord = wordsIn.begin(); iWord != wordsIn.end(); ++iWord)
         words[idx++] = *iWord;
      GetIndexChars(words, numSectionsIn, sectionMarkersOut);
   }

}

extern "C" { // std::qsort on solaris wants the sorter to be extern "C"

   /////////////////////////////////////////////////////////////////////////////
   /// Friend function for sorting strings, case insensitive
   ///
   ///
   /// Input: name1 - pointer to the first string
   ///        name2 - pointer to the second string
   ///
   ///  NOTE: This function compares its arguments and returns an integer less
   ///        than, equal to, or greater than zero, depending on whether name1
   ///        is lexicographically less than, equal to, or greater than name2,
   ///        but characters are forced to lower-case prior to comparison.
   ///
   ///

   static int CaseInsensitiveSort(const void *name1, const void *name2)
   {
      return (strcasecmp(*((char **) name1), *((char **) name2)));
   }
}

namespace {

   // std::list::sort(with_stricmp_predicate) doesn't work with Solaris CC...
   static void sort_strlist_stricmp(std::vector<std::string>& l)
   {
      // sort strings ignoring case - easier for humans
      struct posList {
         const char* str;
         size_t pos;
      };
      posList* carr = new posList[l.size()];
      size_t idx = 0;
      for (size_t iS = 0, iSE = l.size(); iS < iSE; ++iS) {
         carr[idx].pos = iS;
         carr[idx++].str = l[iS].c_str();
      }
      qsort(&carr[0].str, idx, sizeof(posList), CaseInsensitiveSort);
      std::vector<std::string> lsort(l.size());
      for (size_t iS = 0, iSE = l.size(); iS < iSE; ++iS) {
         lsort[iS].swap(l[carr[iS].pos]);
      }
      delete [] carr;
      l.swap(lsort);
   }

}

//______________________________________________________________________________
//
// THtml generated documentation is written to file by TDocOutput. So far only
// output of HTML is implemented. Customization of the output should be done
// with THtml's interfaces - TDocOutput should not be used nor re-implemented
// directly.
//
// TDocOutput generates the index tables:
// * classes (THtml invokes TClassDocOutput for each),
// * inheritance hierarchy,
// * types and typedefs,
// * libraries,
// * the product index, and
// * the module index (including the links to per-module documentation).
// It invokes AT&T's GraphViz tool (dot) if available; charts benefit a lot
// from it.
//
// TDocOutput also writes all pages' header and footer, which can be customized
// by calling THtml::SetHeader(), THtml::SetFooter().
//______________________________________________________________________________

ClassImp(TDocOutput);

////////////////////////////////////////////////////////////////////////////////

TDocOutput::TDocOutput(THtml& html): fHtml(&html)
{}

////////////////////////////////////////////////////////////////////////////////

TDocOutput::~TDocOutput()
{}

////////////////////////////////////////////////////////////////////////////////
/// Add a link around str, with title comment.
/// Update str so it surrounds the link.

void TDocOutput::AddLink(TSubString& str, TString& link, const char* comment)
{
   // prepend "./" to allow callers to replace a different relative directory
   if (ReferenceIsRelative(link) && !link.BeginsWith("./"))
      link.Prepend("./");
   link.Prepend("<a href=\"");
   link += "\"";
   if (comment && strlen(comment)) {
      link += " title=\"";
      TString description(comment);
      ReplaceSpecialChars(description);
      description.ReplaceAll("\"", "&quot;");
      link += description;
      link += "\"";
   }
   link += ">";

   str.String().Insert(str.Start() + str.Length(), "</a>");
   str.String().Insert(str.Start(), link);

   TString &strString = str.String();
   TSubString update = strString(str.Start(), str.Length() + link.Length() + 4);
   str = update;
}

////////////////////////////////////////////////////////////////////////////////
/// adjust the path of links for source files, which are in src/, but need
/// to point to relpath (usually "../"). Simply replaces "=\"./" by "=\"../"

void TDocOutput::AdjustSourcePath(TString& line, const char* relpath /*= "../"*/)
{
   TString replWithRelPath("=\"@!@");
   line.ReplaceAll("=\"../", replWithRelPath + "../" + relpath);
   line.ReplaceAll("=\"./", replWithRelPath + relpath);
   line.ReplaceAll("=\"@!@","=\"");
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a text file into a html file.
/// outfilename doesn't have an extension yet; up to us to decide.
/// We generate HTML, so our extension is ".html".
/// See THtml::Convert() for the other parameters.

void TDocOutput::Convert(std::istream& in, const char* infilename,
                         const char* outfilename, const char *title,
                         const char *relpath /*= "../"*/, Int_t includeOutput /*=0*/,
                         const char* context /*= ""*/,
                         TGClient* gclient /*= 0*/)
{
   TString htmlFilename(outfilename);
   htmlFilename += ".html";

   std::ofstream out(htmlFilename);

   if (!out.good()) {
      Error("Convert", "Can't open file '%s' !", htmlFilename.Data());
      return;
   }

   // write a HTML header
   WriteHtmlHeader(out, title, relpath);

   if (context && context[0])
      out << context << std::endl;
   else if (title && title[0])
      out << "<h1 class=\"convert\">" << title << "</h1>" << std::endl;

   Int_t numReuseCanvases = 0;
   if (includeOutput && !(includeOutput & THtml::kForceOutput)) {
      void* dirHandle = gSystem->OpenDirectory(gSystem->GetDirName(htmlFilename));
      if (dirHandle) {
         FileStat_t infile_stat;
         if (!gSystem->GetPathInfo(infilename, infile_stat)) {
            // can stat.
            const char* outfile = 0;
            TString firstCanvasFileBase(gSystem->BaseName(outfilename));
            firstCanvasFileBase += "_0.png";
            // first check whether the firstCanvasFile exists:
            Bool_t haveFirstCanvasFile = false;
            while ((outfile = gSystem->GetDirEntry(dirHandle))) {
               if (firstCanvasFileBase == outfile) {
                  haveFirstCanvasFile = true;
                  break;
               }
            }
            gSystem->FreeDirectory(dirHandle);

            FileStat_t outfile_stat;
            TString firstCanvasFile = outfilename;
            firstCanvasFile += "_0.png";
            Int_t maxIdx = -1;
            if (haveFirstCanvasFile && !gSystem->GetPathInfo(firstCanvasFile, outfile_stat)
                && outfile_stat.fMtime > infile_stat.fMtime) {
               // the first canvas file exists and it is newer than the script, so we reuse
               // the canvas files. We need to know how many there are:
               dirHandle = gSystem->OpenDirectory(gSystem->GetDirName(htmlFilename));
               TString stem(gSystem->BaseName(outfilename));
               stem += "_";
               TString dir = gSystem->GetDirName(htmlFilename);
               while ((outfile = gSystem->GetDirEntry(dirHandle))) {
                  if (strncmp(outfile, stem, stem.Length()))
                     continue;
                  const char* posext = strrchr(outfile, '.');
                  if (!posext || strcmp(posext, ".png"))
                     continue;

                  // extract the mod time of the PNG file
                  if (gSystem->GetPathInfo(dir + "/" + outfile, outfile_stat))
                     // can't stat!
                     continue;

                  if (outfile_stat.fMtime > infile_stat.fMtime) {
                     ++numReuseCanvases;
                     // The canvas PNG is newer than the script, so
                     // extract the index of the canvas
                     TString idxStr(outfile + stem.Length());
                     idxStr.Remove(idxStr.Length() - 4);
                     Int_t idx = idxStr.Atoi();
                     if (maxIdx < idx)
                        maxIdx = idx;
                  }
               }
               gSystem->FreeDirectory(dirHandle);
               if (maxIdx + 1 != numReuseCanvases)
                  // bad: the number of canvases to reuse noes not correspond to the highest index we saw.
                  // we will need to regenerate everything.
                  numReuseCanvases = 0;
            }
         } // infile can be stat'ed
      } // can open output directory
   } // canvases wanted

   if (numReuseCanvases)
      Printf("Convert: %s (reusing %d saved canvas%s)", htmlFilename.Data(), numReuseCanvases, (numReuseCanvases > 1 ? "es" : ""));
   else
      Printf("Convert: %s", htmlFilename.Data());

   UInt_t nCanvases = numReuseCanvases;
   if (includeOutput) {
      if (!numReuseCanvases) {
         // need to run the script
         if (includeOutput & THtml::kSeparateProcessOutput) {
            TString baseInFileName = gSystem->BaseName(infilename);
            TPMERegexp reOutFile(baseInFileName + "_[[:digit:]]+\\.png");

            // remove all files matching what saveScriptOutput.C could produce:
            void* outdirH = gSystem->OpenDirectory(gSystem->GetDirName(outfilename));
            if (outdirH) {
               // the directory exists.
               const char* outdirE = 0;
               while ((outdirE = gSystem->GetDirEntry(outdirH))) {
                  if (reOutFile.Match(outdirE)) {
                     gSystem->Unlink(outdirE);
                  }
               }
               gSystem->FreeDirectory(outdirH);
            }

            gSystem->Exec(TString::Format("ROOT_HIST=0 root.exe -l -q %s $ROOTSYS/etc/html/saveScriptOutput.C\\(\\\"%s\\\",\\\"%s\\\",%d\\)",
                          gROOT->IsBatch() ? "-b" : "",
                          infilename,
                          gSystem->GetDirName(outfilename).Data(),
                          includeOutput & THtml::kCompiledOutput));

            // determine how many output files were created:
            outdirH = gSystem->OpenDirectory(gSystem->GetDirName(outfilename));
            if (outdirH) {
               // the directory exists.
               const char* outdirE = 0;
               while ((outdirE = gSystem->GetDirEntry(outdirH))) {
                  if (reOutFile.Match(outdirE)) {
                     ++nCanvases;
                  }
               }
               gSystem->FreeDirectory(outdirH);
            }
         } else {
            // run in this ROOT process
            TString pwd(gSystem->pwd());
            gSystem->cd(gSystem->GetDirName(infilename));

            TList* gClientGetListOfWindows = nullptr;
            TObject* gClientGetDefaultRoot = nullptr;
            std::set<TObject*> previousWindows;
            if (gclient) {
               gROOT->ProcessLine(TString::Format("*((TList**)0x%lx) = ((TGClient*)0x%lx)->GetListOfWindows();",
                                                  (ULong_t)&gClientGetListOfWindows, (ULong_t)gclient));
               gROOT->ProcessLine(TString::Format("*((TObject**)0x%lx) = ((TGClient*)0x%lx)->GetDefaultRoot();",
                                                  (ULong_t)&gClientGetDefaultRoot, (ULong_t)gclient));
               TObject *win = nullptr;
               TIter iWin(gClientGetListOfWindows);
               while((win = iWin())) {
                  TObject *winGetParent = nullptr;
                  gROOT->ProcessLine(TString::Format("*((TObject**)0x%lx) = ((TGWindow*)0x%lx)->GetParent();",
                                                     (ULong_t)&winGetParent, (ULong_t)win));
                  if (winGetParent == gClientGetDefaultRoot)
                     previousWindows.insert(win);
               }
            } else {
               if (gROOT->GetListOfCanvases()->GetSize())
                  previousWindows.insert(gROOT->GetListOfCanvases()->Last());
            }
            TIter iTimer(gSystem->GetListOfTimers());
            std::set<TObject*> timersBefore;
            TObject* timerOld = 0;
            while ((timerOld = iTimer()))
               timersBefore.insert(timerOld);

            TString cmd(".x ");
            cmd += gSystem->BaseName(infilename);
            if (includeOutput & THtml::kCompiledOutput)
               cmd += "+";
            gInterpreter->SaveContext();
            gInterpreter->SaveGlobalsContext();
            Int_t err;
            gROOT->ProcessLine(cmd, &err);
            gSystem->ProcessEvents();
            gSystem->cd(pwd);

            if (err == TInterpreter::kNoError) {
               if (gclient) {
                  TClass* clRootCanvas = TClass::GetClass("TRootCanvas");
                  TClass* clGMainFrame = TClass::GetClass("TGMainFrame");
                  TObject* win = 0;
                  TIter iWin(gClientGetListOfWindows);
                  while((win = iWin())) {
                     TObject* winGetParent = 0;
                     gROOT->ProcessLine(TString::Format("*((TObject**)0x%lx) = ((TGWindow*)0x%lx)->GetParent();",
                                                        (ULong_t)&winGetParent, (ULong_t)win));
                     Bool_t winIsMapped = kFALSE;
                     if (winGetParent == gClientGetDefaultRoot)
                        gROOT->ProcessLine(TString::Format("*((Bool_t*)0x%lx) = ((TGWindow*)0x%lx)->IsMapped();",
                                                           (ULong_t)&winIsMapped, (ULong_t)win));
                     if (winIsMapped && previousWindows.find(win) == previousWindows.end()
                         && win->InheritsFrom(clGMainFrame)) {
                        gROOT->ProcessLine(TString::Format("((TGWindow*)0x%lx)->MapRaised();", (ULong_t)win));
                        Bool_t isRootCanvas = win->InheritsFrom(clRootCanvas);
                        Bool_t hasEditor = false;
                        if (isRootCanvas) {
                           gROOT->ProcessLine(TString::Format("*((Bool_t*)0x%lx) = ((TRootCanvas*)0x%lx)->HasEditor();",
                                                              (ULong_t)&hasEditor, (ULong_t)win));
                        }
                        if (isRootCanvas && !hasEditor) {
                           TVirtualPad* pad = 0;
                           gROOT->ProcessLine(TString::Format("*((TVirtualPad**)0x%lx) = ((TRootCanvas*)0x%lx)->Canvas();",
                                                              (ULong_t)&pad, (ULong_t)win));
                           if (!pad->HasViewer3D() || pad->GetViewer3D()->InheritsFrom("TViewer3DPad")) {
                              pad->SaveAs(TString::Format("%s_%d.png", outfilename, nCanvases++));
                           }
                        } else
                           gROOT->ProcessLine(TString::Format("((TGWindow*)0x%lx)->SaveAs(\"%s_%d.png\");",
                                                              (ULong_t)win, outfilename, nCanvases++));
                     }
                  }
               } else {
                  // no gClient
                  TVirtualPad* pad = 0;
                  TVirtualPad* last = 0;
                  if (!previousWindows.empty())
                     last = (TVirtualPad*) *previousWindows.begin();
                  TIter iCanvas(gROOT->GetListOfCanvases());
                  while ((pad = (TVirtualPad*) iCanvas())) {
                     if (last) {
                        if (last == pad) last = 0;
                        continue;
                     }
                     pad->SaveAs(TString::Format("%s_%d.png", outfilename, nCanvases++));
                  }
               }
               gInterpreter->Reset();
               gInterpreter->ResetGlobals();
               TIter iTimerRemove(gSystem->GetListOfTimers());
               TTimer* timer = 0;
               while ((timer = (TTimer*) iTimerRemove()))
                  if (timersBefore.find(timer) == timersBefore.end())
                     gSystem->RemoveTimer(timer);
            }
         } // run script in this ROOT process
      }
      out << "<table><tr><td style=\"vertical-align:top;padding-right:2em;\">" << std::endl;
   }
   out << "<div class=\"listing\"><pre class=\"listing\">" << std::endl;

   TDocParser parser(*this);
   parser.Convert(out, in, relpath, (includeOutput) /* determines whether it's code or not */,
                  kFALSE /*interpretDirectives*/);

   out << "</pre></div>" << std::endl;

   WriteLineNumbers(out, parser.GetLineNumber(), gSystem->BaseName(infilename));

   if (includeOutput) {
      out << "</td><td style=\"vertical-align:top;\">" << std::endl;
      out << "<table>" << std::endl;
      for (UInt_t i = 0; i < nCanvases; ++i) {
         TString pngname = TString::Format("%s_%d.png", gSystem->BaseName(outfilename), i);
         out << "<tr><td><a href=\"" << pngname << "\">" << std::endl
             << "<img src=\"" << pngname << "\" id=\"canv" << i << "\" alt=\"thumb\" style=\"border:none;width:22em;\" "
            "onmouseover=\"javascript:canv" << i << ".style.width='auto';\" />" << std::endl
             << "</a></td></tr>" << std::endl;
         }
      out << "</table>" << std::endl;
      out << "</td></tr></table>" << std::endl;
   }

   // write a HTML footer
   WriteHtmlFooter(out, relpath);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy file to HTML directory
///
///
///  Input: sourceName - source file name (fully qualified i.e. file system path)
///         destName   - optional destination name, if not
///                      specified it would be the same
///                      as the source file name
///
/// Output: TRUE if file is successfully copied, or
///         FALSE if it's not
///
///
///   NOTE: The destination directory is always fHtml->GetOutputDir()
///

Bool_t TDocOutput::CopyHtmlFile(const char *sourceName, const char *destName)
{
   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   TString sourceFile(sourceName);

   if (!sourceFile.Length()) {
      Error("Copy", "Can't copy file '%s' to '%s' directory - source file name invalid!", sourceName,
            fHtml->GetOutputDir().Data());
      return kFALSE;
   }

   // destination file name
   TString destFile;
   if (!destName || !*destName)
      destFile = gSystem->BaseName(sourceFile);
   else
      destFile = gSystem->BaseName(destName);

   gSystem->PrependPathName(fHtml->GetOutputDir(), destFile);

   // Get info about a file
   Long64_t size;
   Long_t id, flags, sModtime, dModtime;
   sModtime = 0;
   dModtime = 0;
   if (gSystem->GetPathInfo(sourceFile, &id, &size, &flags, &sModtime)
      || gSystem->GetPathInfo(destFile, &id, &size, &flags, &dModtime)
      || sModtime > dModtime)
      if (gSystem->CopyFile(sourceFile, destFile, kTRUE) < 0) {
         Error("Copy", "Can't copy file '%s' to '%s'!",
               sourceFile.Data(), destFile.Data());
         return kFALSE;
      }

   return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a hierarchical class list
/// The algorithm descends from the base classes and branches into
/// all derived classes. Mixing classes are displayed several times.
///
///

void TDocOutput::CreateHierarchy()
{
   // if (CreateHierarchyDot()) return;

   TString filename("ClassHierarchy.html");
   gSystem->PrependPathName(fHtml->GetOutputDir(), filename);

   // open out file
   std::ofstream out(filename);

   if (!out.good()) {
      Error("CreateHierarchy", "Can't open file '%s' !", filename.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), filename.Data());

   // write out header
   WriteHtmlHeader(out, "Class Hierarchy");

   WriteTopLinks(out, 0);

   out << "<h1>Class Hierarchy</h1>" << std::endl;


   // loop on all classes
   TClassDocInfo* cdi = 0;
   TIter iClass(fHtml->GetListOfClasses());
   while ((cdi = (TClassDocInfo*)iClass())) {
      if (!cdi->HaveSource())
         continue;

      // get class
      TDictionary *dictPtr = cdi->GetClass();
      TClass *basePtr = dynamic_cast<TClass*>(dictPtr);
      if (basePtr == 0) {
         if (!dictPtr)
            Warning("THtml::CreateHierarchy", "skipping class %s\n", cdi->GetName());
         continue;
      }

      TClassDocOutput cdo(*fHtml, basePtr, 0);
      cdo.CreateClassHierarchy(out, cdi->GetHtmlFileName());
   }

   // write out footer
   WriteHtmlFooter(out);
}

////////////////////////////////////////////////////////////////////////////////
/// Create index of all classes
///

void TDocOutput::CreateClassIndex()
{
   // create CSS file, we need it
   fHtml->CreateAuxiliaryFiles();

   TString filename("ClassIndex.html");
   gSystem->PrependPathName(fHtml->GetOutputDir(), filename);

   // open indexFile file
   std::ofstream indexFile(filename.Data());

   if (!indexFile.good()) {
      Error("CreateClassIndex", "Can't open file '%s' !", filename.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), filename.Data());

   // write indexFile header
   WriteHtmlHeader(indexFile, "Class Index");

   WriteTopLinks(indexFile, 0);

   indexFile << "<h1>Class Index</h1>" << std::endl;

   WriteModuleLinks(indexFile);

   std::vector<std::string> indexChars;
   if (fHtml->GetListOfClasses()->GetSize() > 10) {
      std::vector<std::string> classNames;
      {
         TIter iClass(fHtml->GetListOfClasses());
         TClassDocInfo* cdi = 0;
         while ((cdi = (TClassDocInfo*)iClass()))
            if (cdi->IsSelected() && cdi->HaveSource())
               classNames.push_back(cdi->GetName());
      }

      if (classNames.size() > 10) {
         indexFile << "<div id=\"indxShortX\"><h4>Jump to</h4>" << std::endl;
         // find index chars
         GetIndexChars(classNames, 50 /*sections*/, indexChars);
         for (UInt_t iIdxEntry = 0; iIdxEntry < indexChars.size(); ++iIdxEntry) {
            indexFile << "<a href=\"#idx" << iIdxEntry << "\">";
            ReplaceSpecialChars(indexFile, indexChars[iIdxEntry].c_str());
            indexFile << "</a>" << std::endl;
         }
         indexFile << "</div><br />" << std::endl;
      }
   }

   indexFile << "<ul id=\"indx\">" << std::endl;

   // loop on all classes
   UInt_t currentIndexEntry = 0;
   TIter iClass(fHtml->GetListOfClasses());
   TClassDocInfo* cdi = 0;
   Int_t i = 0;
   while ((cdi = (TClassDocInfo*)iClass())) {
      if (!cdi->IsSelected() || !cdi->HaveSource())
         continue;

      // get class
      TDictionary *currentDict = cdi->GetClass();
      TClass* currentClass = dynamic_cast<TClass*>(currentDict);
      if (!currentClass) {
         if (!currentDict)
            Warning("THtml::CreateClassIndex", "skipping class %s\n", cdi->GetName());
         continue;
      }

      indexFile << "<li class=\"idxl" << (i++)%2 << "\">";
      if (currentIndexEntry < indexChars.size()
         && !strncmp(indexChars[currentIndexEntry].c_str(), cdi->GetName(),
                     indexChars[currentIndexEntry].length()))
         indexFile << "<a name=\"idx" << currentIndexEntry++ << "\"></a>";

      TString htmlFile(cdi->GetHtmlFileName());
      if (htmlFile.Length()) {
         indexFile << "<a href=\"";
         indexFile << htmlFile;
         indexFile << "\"><span class=\"typename\">";
         ReplaceSpecialChars(indexFile, cdi->GetName());
         indexFile << "</span></a> ";
      } else {
         indexFile << "<span class=\"typename\">";
         ReplaceSpecialChars(indexFile, cdi->GetName());
         indexFile << "</span> ";
      }

      // write title == short doc
      ReplaceSpecialChars(indexFile, currentClass->GetTitle());
      indexFile << "</li>" << std::endl;
   }

   indexFile << "</ul>" << std::endl;

   // write indexFile footer
   WriteHtmlFooter(indexFile);
}


////////////////////////////////////////////////////////////////////////////////
/// Create the class index for each module, picking up documentation from the
/// module's TModuleDocInfo::GetInputPath() plus the (possibly relative)
/// THtml::GetModuleDocPath(). Also creates the library dependency plot if dot
/// exists, see THtml::HaveDot().

void TDocOutput::CreateModuleIndex()
{
   const char* title = "LibraryDependencies";
   TString dotfilename(title);
   gSystem->PrependPathName(fHtml->GetOutputDir(), dotfilename);

   std::ofstream libDepDotFile(dotfilename + ".dot");
   libDepDotFile << "digraph G {" << std::endl
                 << "ratio=compress;" << std::endl
                 << "node [fontsize=22,labeldistance=0.1];" << std::endl
                 << "edge [len=0.01];" << std::endl
                 << "fontsize=22;" << std::endl
                 << "size=\"16,16\";" << std::endl
                 << "overlap=false;" << std::endl
                 << "splines=true;" << std::endl
                 << "K=0.1;" << std::endl;

   TModuleDocInfo* module = 0;
   TIter iterModule(fHtml->GetListOfModules());

   std::stringstream sstrCluster;
   std::stringstream sstrDeps;
   while ((module = (TModuleDocInfo*)iterModule())) {
      if (!module->IsSelected())
         continue;

      std::vector<std::string> indexChars;
      TString filename(module->GetName());
      filename.ToUpper();
      filename.ReplaceAll("/","_");
      filename += "_Index.html";
      gSystem->PrependPathName(fHtml->GetOutputDir(), filename);
      std::ofstream outputFile(filename.Data());
      if (!outputFile.good()) {
         Error("CreateModuleIndex", "Can't open file '%s' !", filename.Data());
         continue;
      }
      Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), filename.Data());

      TString htmltitle("Index of ");
      TString moduletitle(module->GetName());
      moduletitle.ToUpper();
      htmltitle += moduletitle;
      WriteHtmlHeader(outputFile, htmltitle);

      WriteTopLinks(outputFile, module);

      outputFile << "<h2>" << htmltitle << "</h2>" << std::endl;

      // Module doc
      if (GetHtml()->GetModuleDocPath().Length()) {
         TString outdir(module->GetName());
         gSystem->PrependPathName(GetHtml()->GetOutputDir(), outdir);

         TString moduleDocDir;
         GetHtml()->GetPathDefinition().GetDocDir(module->GetName(), moduleDocDir);
         ProcessDocInDir(outputFile, moduleDocDir, outdir, module->GetName());
      }

      WriteModuleLinks(outputFile, module);

      std::list<std::string> classNames;
      {
         TIter iClass(module->GetClasses());
         TClassDocInfo* cdi = 0;
         while ((cdi = (TClassDocInfo*) iClass())) {
            if (!cdi->IsSelected() || !cdi->HaveSource())
               continue;
            classNames.push_back(cdi->GetName());

            if (classNames.size() > 1) continue;

            TClass* cdiClass = dynamic_cast<TClass*>(cdi->GetClass());
            if (!cdiClass)
               continue;

            TString libs(cdiClass->GetSharedLibs());
            Ssiz_t posDepLibs = libs.Index(' ');
            TString thisLib(libs);
            if (posDepLibs != kNPOS)
               thisLib.Remove(posDepLibs, thisLib.Length());

            {
               Ssiz_t posExt = thisLib.First('.');
               if (posExt != kNPOS)
                  thisLib.Remove(posExt, thisLib.Length());
            }

            if (!thisLib.Length())
               continue;

            // allocate entry, even if no dependencies
            TLibraryDocInfo *libdeps =
               (TLibraryDocInfo*)fHtml->GetLibraryDependencies()->FindObject(thisLib);
            if (!libdeps) {
               libdeps = new TLibraryDocInfo(thisLib);
               fHtml->GetLibraryDependencies()->Add(libdeps);
            }
            libdeps->AddModule(module->GetName());
            if (posDepLibs != kNPOS) {
               std::string lib;
               for(Ssiz_t pos = posDepLibs + 1; libs[pos]; ++pos) {
                  if (libs[pos] == ' ') {
                     if (thisLib.Length() && lib.length()) {
                        size_t posExt = lib.find('.');
                        if (posExt != std::string::npos)
                           lib.erase(posExt);
                        libdeps->AddDependency(lib);
                     }
                     lib.erase();
                  } else
                     lib += libs[pos];
               }
               if (lib.length() && thisLib.Length()) {
                  size_t posExt = lib.find('.');
                  if (posExt != std::string::npos)
                     lib.erase(posExt);
                  libdeps->AddDependency(lib);
               }
            } // if dependencies
         } // while next class in module
      } // just a scope block

      TIter iClass(module->GetClasses());
      TClassDocInfo* cdi = 0;
      UInt_t count = 0;
      UInt_t currentIndexEntry = 0;
      while ((cdi = (TClassDocInfo*) iClass())) {
         if (!cdi->IsSelected() || !cdi->HaveSource())
            continue;

         TDictionary *classPtr = cdi->GetClass();
         if (!classPtr) {
            Error("CreateModuleIndex", "Unknown class '%s' !", cdi->GetName());
            continue;
         }

         if (!count) {
            outputFile << "<h2>Class Index</h2>" << std::endl;

            if (classNames.size() > 10) {
               outputFile << "<div id=\"indxShortX\"><h4>Jump to</h4>" << std::endl;
               UInt_t numSections = classNames.size() / 10;
               if (numSections < 10) numSections = 10;
               if (numSections > 50) numSections = 50;
               // find index chars
               GetIndexChars(classNames, numSections, indexChars);
               for (UInt_t iIdxEntry = 0; iIdxEntry < indexChars.size(); ++iIdxEntry) {
                  outputFile << "<a href=\"#idx" << iIdxEntry << "\">";
                  ReplaceSpecialChars(outputFile, indexChars[iIdxEntry].c_str());
                  outputFile << "</a>" << std::endl;
               }
               outputFile << "</div><br />" << std::endl;
            }

            outputFile << "<ul id=\"indx\">" << std::endl;
         }

         // write a classname to an index file
         outputFile << "<li class=\"idxl" << (count++)%2 << "\">";
         if (currentIndexEntry < indexChars.size()
            && !strncmp(indexChars[currentIndexEntry].c_str(), cdi->GetName(),
                        indexChars[currentIndexEntry].length()))
            outputFile << "<a name=\"idx" << currentIndexEntry++ << "\"></a>";

         TString htmlFile(cdi->GetHtmlFileName());
         if (htmlFile.Length()) {
            outputFile << "<a href=\"";
            outputFile << htmlFile;
            outputFile << "\"><span class=\"typename\">";
            ReplaceSpecialChars(outputFile, classPtr->GetName());
            outputFile << "</span></a> ";
         } else {
            outputFile << "<span class=\"typename\">";
            ReplaceSpecialChars(outputFile, classPtr->GetName());
            outputFile << "</span> ";
         }

         // write title
         ReplaceSpecialChars(outputFile, classPtr->GetTitle());
         outputFile << "</li>" << std::endl;
      }


      if (count)
         outputFile << "</ul>" << std::endl;

      // write outputFile footer
      WriteHtmlFooter(outputFile);
   } // while next module

   // libCint is missing as we don't have class doc for it
   // We need it for dependencies nevertheless, so add it by hand.
   /*
   sstrCluster << "subgraph clusterlibCint {" << std::endl
      << "style=filled;" << std::endl
      << "color=lightgray;" << std::endl
      << "label=\"libCint\";" << std::endl
      << "\"CINT\" [style=filled,color=white,fontsize=10]" << std::endl
      << "}" << std::endl;
   */

   // simplify the library dependencies, by removing direct links
   // that are equivalent to indirect ones, e.g. instead of having both
   // A->C, A->B->C, keep only A->B->C.

   TIter iLib(fHtml->GetLibraryDependencies());
   TLibraryDocInfo* libinfo = 0;
   while ((libinfo = (TLibraryDocInfo*)iLib())) {
      if (!libinfo->GetName() || !libinfo->GetName()[0]) continue;

      std::set<std::string>& deps = libinfo->GetDependencies();
      for (std::set<std::string>::iterator iDep = deps.begin();
           iDep != deps.end(); ) {
         Bool_t already_indirect = kFALSE;
         for (std::set<std::string>::const_iterator iDep2 = deps.begin();
              !already_indirect && iDep2 != deps.end(); ++iDep2) {
            if (iDep == iDep2) continue;
            TLibraryDocInfo* libinfo2 = (TLibraryDocInfo*)
               fHtml->GetLibraryDependencies()->FindObject(iDep2->c_str());
            if (!libinfo2) continue;
            const std::set<std::string>& deps2 = libinfo2->GetDependencies();
            already_indirect |= deps2.find(*iDep) != deps2.end();
         }
         if (already_indirect) {
            std::set<std::string>::iterator iRemove = iDep;
            // Advance the iterator before erasing the element which invalidates the iterator.
            ++iDep;
            deps.erase(iRemove);
         } else {
            ++iDep;
         }
      } // for library dependencies of module in library
   } // for libraries

   iLib.Reset();
   while ((libinfo = (TLibraryDocInfo*)iLib())) {
      if (!libinfo->GetName() || !libinfo->GetName()[0]) continue;

      const std::set<std::string>& modules = libinfo->GetModules();
      if (modules.size() > 1) {
         sstrCluster << "subgraph cluster" << libinfo->GetName() << " {" << std::endl
                     << "style=filled;" << std::endl
                     << "color=lightgray;" << std::endl
                     << "label=\"";
         if (!strcmp(libinfo->GetName(), "libCore"))
            sstrCluster << "Everything depends on ";
         sstrCluster << libinfo->GetName() << "\";" << std::endl;

         for (std::set<std::string>::const_iterator iModule = modules.begin();
              iModule != modules.end(); ++iModule) {
            TString modURL(*iModule);
            modURL.ReplaceAll("/", "_");
            modURL.ToUpper();
            sstrCluster << "\"" << *iModule << "\" [style=filled,color=white,URL=\""
                        << modURL << "_Index.html\"];" << std::endl;
         }
         sstrCluster << std::endl
                     << "}" << std::endl;
      } else {
         // only one module
         TString modURL(*modules.begin());
         modURL.ReplaceAll("/", "_");
         modURL.ToUpper();
         sstrCluster << "\"" << *modules.begin()
                     << "\" [label=\"" << libinfo->GetName()
                     << "\",style=filled,color=lightgray,shape=box,URL=\""
                     << modURL << "_Index.html\"];" << std::endl;
      }

      // GetSharedLib doesn't mention libCore or libCint; add them by hand
      /*
        if (iLibDep->first != "libCore")
        sstrDeps << "\"" << iModule->first << "\" -> \"BASE\" [lhead=clusterlibCore];" << std::endl;
        sstrDeps << "\"" << iModule->first << "\" -> \"CINT\" [lhead=clusterlibCint];" << std::endl;
      */

      const std::string& mod = *(modules.begin());
      const std::set<std::string>& deps = libinfo->GetDependencies();
      for (std::set<std::string>::const_iterator iDep = deps.begin();
            iDep != deps.end(); ++iDep) {
         // cannot create dependency on iDep directly, use its first module instead.
         TLibraryDocInfo* depLibInfo = (TLibraryDocInfo*)
            fHtml->GetLibraryDependencies()->FindObject(iDep->c_str());
         if (!depLibInfo || depLibInfo->GetModules().empty())
            continue; // ouch!

         const std::string& moddep = *(depLibInfo->GetModules().begin());
         sstrDeps << "\"" << mod << "\" -> \"" << moddep << "\";" << std::endl;
      }
      // make sure libCore ends up at the bottom
      sstrDeps << "\"" << mod <<  "\" -> \"CONT\" [style=invis];" << std::endl;
   } // for libs

   libDepDotFile << sstrCluster.str() << std::endl
      << sstrDeps.str();
   libDepDotFile << "}" << std::endl;
   libDepDotFile.close();

   std::ofstream out(dotfilename + ".html");
   if (!out.good()) {
      Error("CreateModuleIndex", "Can't open file '%s.html' !",
            dotfilename.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), (dotfilename + ".html").Data());
   // write out header
   WriteHtmlHeader(out, "Library Dependencies");

   WriteTopLinks(out, 0);

   out << "<h1>Library Dependencies</h1>" << std::endl;

   RunDot(dotfilename, &out, kFdp);

   out << "<img alt=\"Library Dependencies\" class=\"classcharts\" usemap=\"#Map" << title << "\" src=\"" << title << ".png\"/>" << std::endl;

   // write out footer
   WriteHtmlFooter(out);
}

////////////////////////////////////////////////////////////////////////////////
/// Fetch documentation from THtml::GetProductDocDir() and put it into the
/// product index page.

void TDocOutput::CreateProductIndex()
{
   //TString outFile(GetHtml()->GetProductName());
   //outFile += ".html";
   TString outFile("index.html");
   gSystem->PrependPathName(GetHtml()->GetOutputDir(), outFile);
   std::ofstream out(outFile);

   if (!out.good()) {
      Error("CreateProductIndex", "Can't open file '%s' !", outFile.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", "", outFile.Data());

   WriteHtmlHeader(out, GetHtml()->GetProductName() + " Reference Guide");

   WriteTopLinks(out, 0);

   out << "<h1>" << GetHtml()->GetProductName() + " Reference Guide</h1>" << std::endl;

   TString prodDoc;
   if (GetHtml()->GetPathDefinition().GetDocDir("", prodDoc))
      ProcessDocInDir(out, prodDoc, GetHtml()->GetOutputDir(), "./");

   WriteModuleLinks(out);

   out << "<h2>Chapters</h2>" << std::endl
      << "<h3><a href=\"./ClassIndex.html\">Class Index</a></h3>" << std::endl
      << "<p>A complete list of all classes defined in " << GetHtml()->GetProductName() << "</p>" << std::endl
      << "<h3><a href=\"./ClassHierarchy.html\">Class Hierarchy</a></h3>" << std::endl
      << "<p>A hierarchy graph of all classes, showing each class's base and derived classes</p>" << std::endl
      << "<h3><a href=\"./ListOfTypes.html\">Type Index</a></h3>" << std::endl
      << "<p>A complete list of all types</p>" << std::endl
      << "<h3><a href=\"./LibraryDependencies.html\">Library Dependency</a></h3>" << std::endl
      << "<p>A diagram showing all of " << GetHtml()->GetProductName() << "'s libraries and their dependencies</p>" << std::endl;

   WriteHtmlFooter(out);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a forwarding page for each typedef pointing to a class.

void TDocOutput::CreateClassTypeDefs()
{
   TDocParser parser(*this);

   TIter iClass(GetHtml()->GetListOfClasses());
   TClassDocInfo* cdi = 0;
   while ((cdi = (TClassDocInfo*) iClass())) {
      if (cdi->GetListOfTypedefs().IsEmpty())
         continue;
      TIter iTypedefs(&cdi->GetListOfTypedefs());
      TDataType* dt = 0;
      while ((dt = (TDataType*) iTypedefs())) {
         if (gDebug > 0)
            Info("CreateClassTypeDefs", "Creating typedef %s to class %s",
                 dt->GetName(), cdi->GetName());
         // create a filename
         TString filename(dt->GetName());
         NameSpace2FileName(filename);

         gSystem->PrependPathName(fHtml->GetOutputDir(), filename);

         filename += ".html";

         // open class file
         std::ofstream outfile(filename);

         if (!outfile.good()) {
            Error("CreateClassTypeDefs", "Can't open file '%s' !", filename.Data());
            continue;
         }

         WriteHtmlHeader(outfile, dt->GetName());

         outfile << "<a name=\"TopOfPage\"></a>" << std::endl;

         TString dtName(dt->GetName());
         ReplaceSpecialChars(dtName);
         TString sTitle("typedef ");
         sTitle += dtName;

         TClass* cls = dynamic_cast<TClass*>(cdi->GetClass());
         if (cls) {
            // show box with lib, include
            // needs to go first to allow title on the left
            TString sInclude;
            TString sLib;
            const char* lib=cls->GetSharedLibs();
            GetHtml()->GetPathDefinition().GetIncludeAs(cls, sInclude);
            if (lib) {
               char* libDup=StrDup(lib);
               char* libDupSpace=strchr(libDup,' ');
               if (libDupSpace) *libDupSpace = 0;
               char* libDupEnd=libDup+strlen(libDup);
               while (libDupEnd!=libDup)
                  if (*(--libDupEnd)=='.') {
                     *libDupEnd=0;
                     break;
                  }
               sLib = libDup;
               delete[] libDup;
            }
            outfile << "<script type=\"text/javascript\">WriteFollowPageBox('"
                    << sTitle << "','" << sLib << "','" << sInclude << "');</script>" << std::endl;
         }

         TString modulename;
         fHtml->GetModuleNameForClass(modulename, cls);
         TModuleDocInfo* module = (TModuleDocInfo*) fHtml->GetListOfModules()->FindObject(modulename);
         WriteTopLinks(outfile, module, dt->GetName());

         outfile << "<div class=\"dropshadow\"><div class=\"withshadow\">";
         outfile << "<h1>" << sTitle << "</h1>" << std::endl
            << "<div class=\"classdescr\">" << std::endl;

         outfile << dtName << " is a typedef to ";
         std::string shortClsName(fHtml->ShortType(cdi->GetName()));
         parser.DecorateKeywords(outfile, shortClsName.c_str());
         outfile << std::endl
            << "</div>" << std::endl
            << "</div></div><div style=\"clear:both;\"></div>" << std::endl;

         // the typedef isn't a data member, but the CSS is applicable nevertheless
         outfile << std::endl << "<div id=\"datamembers\">" << std::endl
            << "<table class=\"data\" cellspacing=\"0\">" << std::endl;
         outfile << "<tr class=\"data";
         outfile << "\"><td class=\"datatype\">typedef ";
         parser.DecorateKeywords(outfile, dt->GetFullTypeName());
         outfile << "</td><td class=\"dataname\">";
         ReplaceSpecialChars(outfile, dt->GetName());
         if (dt->GetTitle() && dt->GetTitle()[0]) {
            outfile << "</td><td class=\"datadesc\">";
            ReplaceSpecialChars(outfile, dt->GetTitle());
         } else outfile << "</td><td>";
         outfile << "</td></tr>" << std::endl
            << "</table></div>" << std::endl;

         // write footer
         WriteHtmlFooter(outfile);

      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create index of all data types

void TDocOutput::CreateTypeIndex()
{
   // open file
   TString outFile("ListOfTypes.html");
   gSystem->PrependPathName(fHtml->GetOutputDir(), outFile);
   std::ofstream typesList(outFile);

   if (!typesList.good()) {
      Error("CreateTypeIndex", "Can't open file '%s' !", outFile.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", "", outFile.Data());

   // write typesList header
   WriteHtmlHeader(typesList, "List of data types");
   typesList << "<h2> List of data types </h2>" << std::endl;

   typesList << "<dl><dd>" << std::endl;

   // make loop on data types
   std::vector<std::string> typeNames(gROOT->GetListOfTypes()->GetSize());

   {
      TDataType *type;
      TIter nextType(gROOT->GetListOfTypes());
      size_t tnIdx = 0;

      while ((type = (TDataType *) nextType()))
         // no templates ('<' and '>'), no idea why the '(' is in here...
         if (*type->GetTitle() && !strchr(type->GetName(), '(')
             && !( strchr(type->GetName(), '<') && strchr(type->GetName(),'>'))
             && type->GetName())
            typeNames[tnIdx++] = type->GetName();
      typeNames.resize(tnIdx);
   }

   sort_strlist_stricmp(typeNames);

   std::vector<std::string> indexChars;
   if (typeNames.size() > 10) {
      typesList << "<div id=\"indxShortX\"><h4>Jump to</h4>" << std::endl;
      // find index chars
      GetIndexChars(typeNames, 10 /*sections*/, indexChars);
      for (UInt_t iIdxEntry = 0; iIdxEntry < indexChars.size(); ++iIdxEntry) {
         typesList << "<a href=\"#idx" << iIdxEntry << "\">";
         ReplaceSpecialChars(typesList, indexChars[iIdxEntry].c_str());
         typesList << "</a>" << std::endl;
      }
      typesList << "</div><br />" << std::endl;
   }

   typesList << "<ul id=\"indx\">" << std::endl;

   int idx = 0;
   UInt_t currentIndexEntry = 0;

   for (std::vector<std::string>::iterator iTypeName = typeNames.begin();
      iTypeName != typeNames.end(); ++iTypeName) {
      TDataType* type = gROOT->GetType(iTypeName->c_str(), kFALSE);
      typesList << "<li class=\"idxl" << idx%2 << "\">";
      if (currentIndexEntry < indexChars.size()
         && !strncmp(indexChars[currentIndexEntry].c_str(), iTypeName->c_str(),
                     indexChars[currentIndexEntry].length()))
         typesList << "<a name=\"idx" << currentIndexEntry++ << "\"></a>" << std::endl;
      typesList << "<a name=\"";
      ReplaceSpecialChars(typesList, iTypeName->c_str());
      typesList << "\"><span class=\"typename\">";
      ReplaceSpecialChars(typesList, iTypeName->c_str());
      typesList << "</span></a> ";
      ReplaceSpecialChars(typesList, type->GetTitle());
      typesList << "</li>" << std::endl;
      ++idx;
   }
   typesList << "</ul>" << std::endl;

   // write typesList footer
   WriteHtmlFooter(typesList);

   // close file
   typesList.close();

}


////////////////////////////////////////////////////////////////////////////////
/// Add some colors etc to a source entity, contained in str.
/// The type of what's contained in str is given by type.
/// It's called e.g. by TDocParser::BeautifyLine().
/// This function should assume that only str.Begin() is valid.
/// When inserting into str.String(), str.Begin() must be updated.

void TDocOutput::DecorateEntityBegin(TString& str, Ssiz_t& pos, TDocParser::EParseContext type)
{
   Ssiz_t originalLen = str.Length();

   switch (type) {
      case TDocParser::kCode: break;
      case TDocParser::kComment:
         str.Insert(pos, "<span class=\"comment\">");
         break;
      case TDocParser::kDirective:
         break;
      case TDocParser::kString:
         str.Insert(pos, "<span class=\"string\">");
         break;
      case TDocParser::kKeyword:
         str.Insert(pos, "<span class=\"keyword\">");
         break;
      case TDocParser::kCPP:
         str.Insert(pos, "<span class=\"cpp\">");
         break;
      case TDocParser::kVerbatim:
         str.Insert(pos, "<pre>");
         break;
      default:
         Error("DecorateEntityBegin", "Unhandled / invalid entity type %d!", (Int_t)type);
         return;
   }

   Ssiz_t addedLen = str.Length() - originalLen;
   pos += addedLen;
}

////////////////////////////////////////////////////////////////////////////////
/// Add some colors etc to a source entity, contained in str.
/// The type of what's contained in str is given by type.
/// It's called e.g. by TDocParser::BeautifyLine().
/// This function should assume that only str."End()"
/// (i.e. str.Begin()+str.Length()) is valid.
/// When inserting into str.String(), str.Length() must be updated.

void TDocOutput::DecorateEntityEnd(TString& str, Ssiz_t& pos, TDocParser::EParseContext type)
{
   Ssiz_t originalLen = str.Length();

   switch (type) {
      case TDocParser::kCode: break;
      case TDocParser::kComment:
         str.Insert(pos, "</span>");
         break;
      case TDocParser::kDirective:
         break;
      case TDocParser::kString:
         str.Insert(pos, "</span>");
         break;
      case TDocParser::kKeyword:
         str.Insert(pos, "</span>");
         break;
      case TDocParser::kCPP:
         str.Insert(pos, "</span>");
         break;
      case TDocParser::kVerbatim:
         str.Insert(pos, "</pre>");
         break;
      default:
         Error("DecorateEntityBegin", "Unhandled / invalid entity type %d!", (Int_t)type);
         return;
   }
   Ssiz_t addedLen = str.Length() - originalLen;
   pos += addedLen;
}

////////////////////////////////////////////////////////////////////////////////
/// Special author treatment; called when TDocParser::fSourceInfo[kInfoAuthor] is set.
/// Modifies the author(s) description, which is a comma separated list of tokens
/// either in the format
/// (i) "FirstName LastName " or
/// (ii) "FirstName LastName <link> more stuff"
/// The first one generates an XWho link (CERN compatible),
/// the second a http link (WORLD compatible), <link> being e.g.
/// <mailto:user@host.bla> or <http://www.host.bla/page>.

void TDocOutput::FixupAuthorSourceInfo(TString& authors)
{
   TString original(authors);
   authors = "";

   TString author;
   Ssiz_t pos = 0;
   Bool_t firstAuthor = kTRUE;
   while (original.Tokenize(author, pos, ",")) {
      author.Strip(TString::kBoth);

      if (!firstAuthor)
         authors += ", ";
      firstAuthor = kFALSE;

      // do we have a link for the current name?
      Ssiz_t cLink = author.First('<'); // look for link start tag
      if (cLink != kNPOS) {
         // split NAME <LINK> POST
         // into  <a href="LINK">NAME</a> POST
         Ssiz_t endLink = author.Index(">", cLink + 1);
         if(endLink == kNPOS)
            endLink = author.Length();
         authors += "<a href=\"";
         authors += author(cLink + 1,  endLink - (cLink + 1));
         authors += "\">";
         authors += author(0, cLink);
         authors += "</a>";
         if (endLink != author.Length())
            authors += author(endLink + 1, author.Length());
      } else {
         authors += "<a href=\"";
         authors += fHtml->GetXwho();

         // separate Firstname Middlename Lastname by '+'
         TString namePart;
         Ssiz_t posNamePart = 0;
         Bool_t firstNamePart = kTRUE;
         while (author.Tokenize(namePart, posNamePart, " ")) {
            namePart.Strip(TString::kBoth);
            if (!namePart.Length())
               continue;
            if (isdigit(namePart[0])) continue; //likely a date
            if (!firstNamePart)
               authors += '+';
            firstNamePart = kFALSE;
            authors += namePart;
         }
         authors += "\">";
         authors += author;
         authors += "</a>";
      }
   } // while next author
}

////////////////////////////////////////////////////////////////////////////////
/// Check if file is modified
///
///
///  Input: classPtr - pointer to the class
///         type     - file type to compare with
///                    values: kSource, kInclude, kTree
///
/// Output: TRUE     - if file is modified since last time
///         FALSE    - if file is up to date
///

Bool_t TDocOutput::IsModified(TClass * classPtr, EFileType type)
{
   TString sourceFile;
   TString classname(classPtr->GetName());
   TString filename;
   TString dir;

   switch (type) {
   case kSource:
      {
         TString declFile;
         if (classPtr->GetImplFileLine()) {
            fHtml->GetImplFileName(classPtr, kTRUE, sourceFile);
         }
         fHtml->GetDeclFileName(classPtr, kTRUE, declFile);
         Long64_t size;
         Long_t id, flags, iModtime, dModtime;
         if (!(gSystem->GetPathInfo(sourceFile, &id, &size, &flags, &iModtime))) {
            if (!(gSystem->GetPathInfo(declFile, &id, &size, &flags, &dModtime))) {
               if (iModtime < dModtime) {
                  // decl is newer than impl
                  sourceFile = declFile;
               }
            }
         }
         dir = "src";
         gSystem->PrependPathName(fHtml->GetOutputDir(), dir);
         filename = classname;
         NameSpace2FileName(filename);
         gSystem->PrependPathName(dir, filename);
         if (classPtr->GetImplFileLine())
            filename += ".cxx.html";
         else
            filename += ".h.html";
         break;
      }

   case kInclude:
      fHtml->GetDeclFileName(classPtr, kFALSE, filename);
      filename = gSystem->BaseName(filename);
      fHtml->GetDeclFileName(classPtr, kTRUE, sourceFile);
      gSystem->PrependPathName(fHtml->GetOutputDir(), filename);
      break;

   case kTree:
      fHtml->GetDeclFileName(classPtr, kTRUE, sourceFile);
      NameSpace2FileName(classname);
      gSystem->PrependPathName(fHtml->GetOutputDir(), classname);
      filename = classname;
      filename += "_Tree.pdf";
      break;

   case kDoc:
      {
         TString declFile;
         if (classPtr->GetImplFileLine()) {
            fHtml->GetImplFileName(classPtr, kTRUE, sourceFile);
         }
         fHtml->GetDeclFileName(classPtr, kTRUE, declFile);
         Long64_t size;
         Long_t id, flags, iModtime, dModtime;
         if (!(gSystem->GetPathInfo(sourceFile, &id, &size, &flags, &iModtime))) {
            if (!(gSystem->GetPathInfo(declFile, &id, &size, &flags, &dModtime))) {
               if (iModtime < dModtime) {
                  // decl is newer than impl
                  sourceFile = declFile;
               }
            }
         }
         filename = classname;
         NameSpace2FileName(filename);
         gSystem->PrependPathName(fHtml->GetOutputDir(), filename);
         filename += ".html";
         break;
      }

   default:
      Error("IsModified", "Unknown file type !");
   }

   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   // Get info about a file
   Long64_t size;
   Long_t id, flags, sModtime, dModtime;

   if (!(gSystem->GetPathInfo(sourceFile, &id, &size, &flags, &sModtime))) {
      if (!(gSystem->GetPathInfo(filename, &id, &size, &flags, &dModtime))) {
         return (sModtime > dModtime);
      }
   }

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace "::" in name by "__"
/// Replace "<", ">", " ", ",", "~", "=" in name by "_"
/// Replace "A::X<A::Y>" by "A::X<-p0Y>",
///         "A::B::X<A::B::Y>" by "A::B::X<-p1Y>", etc

void TDocOutput::NameSpace2FileName(TString& name)
{
   TString encScope(name);
   Ssiz_t posTemplate = encScope.Index('<');
   if (posTemplate != kNPOS) {
      // strip default template params
      name = fHtml->ShortType(name);
      TString templateArgs = encScope(posTemplate, encScope.Length());
      encScope.Remove(posTemplate, encScope.Length());
      // shorten the name a bit:
      // convert A::B::X<A::B::Y> to A::X<-p1Y>, i.e.
      // the filename A__X_A__Y_ to A__X_-p1Y_
      // The rule: if the enclosing scope up to the N-th scope matches,
      // the name becomes -pN
      Ssiz_t posName = encScope.Last(':');
      if (posName != kNPOS) {
         Int_t numDblColumn = encScope.CountChar(':');
         while (numDblColumn > 1) {
            encScope.Remove(posName + 1, encScope.Length());
            numDblColumn -= 2;
            templateArgs.ReplaceAll(encScope, TString::Format("-p%d", numDblColumn / 2));
            encScope.Remove(encScope.Length() - 2, 2);
            posName = encScope.Last(':');
            if (posName == kNPOS)
               break; // should be handled by numDblColumn...
         }
         name.Replace(posTemplate, name.Length(), templateArgs);
      }
   }

   if (name.Length() > 240) { // really 240! It might get some extra prefix or extension
      // 8.3 is dead, but e.g. ext2 can only hold 255 chars in a file name.
      // So mangle name to "beginning_of_name"-h"hash"."extension", where
      // beginning_of_name is short enough such that the full name is <255 characters.

      TString hash;
      TDocParser::AnchorFromLine(name, hash);
      hash.Prepend("-h");
      Ssiz_t posDot = name.Last('.');
      TString ext;
      if (posDot != kNPOS)
         ext = name(posDot, name.Length());
      Ssiz_t namelen = 240 - hash.Length() - ext.Length();
      name = name(0, namelen) + hash + ext;
   }

   const char* replaceWhat = ":<> ,~=";
   for (Ssiz_t i=0; i < name.Length(); ++i)
      if (strchr(replaceWhat, name[i]))
         name[i] = '_';
}

////////////////////////////////////////////////////////////////////////////////
/// Write links to files indir/*.txt, indir/*.html (non-recursive) to out.
/// If one of the files is called "index.{html,txt}" it will be
/// included in out (instead of copying it to outdir and generating a link
/// to linkdir). txt files are passed through Convert().
/// The files' links are sorted alphabetically.

void TDocOutput::ProcessDocInDir(std::ostream& out, const char* indir,
                                 const char* outdir, const char* linkdir)
{
   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   void * dirHandle = gSystem->OpenDirectory(indir);
   if (!dirHandle) return;

   const char* entry = 0;
   std::list<std::string> files;
   while ((entry = gSystem->GetDirEntry(dirHandle))) {
      FileStat_t stat;
      TString filename(entry);
      gSystem->PrependPathName(indir, filename);
      if (gSystem->GetPathInfo(filename, stat)) // funny ret
         continue;
      if (!R_ISREG(stat.fMode)) continue;

      if (TString(entry).BeginsWith("index.", TString::kIgnoreCase)) {
         // This is the part we put directly (verbatim) into the module index.
         // If it ends on ".txt" we run Convert first.
         if (filename.EndsWith(".txt", TString::kIgnoreCase)) {
            std::ifstream in(filename);
            if (in) {
               out << "<pre>"; // this is what e.g. the html directive expects
               TDocParser parser(*this);
               parser.Convert(out, in, "./", kFALSE /* no code */, kTRUE /*process Directives*/);
               out << "</pre>";
            }
         } else if (filename.EndsWith(".html", TString::kIgnoreCase)) {
            std::ifstream in(filename);
            TString line;
            while (in) {
               if (!line.ReadLine(in)) break;
               out << line << std::endl;
            }
         } else
            files.push_back(filename.Data());
      } else
         files.push_back(filename.Data());
   }

   std::stringstream furtherReading;
   files.sort();
   for (std::list<std::string>::const_iterator iFile = files.begin();
      iFile != files.end(); ++iFile) {
      TString filename(iFile->c_str());
      if (gSystem->AccessPathName(outdir))
         if (gSystem->mkdir(outdir, kTRUE) == -1)
            // bad - but let's still try to create the output
            Error("CreateModuleIndex", "Cannot create output directory %s", outdir);

      TString outfile(gSystem->BaseName(filename));
      gSystem->PrependPathName(outdir, outfile);

      if (!filename.EndsWith(".txt", TString::kIgnoreCase)
          && !filename.EndsWith(".html", TString::kIgnoreCase)) {
         // copy to outdir, who know whether it's needed...
         if (gSystem->CopyFile(filename, outfile, kTRUE) == -1) {
            Error("CreateModuleIndex", "Cannot copy file %s to %s",
                  filename.Data(), outfile.Data());
            continue;
         }
         continue;
      }

      // Just copy and link this page.
      if (outfile.EndsWith(".txt", TString::kIgnoreCase)) {
         // convert first
         outfile.Remove(outfile.Length()-3, 3);
         outfile += "html";
         std::ifstream inFurther(filename);
         std::ofstream outFurther(outfile);
         if (inFurther && outFurther) {
            outFurther << "<pre>"; // this is what e.g. the html directive expects
            TDocParser parser(*this);
            parser.Convert(outFurther, inFurther, "../", kFALSE /*no code*/, kTRUE /*process Directives*/);
            outFurther << "</pre>";
         }
      } else {
         if (gSystem->CopyFile(filename, outfile, kTRUE) == -1)
            continue;
      }
      TString showname(gSystem->BaseName(outfile));
      furtherReading << "<a class=\"linkeddoc\" href=\"" << linkdir << "/" << showname << "\">";
      showname.Remove(showname.Length() - 5, 5); // .html
      showname.ReplaceAll("_", " ");
      ReplaceSpecialChars(furtherReading, showname);
      furtherReading << "</a> " << std::endl;
   }

   gSystem->FreeDirectory(dirHandle);
   if (furtherReading.str().length())
      out << "<h3>Further Reading</h3><div id=\"furtherreading\">" << std::endl
          << furtherReading.str() << "</div><h3>List of Classes</h3>" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a reference to a class documentation page.
/// str encloses the text to create the reference for (e.g. name of instance).
/// comment will be added e.g. as tooltip text.
/// After the reference is put into str.String(), str will enclose the reference
/// and the original text. Example:
/// Input:
///  str.String(): "a gHtml test"
///  str.Begin():  2
///  str.Length(): 5
/// Output:
///  str.String(): "a <a href="THtml.html">gHtml</a> test"
///  str.Begin():  2
///  str.Length(): 30

void TDocOutput::ReferenceEntity(TSubString& str, TClass* entity, const char* comment /*= 0*/)
{
   TString link;
   fHtml->GetHtmlFileName(entity, link);

   if (comment && !strcmp(comment, entity->GetName()))
      comment = "";

   AddLink(str, link, comment);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a reference to a data member documentation page.
/// str encloses the text to create the reference for (e.g. name of instance).
/// comment will be added e.g. as tooltip text.
/// After the reference is put into str.String(), str will enclose the reference
/// and the original text. Example:
/// Input:
///  str.String(): "a gHtml test"
///  str.Begin():  2
///  str.Length(): 5
/// Output:
///  str.String(): "a <a href="THtml.html">gHtml</a> test"
///  str.Begin():  2
///  str.Length(): 30

void TDocOutput::ReferenceEntity(TSubString& str, TDataMember* entity, const char* comment /*= 0*/)
{
   TString link;
   TClass* scope = entity->GetClass();
   fHtml->GetHtmlFileName(scope, link);
   link += "#";

   TString mangledName;
   if (scope) {
      mangledName = scope->GetName();
      NameSpace2FileName(mangledName);
      link += mangledName;
      link += ":";
   }

   mangledName = entity->GetName();
   NameSpace2FileName(mangledName);
   link += mangledName;

   TString description;
   if (!comment) {
      description = entity->GetFullTypeName();
      description += " ";
      if (scope) {
         description += scope->GetName();
         description += "::";
      }
      description += entity->GetName();
      comment = description.Data();
   }

   if (comment && !strcmp(comment, entity->GetName()))
      comment = "";

   AddLink(str, link, comment);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a reference to a type documentation page.
/// str encloses the text to create the reference for (e.g. name of instance).
/// comment will be added e.g. as tooltip text.
/// After the reference is put into str.String(), str will enclose the reference
/// and the original text. Example:
/// Input:
///  str.String(): "a gHtml test"
///  str.Begin():  2
///  str.Length(): 5
/// Output:
///  str.String(): "a <a href="THtml.html">gHtml</a> test"
///  str.Begin():  2
///  str.Length(): 30

void TDocOutput::ReferenceEntity(TSubString& str, TDataType* entity, const char* comment /*= 0*/)
{
   TString mangledEntity(entity->GetName());
   NameSpace2FileName(mangledEntity);

   TString link;
   TClassDocInfo* cdi = 0;
   bool isClassTypedef = entity->GetType() == -1;
   if (isClassTypedef)
      /* is class/ struct / union */
      isClassTypedef = isClassTypedef && (entity->Property() & 7);
   if (isClassTypedef) {
      std::string shortTypeName(fHtml->ShortType(entity->GetFullTypeName()));
      cdi = (TClassDocInfo*) GetHtml()->GetListOfClasses()->FindObject(shortTypeName.c_str());
   }
   if (cdi) {
      link = mangledEntity + ".html";
   } else {
      link = "ListOfTypes.html#";
      link += mangledEntity;
   }

   if (comment && !strcmp(comment, entity->GetName()))
      comment = "";

   AddLink(str, link, comment);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a reference to a method documentation page.
/// str encloses the text to create the reference for (e.g. name of instance).
/// comment will be added e.g. as tooltip text.
/// After the reference is put into str.String(), str will enclose the reference
/// and the original text. Example:
/// Input:
///  str.String(): "a gHtml test"
///  str.Begin():  2
///  str.Length(): 5
/// Output:
///  str.String(): "a <a href="THtml.html">gHtml</a> test"
///  str.Begin():  2
///  str.Length(): 30

void TDocOutput::ReferenceEntity(TSubString& str, TMethod* entity, const char* comment /*= 0*/)
{
   TString link;
   TClass* scope = entity->GetClass();
   fHtml->GetHtmlFileName(scope, link);
   link += "#";

   TString mangledName(scope->GetName());
   NameSpace2FileName(mangledName);
   link += mangledName;
   link += ":";

   mangledName = entity->GetName();
   NameSpace2FileName(mangledName);
   link += mangledName;

   TString description;
   if (!comment && entity->GetClass()) {
      TIter iMeth(scope->GetListOfMethods());
      TMethod* mCand = 0;
      while ((mCand = (TMethod*)iMeth()))
         if (!strcmp(mCand->GetName(), entity->GetName())) {
            if (description.Length()) {
               description += " or overloads";
               break;
            }
            description = mCand->GetPrototype();
         }
      comment = description.Data();
   }

   if (comment && !strcmp(comment, entity->GetName()))
      comment = "";

   AddLink(str, link, comment);
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether reference is a relative reference, and can (or should)
/// be prependen by relative paths. For HTML, check that it doesn't start
/// with "http://" or "https://"

Bool_t TDocOutput::ReferenceIsRelative(const char* reference) const
{
   return !reference ||
      strncmp(reference, "http", 4) ||
      (strncmp(reference + 4, "://", 3) && strncmp(reference + 4, "s://", 4));
}

////////////////////////////////////////////////////////////////////////////////
/// Replace ampersand, less-than and greater-than character, writing to out.
/// If 0 is returned, no replacement needs to be done.

const char* TDocOutput::ReplaceSpecialChars(char c)
{
   /*
   if (fEscFlag) {
      fEscFlag = kFALSE;
      return buf;
   } else if (c == fEsc) {
      // text.Remove(pos, 1); - NO! we want to keep it nevertheless!
      fEscFlag = kTRUE;
      return buf;
   }

   */
   switch (c) {
      case '<': return "&lt;";
      case '&': return "&amp;";
      case '>': return "&gt;";
   };
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace ampersand, less-than and greater-than character
///
///
/// Input: text - text where replacement will happen,
///        pos  - index of char to be replaced; will point to next char to be
///               replaced when function returns
///

void TDocOutput::ReplaceSpecialChars(TString& text, Ssiz_t &pos)
{
   const char c = text[pos];
   const char* replaced = ReplaceSpecialChars(c);
   if (replaced) {
         text.Replace(pos, 1, replaced);
         pos += strlen(replaced) - 1;
   }
   ++pos;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace ampersand, less-than and greater-than character
///
///
/// Input: text - text where replacement will happen,
///

void TDocOutput::ReplaceSpecialChars(TString& text) {
   Ssiz_t pos = 0;
   while (pos < text.Length())
         ReplaceSpecialChars(text, pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Replace ampersand, less-than and greater-than characters, writing to out
///
///
/// Input: out    - output file stream
///        string - pointer to an array of characters
///

void TDocOutput::ReplaceSpecialChars(std::ostream& out, const char *string)
{
   while (string && *string) {
      const char* replaced = ReplaceSpecialChars(*string);
      if (replaced)
         out << replaced;
      else
         out << *string;
      string++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Run filename".dot", creating filename".png", and - if outMap is !=0,
/// filename".map", which gets then included literally into outMap.

Bool_t TDocOutput::RunDot(const char* filename, std::ostream* outMap /* =0 */,
                          EGraphvizTool gvwhat /*= kDot*/) {
   if (!fHtml->HaveDot())
      return kFALSE;

   TString runDot;
   switch (gvwhat) {
   case kNeato: runDot = "neato"; break;
   case kFdp: runDot = "fdp"; break;
   case kCirco: runDot = "circo"; break;
   default: runDot = "dot";
   };
   if (fHtml->GetDotDir() && *fHtml->GetDotDir())
      gSystem->PrependPathName(fHtml->GetDotDir(), runDot);
   runDot += " -q1 -Tpng -o";
   runDot += filename;
   runDot += ".png ";
   if (outMap) {
      runDot += "-Tcmap -o";
      runDot += filename;
      runDot += ".map ";
   }
   runDot += filename;
   runDot += ".dot";

   if (gDebug > 3)
      Info("RunDot", "Running: %s", runDot.Data());
   Int_t retDot = gSystem->Exec(runDot);
   if (gDebug < 4 && !retDot)
      gSystem->Unlink(Form("%s.dot", filename));

   if (!retDot && outMap) {
      std::ifstream inmap(Form("%s.map", filename));
      std::string line;
      std::getline(inmap, line);
      if (inmap && !inmap.eof()) {
         *outMap << "<map name=\"Map" << gSystem->BaseName(filename)
            << "\" id=\"Map" << gSystem->BaseName(filename) << "\">" << std::endl;
         while (inmap && !inmap.eof()) {
            if (line.compare(0, 6, "<area ") == 0) {
               size_t posEndTag = line.find('>');
               if (posEndTag != std::string::npos)
                  line.replace(posEndTag, 1, "/>");
            }
            *outMap << line << std::endl;
            std::getline(inmap, line);
         }
         *outMap << "</map>" << std::endl;
      }
      inmap.close();
      if (gDebug < 7)
         gSystem->Unlink(Form("%s.map", filename));
   }

   if (retDot) {
      Error("RunDot", "Error running %s!", runDot.Data());
      fHtml->SetFoundDot(kFALSE);
      return kFALSE;
   }

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Write HTML header
///
/// Internal method invoked by the overload

void TDocOutput::WriteHtmlHeader(std::ostream& out, const char *titleNoSpecial,
                            const char* dir /*=""*/, TClass *cls /*=0*/,
                            const char* header)
{
   std::ifstream addHeaderFile(header);

   if (!addHeaderFile.good()) {
      Warning("THtml::WriteHtmlHeader",
              "Can't open html header file %s\n", header);
      return;
   }

   TString declFileName;
   if (cls) fHtml->GetDeclFileName(cls, kFALSE, declFileName);
   TString implFileName;
   if (cls) fHtml->GetImplFileName(cls, kFALSE, implFileName);

   const TString& charset = GetHtml()->GetCharset();
   TDatime date;
   TString strDate(date.AsString());
   TString line;

   while (!addHeaderFile.eof()) {

      line.ReadLine(addHeaderFile, kFALSE);
      if (addHeaderFile.eof())
         break;

      if (line) {

         if (!cls && (
            line.Index("%CLASS%") != kNPOS ||
            line.Index("%INCFILE%") != kNPOS ||
            line.Index("%SRCFILE%") != kNPOS))
            continue; // skip class line for non-class files

         TString txt(line);

         txt.ReplaceAll("%TITLE%", titleNoSpecial);
         txt.ReplaceAll("%DATE%", strDate);
         txt.ReplaceAll("%RELDIR%", dir);
         txt.ReplaceAll("%CHARSET%", charset);

         if (cls) {
            txt.ReplaceAll("%CLASS%", cls->GetName());
            txt.ReplaceAll("%INCFILE%", declFileName);
            txt.ReplaceAll("%SRCFILE%", implFileName);
         }

         out << txt << std::endl;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write HTML header
///
///
/// Input: out   - output file stream
///        title - title for the HTML page
///        cls   - current class
///        dir   - relative directory to reach the top
///                ("" for html doc, "../" for src/*cxx.html etc)
///
/// evaluates the Root.Html.Header setting:
/// * if not set, the standard header is written. (ROOT)
/// * if set, and ends with a "+", the standard header is written and this file
///   included afterwards. (ROOT, USER)
/// * if set but doesn't end on "+" the file specified will be written instead
///   of the standard header (USER)
///
/// Any occurrence of "%TITLE%" (without the quotation marks) in the user
/// provided header file will be replaced by the value of this method's
/// parameter "title" before written to the output file. %CLASS% is replaced by
/// the class name, %INCFILE% by the header file name as given by
/// TClass::GetDeclFileName() and %SRCFILE% by the source file name as given by
/// TClass::GetImplFileName(). If the header is written for a non-class page,
/// i.e. cls==0, lines containing %CLASS%, %INCFILE%, or %SRCFILE% will be
/// skipped.

void TDocOutput::WriteHtmlHeader(std::ostream& out, const char *title,
                            const char* dir /*=""*/, TClass *cls/*=0*/)
{
   TString userHeader = GetHtml()->GetHeader();
   TString noSpecialCharTitle(title);
   ReplaceSpecialChars(noSpecialCharTitle);

   Ssiz_t lenUserHeader = userHeader.Length();
   // standard header output if Root.Html.Header is not set, or it's set and it ends with a "+".
   Bool_t bothHeaders = lenUserHeader > 0 && userHeader[lenUserHeader - 1] == '+';
   if (lenUserHeader == 0 || bothHeaders) {
      TString header("header.html");
      gSystem->PrependPathName(fHtml->GetEtcDir(), header);
      WriteHtmlHeader(out, noSpecialCharTitle, dir, cls, header);
   }

   if (lenUserHeader != 0) {
      if (bothHeaders)
         userHeader.Remove(lenUserHeader - 1);
      WriteHtmlHeader(out, noSpecialCharTitle, dir, cls, userHeader);
   };
}

////////////////////////////////////////////////////////////////////////////////
/// Write HTML footer
///
/// Internal method invoked by the overload

void TDocOutput::WriteHtmlFooter(std::ostream& out, const char* /*dir*/,
                                 const char* lastUpdate, const char* author,
                                 const char* copyright, const char* footer)
{
   static const char* templateSITags[TDocParser::kNumSourceInfos] = { "%UPDATE%", "%AUTHOR%", "%COPYRIGHT%", "%CHANGED%", "%GENERATED%"};

   TString today;
   TDatime dtToday;
   today.Form("%d-%02d-%02d %02d:%02d", dtToday.GetYear(), dtToday.GetMonth(), dtToday.GetDay(), dtToday.GetHour(), dtToday.GetMinute());

   TString datimeString;
   if (!lastUpdate || !lastUpdate[0]) {
      lastUpdate = today;
   }
   const char* siValues[TDocParser::kNumSourceInfos] = { lastUpdate, author, copyright, lastUpdate, today };

   std::ifstream addFooterFile(footer);

   if (!addFooterFile.good()) {
      Warning("THtml::WriteHtmlFooter",
              "Can't open html footer file %s\n", footer);
      return;
   }

   TString line;
   while (!addFooterFile.eof()) {

      line.ReadLine(addFooterFile, kFALSE);
      if (addFooterFile.eof())
         break;

      if (!line)
         continue;

      for (Int_t siTag = 0; siTag < (Int_t) TDocParser::kNumSourceInfos; ++siTag) {
         Ssiz_t siPos = line.Index(templateSITags[siTag]);
         if (siPos != kNPOS) {
            if (siValues[siTag] && siValues[siTag][0])
               line.Replace(siPos, strlen(templateSITags[siTag]), siValues[siTag]);
            else
               line = ""; // skip e.g. %AUTHOR% lines if no author is set
         }
      }

      out << line << std::endl;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Write HTML footer
///
///
/// Input: out        - output file stream
///        dir        - usually equal to "" or "../", depends of
///                     current file directory position, i.e. if
///                     file is in the fHtml->GetOutputDir(), then dir will be ""
///        lastUpdate - last update string
///        author     - author's name
///        copyright  - copyright note
///
/// Allows optional user provided footer to be written. Root.Html.Footer holds
/// the file name for this footer. For details see THtml::WriteHtmlHeader (here,
/// the "+" means the user's footer is written in front of Root's!) Occurrences
/// of %AUTHOR%, %CHANGED%, %GENERATED%, and %COPYRIGHT% in the user's file are replaced by
/// their corresponding values (author, lastUpdate, today, and copyright) before
/// written to out.
/// If no author is set (author == "", e.g. for ClassIndex.html") skip the whole
/// line of the footer template containing %AUTHOR%. Accordingly for %COPYRIGHT%.

void TDocOutput::WriteHtmlFooter(std::ostream& out, const char *dir,
                            const char *lastUpdate, const char *author,
                            const char *copyright)
{
   out << std::endl;

   TString userFooter = GetHtml()->GetFooter();

   if (userFooter.Length() != 0) {
      TString footer(userFooter);
      if (footer.EndsWith("+"))
         footer.Remove(footer.Length() - 1);
      WriteHtmlFooter(out, dir, lastUpdate, author, copyright, footer);
   };

   if (userFooter.Length() == 0 || userFooter.EndsWith("+")) {
      TString footer("footer.html");
      gSystem->PrependPathName(fHtml->GetEtcDir(), footer);
      WriteHtmlFooter(out, dir, lastUpdate, author, copyright, footer);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a div containing links to all topmost modules

void TDocOutput::WriteModuleLinks(std::ostream& out)
{
   if (fHtml->GetListOfModules()->GetSize()) {
      out << "<div id=\"indxModules\"><h4>Modules</h4>" << std::endl;
      // find index chars
      fHtml->SortListOfModules();
      TIter iModule(fHtml->GetListOfModules());
      TModuleDocInfo* module = 0;
      while ((module = (TModuleDocInfo*) iModule())) {
         if (!module->GetName() || strchr(module->GetName(), '/'))
            continue;
         if (module->IsSelected()) {
            TString name(module->GetName());
            name.ToUpper();
            out << "<a href=\"" << name << "_Index.html\">"
                << name << "</a>" << std::endl;
         }
      }
      out<< "</div><br />" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a div containing the line numbers (for a source listing) 1 to nLines.
/// Create links to the source file's line number and anchors, such that one can
/// jump to SourceFile.cxx.html#27 (using the anchor), and one can copy and paste
/// the link into e.g. gdb to get the text "SourceFile.cxx:27".

void TDocOutput::WriteLineNumbers(std::ostream& out, Long_t nLines, const TString& infileBase) const
{
   out << "<div id=\"linenums\">";
   for (Long_t i = 0; i < nLines; ++i) {
      // &nbsp; to force correct line height
      out << "<div class=\"ln\">&nbsp;<span class=\"lnfile\">" << infileBase
          << ":</span><a name=\"" << i + 1 << "\" href=\"#" << i + 1
          << "\" class=\"ln\">" << i + 1 << "</a></div>";
   }
   out << "</div>" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
/// Create a div containing links to all modules

void TDocOutput::WriteModuleLinks(std::ostream& out, TModuleDocInfo* super)
{
   if (super->GetSub().GetSize()) {
      TString superName(super->GetName());
      superName.ToUpper();
      out << "<div id=\"indxModules\"><h4>" << superName << " Modules</h4>" << std::endl;
      // find index chars
      super->GetSub().Sort();
      TIter iModule(&super->GetSub());
      TModuleDocInfo* module = 0;
      while ((module = (TModuleDocInfo*) iModule())) {
         if (module->IsSelected()) {
            TString name(module->GetName());
            name.ToUpper();
            TString link(name);
            link.ReplaceAll("/", "_");
            Ssiz_t posSlash = name.Last('/');
            if (posSlash != kNPOS)
               name.Remove(0, posSlash + 1);
            out << "<a href=\"" << link << "_Index.html\">" << name << "</a>" << std::endl;
         }
      }
      out<< "</div><br />" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a search link or a search box, based on THtml::GetSearchStemURL()
/// and THtml::GetSearchEngine(). The first one is preferred.

void TDocOutput::WriteSearch(std::ostream& out)
{
   // e.g. searchCmd = "http://www.google.com/search?q=%s+site%3A%u+-site%3A%u%2Fsrc%2F+-site%3A%u%2Fexamples%2F";
   const TString& searchCmd = GetHtml()->GetSearchStemURL();
   const TString& searchEngine = GetHtml()->GetSearchEngine();

   if (!searchCmd.Length() && !searchEngine.Length())
      return;

   if (searchCmd.Length()) {
      TUrl url(searchCmd);
      TString serverName(url.GetHost());
      if (serverName.Length()) {
         serverName.Prepend(" title=\"");
         serverName += "\" ";
      }
      // create search input
      out << "<script type=\"text/javascript\">" << std::endl
          << "function onSearch() {" << std::endl
          << "var s='" << searchCmd <<"';" << std::endl
          << "var ref=String(document.location.href).replace(/https?:\\/\\//,'').replace(/\\/[^\\/]*$/,'').replace(/\\//g,'%2F');" << std::endl
          << "window.location.href=s.replace(/%u/ig,ref).replace(/%s/ig,escape(document.searchform.t.value));" << std::endl
          << "return false;}" << std::endl
          << "</script>" << std::endl
          << "<form id=\"searchform\" name=\"searchform\" onsubmit=\"return onSearch()\" action=\"javascript:onSearch();\" method=\"post\">" << std::endl
          << "<input name=\"t\" size=\"30\" value=\"Search documentation...\" onfocus=\"if (document.searchform.t.value=='Search documentation...') document.searchform.t.value='';\"></input>" << std::endl
          << "<a id=\"searchlink\" " << serverName << " href=\"javascript:onSearch();\" onclick=\"return onSearch()\">Search</a></form>" << std::endl;
   } else if (searchEngine.Length())
      // create link to search engine page
      out << "<a class=\"descrheadentry\" href=\"" << searchEngine
          << "\">Search the Class Reference Guide</a>" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// make a link to the description

void TDocOutput::WriteLocation(std::ostream& out, TModuleDocInfo* module, const char* classname)
{
   out << "<div class=\"location\">" << std::endl; // location
   const char *productName = fHtml->GetProductName();
   out << "<a class=\"locationlevel\" href=\"index.html\">" << productName << "</a>" << std::endl;

   if (module) {
      TString modulename(module->GetName());
      modulename.ToUpper();
      TString modulePart;
      TString modulePath;
      Ssiz_t pos = 0;
      while (modulename.Tokenize(modulePart, pos, "/")) {
         if (pos == kNPOS && !classname)
            // we are documenting the module itself, no need to link it:
            break;
         if (modulePath.Length()) modulePath += "_";
         modulePath += modulePart;
         out << " &#187; <a class=\"locationlevel\" href=\"./" << modulePath << "_Index.html\">" << modulePart << "</a>" << std::endl;
      }
   }

   TString entityName;
   if (classname) entityName = classname;
   else if (module) {
      entityName = module->GetName();
      Ssiz_t posSlash = entityName.Last('/');
      if (posSlash != kNPOS)
         entityName.Remove(0, posSlash + 1);
      entityName.ToUpper();
   }
   if (entityName.Length()) {
      out << " &#187; <a class=\"locationlevel\" href=\"#TopOfPage\">";
      ReplaceSpecialChars(out, entityName);
      out << "</a>" << std::endl;
   }
   out << "</div>" << std::endl; // location
}


////////////////////////////////////////////////////////////////////////////////
/// Write the first part of the links shown ontop of each doc page;
/// one \<div\> has to be closed by caller so additional items can still
/// be added.

void TDocOutput::WriteTopLinks(std::ostream& out, TModuleDocInfo* module, const char* classname,
                               Bool_t withLocation)
{
   out << "<div id=\"toplinks\">" << std::endl;

   out << "<div class=\"descrhead\"><div class=\"descrheadcontent\">" << std::endl // descrhead line 1
      << "<span class=\"descrtitle\">Quick Links:</span>" << std::endl;

   // link to the user home page (if exist)
   const char* userHomePage = GetHtml()->GetHomepage();
   const char* productName = fHtml->GetProductName();
   if (!productName) {
      productName = "";
   } else if (!strcmp(productName, "ROOT")) {
      userHomePage = "";
   }
   if (userHomePage && *userHomePage)
      out << "<a class=\"descrheadentry\" href=\"" << userHomePage << "\">" << productName << "</a>" << std::endl;
   out << "<a class=\"descrheadentry\" href=\"http://root.cern.ch\">ROOT Homepage</a>" << std::endl
      << "<a class=\"descrheadentry\" href=\"./ClassIndex.html\">Class Index</a>" << std::endl
      << "<a class=\"descrheadentry\" href=\"./ClassHierarchy.html\">Class Hierarchy</a></div>" << std::endl;
   WriteSearch(out);
   out << "</div>" << std::endl; // descrhead, line 1

   if (withLocation) {
      out << "</div>" << std::endl; //toplinks
      WriteLocation(out, module, classname); // descrhead line 2
   }
   // else {
   //    Closed by caller!
   //    out << "</div>" << std::endl; // toplinks
   // }

}
