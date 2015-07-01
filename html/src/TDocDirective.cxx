#include "TDocDirective.h"

#include "TApplication.h"
#include "TClass.h"
#include "TDocInfo.h"
#include "TDocOutput.h"
#include "TDocParser.h"
#include "TError.h"
#include "THtml.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualPad.h"
#include "TVirtualMutex.h"
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <stdlib.h>

//______________________________________________________________________________
//
// When THtml parses documentation (through TDocParser), it checks for special
// words ("begin_something", "end_something", where the begin and end are the
// significant part). THtml then searches for a TDocDirective which can handle
// these tags ("whatever" in the example), passes the text enclosed by these
// tags to the directive, which in turn processes it.
//
// That way, HTML, latex, and C++ macros can be processed by THtml, e.g. to
// generate plain HTML or GIF pictures. The classes reposinsible for parsing
// that are TDocHtmlDirective, TDocLatexDirective, and TDocMacroDirective,
// respecively.
//
// Directives can have optional parameters; these are passed as paranthesis
// enclosed, comma delimited name=value pairs; see SetParameters().
//
// You can implement your own directive simply by deriving from TDocDirective;
// the tag corresponds to TDocDirective's name (e.g. "HTML" for "begin_html" /
// "end_html").
//______________________________________________________________________________

ClassImp(TDocDirective);

////////////////////////////////////////////////////////////////////////////////
/// Delete all output generated by the directive beginning
/// with Name() and ending with ext

void TDocDirective::DeleteOutputFiles(const char* ext) const
{
   TString basename;
   GetName(basename);
   basename += "_";
   TString dirname(GetOutputDir());
   void* hDir = gSystem->OpenDirectory(dirname);
   const char* entry = 0;
   while ((entry = gSystem->GetDirEntry(hDir))) {
      TString sEntry(entry);
      if (sEntry.BeginsWith(basename) && isdigit(sEntry[basename.Length()]) && (!ext || sEntry.EndsWith(ext)))
         gSystem->Unlink((dirname + "/" + entry).Data());
   }
   gSystem->FreeDirectory(hDir);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the full name, based on fName, fTitle, fDocParser's tag.

void TDocDirective::GetName(TString& name) const
{
   name = fName;
   if (fDocParser && fDocParser->GetCurrentClass()) {
      name += "_";
      TString outfilename;
      GetHtml()->GetHtmlFileName(fDocParser->GetCurrentClass(), outfilename);
      outfilename = gSystem->BaseName(outfilename);
      Ssiz_t posExt = outfilename.Last('.');
      outfilename.Remove(posExt, outfilename.Length() - posExt);
      name += outfilename;
   }
   if (GetTitle() && strlen(GetTitle())) {
      name += "_";
      name += GetTitle();
   }
   if (fCounter != -1) {
      name += "_";
      name += fCounter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the directory for documentation output.

const char* TDocDirective::GetOutputDir() const
{
   return fHtml ? fHtml->GetOutputDir().Data() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a string containing parameters in params,
/// we call AddParameter() for each of them.
/// This function splits the parameter names and
/// extracts their values if they are given.
/// Parameters are separated by ",", values are
/// separated from parameter names by "=".
/// params being
///    a = "a, b, c", b='d,e'
/// will issue two calls to AddParameter(), one for
/// a with value "a, b, c" and one for b with value
/// "d,e" (each without the quotation marks).

void TDocDirective::SetParameters(const char* params)
{
   fParameters = params;

   if (!fParameters.Length())
      return;

   TString param;
   Ssiz_t pos = 0;
   while (fParameters.Tokenize(param, pos, ",")) {
      param = param.Strip(TString::kBoth);
      if (!param.Length())
         continue;

      Ssiz_t posAssign = param.Index('=');
      if (posAssign != kNPOS) {
         TString value(param(posAssign + 1, param.Length()));
         value = value.Strip(TString::kBoth);
         if (value[0] == '\'')
            value = value.Strip(TString::kBoth, '\'');
         else if (value[0] == '"')
            value = value.Strip(TString::kBoth, '"');
         param.Remove(posAssign, param.Length());
         param = param.Strip(TString::kBoth);
         AddParameter(param, value);
      } else {
         param = param.Strip(TString::kBoth);
         AddParameter(param, 0);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the parser, and fDocOutput, fHtml from that

void TDocDirective::SetParser(TDocParser* parser)
{
   fDocParser    = parser;
   fDocOutput = parser ? parser->GetDocOutput() : 0;
   fHtml      = fDocOutput? fDocOutput->GetHtml() : 0;
}


//______________________________________________________________________________
//
// Process a "begin_html" / "end_html" block. Stop linking keywords and simply
// copy the text enclosed by the directive to the output HTML file.
//______________________________________________________________________________

ClassImp(TDocHtmlDirective);

////////////////////////////////////////////////////////////////////////////////
/// Add a line of HTML

void TDocHtmlDirective::AddLine(const TSubString& line)
{
   if (line.Start() == -1) return;

   TPRegexp pretag("</?[pP][rR][eE][ >]");
   TSubString iLine(line);
   Ssiz_t posPre = iLine.String().Index(pretag, iLine.Start());
   if (posPre == kNPOS)
      fText += line;
   else {
      // remove <pre> in fVerbatim environments, and
      // </pre> in !fVerbatim environments.
      while (posPre != kNPOS && posPre > 0) {
         Bool_t isOpen = line[posPre + 1 - line.Start()] != '/';
         Ssiz_t posClose = iLine.String().Index(">", posPre);
         if (posClose ==kNPOS) break; // aka oops.
         Ssiz_t len = posClose - posPre;

         if (fVerbatim) {
            if (isOpen) {
               // skip
               fText += iLine.String()(iLine.Start(), posPre - iLine.Start());
            } else {
               // write it out
               fText += iLine.String()(iLine.Start(), posPre + len - iLine.Start());
               fVerbatim = kFALSE;
            }
         } else {
            if (!isOpen) {
               // skip
               fText += iLine.String()(iLine.Start(), posPre - iLine.Start());
            } else {
               // write it out
               fText += iLine.String()(iLine.Start(), posPre + len - iLine.Start());
               fVerbatim = kTRUE;
            }
         }

         iLine = iLine.String()(posPre + len, iLine.Length());
         posPre = iLine.String().Index(pretag, iLine.Start());
      }

      fText += iLine;
   }
   fText += "\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Set result to the HTML code that was passed in via AddLine().
/// Prepend a closing </pre>, append an opening <pre>

Bool_t TDocHtmlDirective::GetResult(TString& result)
{
   result = "</pre><!-- TDocHtmlDirective start -->";
   result += fText + "<!-- TDocHtmlDirective end --><pre>";
   return kTRUE;
}



//______________________________________________________________________________
//
// Process a "begin_macro" / "end_macro" block. The block can be a file name
// or a CINT script (i.e. even ".x file.C" is allowed). See AddParameter() for
// supported options. Example (the quotes prevent THtml from expanding the
// example):
//
// "BEGIN_MACRO"
// .x $ROOTSYS/tutorials/hsimple.C
// "END_MACRO"
//
// The macro is meant to create an object that can be saved as a GIF file by
// calling object->SaveAs(outputfile.gif). The macro is expected to return that
// object as a TObject*; if it does not, gPad is used and saved. The object
// is deleted by TDocMacroDirective once saved.
//______________________________________________________________________________

ClassImp(TDocMacroDirective);

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TDocMacroDirective::~TDocMacroDirective()
{
   delete fMacro;
}

////////////////////////////////////////////////////////////////////////////////

void TDocMacroDirective::SubProcess(const TString& what, const TString& out) {
   Int_t error = TInterpreter::kNoError;
   Long_t ret = gROOT->ProcessLine(TString(".x ") + what, &error);
   Int_t sleepCycles = 50; // 50 = 5 seconds
   while (error == TInterpreter::kProcessing && --sleepCycles > 0)
      gSystem->Sleep(100);

   gSystem->ProcessEvents(); // in case ret needs to handle some events first

   if (error != TInterpreter::kNoError) {
      ::Error("TDocMacroDirective::HandleDirective_Macro",
              "Error processing macro for %s!", out.Data());
      return;
   }
   if (!ret) {
      return;
   }

   // Something with a vtable
   const TObject* objRet = (const TObject*)ret;
   try {
      typeid(*objRet).name(); // needed to test whether ret is indeed an object with a vtable!
      objRet = dynamic_cast<const TObject*>(objRet);
   }
   catch (...) {
      objRet = 0;
   }

   if (!objRet) {
      return;
   }

   if (gDebug > 3)
      ::Info("TDocMacroDirective::HandleDirective_Macro",
             "Saving returned %s to file %s.",
             objRet->IsA()->GetName(), out.Data());

   if (!gROOT->IsBatch()) {
      // to get X11 to sync :-( gVirtualX->Update()/Sync() don't do it
      gSystem->Sleep(1000);
      gVirtualX->Update(0);
      gVirtualX->Update(1);
   }

   gSystem->ProcessEvents();
   if (!gROOT->IsBatch()) {
      gVirtualX->Update(0);
      gVirtualX->Update(1);
   }

   objRet->SaveAs(out);
   gSystem->ProcessEvents(); // SaveAs triggers an event

#ifdef R__BEPAEPSTLICHERALSDERPAPST
   // ensure objRet is not e.g. the TGMainFrame of a new TCanvas: require padSave == gPad
   if (objRet != gPad && padSave == gPad)
      delete objRet;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Add a macro line.
/// Lines ending on "*HIDE*" will be executed as part of the
/// macro, but not shown in the source tab if the parameter
/// source is supplied.

void TDocMacroDirective::AddLine(const TSubString& line)
{
   if (!fMacro) {
      TString name;
      GetName(name);
      fMacro = new TMacro(name);
   }

   // return if no line - or if there was an intentinal line-break,
   // i.e. an empty line
   if (line.Start() == -1 && const_cast<TSubString&>(line).String().Length()) return;

   TString sLine(line);
   fMacro->AddLine(sLine);
   fIsFilename &= !sLine.Contains('{');
}

////////////////////////////////////////////////////////////////////////////////
/// Create the input file for SubProcess().

TString TDocMacroDirective::CreateSubprocessInputFile() {
   if (!fIsFilename) {
      TString fileSysName;
      GetName(fileSysName);
      fileSysName += ".C";
      gSystem->PrependPathName(gSystem->TempDirectory(), fileSysName);
      fMacro->SaveSource(fileSysName);
      return fileSysName;
   }

   // We have a filename; find it and build the invocation.
   TString filename;
   TIter iLine(fMacro->GetListOfLines());
   while (filename.Length() == 0)
      filename = ((TObjString*)iLine())->String().Strip(TString::kBoth);

   TString macroPath;
   TString modulename;
   if (GetHtml() && GetDocParser()) {
      if (GetDocParser()->GetCurrentClass())
         GetHtml()->GetModuleNameForClass(modulename, GetDocParser()->GetCurrentClass());
      else GetDocParser()->GetCurrentModule(modulename);
   }
   if (modulename.Length()) {
      GetHtml()->GetModuleMacroPath(modulename, macroPath);
   } else macroPath = gSystem->pwd();

   const char* pathDelimiter = ":"; // use ":" even on windows
   TObjArray* arrDirs(macroPath.Tokenize(pathDelimiter));
   TIter iDir(arrDirs);
   TObjString* osDir = 0;
   macroPath = "";
   TString filenameDirPart(gSystem->DirName(filename));
   filenameDirPart.Prepend('/'); // as dir delimiter, not as root dir
   while ((osDir = (TObjString*)iDir())) {
      if (osDir->String().EndsWith("\\"))
         osDir->String().Remove(osDir->String().Length() - 1);
      osDir->String() += filenameDirPart;
      macroPath += osDir->String() + pathDelimiter;
   }

   TString plusplus;
   while (filename.EndsWith("+")) {
      plusplus += '+';
      filename.Remove(filename.Length() - 1);
   }

   TString params;
   if (filename.EndsWith(")")) {
      Ssiz_t posOpen = filename.Last('(');
      if (posOpen != kNPOS) {
         params = filename(posOpen, filename.Length());
         filename.Remove(posOpen, filename.Length());
      }
   }

   TString fileSysName(gSystem->BaseName(filename));
   if (!gSystem->FindFile(macroPath, fileSysName)) {
      Error("GetResult", "Cannot find macro '%s' in path '%s'!",
            gSystem->BaseName(filename), macroPath.Data());
      return "";
   }
   fileSysName += params;
   fileSysName += plusplus;


   if (fShowSource) {
      // copy macro into fMacro - before running it, in case the macro blocks its file
      std::ifstream ifMacro(fileSysName);
      fMacro->GetListOfLines()->Delete();
      TString line;
      while (ifMacro) {
         if (!line.ReadLine(ifMacro, kFALSE) || ifMacro.eof())
            break;
         fMacro->AddLine(line);
      }
   }
   return fileSysName;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the result (i.e. an HTML img tag) for the macro invocation.
/// If fShowSource is set, a second tab will be created which shows
/// the source.

Bool_t TDocMacroDirective::GetResult(TString& result)
{
   if (!fMacro)
      return kFALSE;

   if (!fMacro->GetListOfLines()
      || !fMacro->GetListOfLines()->First()) {
      Warning("GetResult", "Empty directive found!");
      return kTRUE;
   }

   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   if (gDebug > 3)
      Info("HandleDirective_Macro", "executing macro \"%s\" with %d lines.",
         fMacro->GetName(), fMacro->GetListOfLines() ? fMacro->GetListOfLines()->GetEntries() + 1 : 0);

   Bool_t wasBatch = gROOT->IsBatch();
   Bool_t wantBatch = kFALSE;
   if (!wasBatch && !fNeedGraphics)
      wantBatch = kTRUE;
   else if (fNeedGraphics) {
      if (fHtml->IsBatch()) {
         Warning("GetResult()", "Will not initialize the graphics system; skipping macro %s!", GetName());
         result = "";
         return kFALSE;
      }
   }

   TString outFileName;
   {
      GetName(outFileName);
      GetDocOutput()->NameSpace2FileName(outFileName);
      outFileName += ".gif";
      outFileName.ReplaceAll(" ", "_");
      gSystem->PrependPathName(GetOutputDir(), outFileName);
   }

   TString subProcInputFile = CreateSubprocessInputFile();
   if (!subProcInputFile.Length()) return kFALSE;

   subProcInputFile.ReplaceAll("\\", "\\\\");
   subProcInputFile.ReplaceAll("\"", "\\\"");
   TString invoc("root.exe -l -q ");
   if (wantBatch) {
      invoc += "-b ";
   }
   invoc += "-e 'TDocMacroDirective::SubProcess(\""
      + subProcInputFile + "\",\"" + outFileName + "\");'";
   gSystem->Unlink(outFileName);
   Int_t exitCode = gSystem->Exec(invoc.Data());

   if (exitCode && gDebug > 0) {
      Info("GetResult()", "Subprocess exited with status %d\n", exitCode);
   } else if (!fIsFilename) {
      // we have created the input file.
      gSystem->Unlink(subProcInputFile);
   }

   if (!gSystem->AccessPathName(outFileName)) {
      // Output file was created
      result = "<span class=\"macro\"><img class=\"macro\" alt=\"output of ";
      result += outFileName;

      result += "\" title=\"MACRO\" src=\"";
      result += gSystem->BaseName(outFileName);
      result += "\" /></span>";
   }

   if (fShowSource) {
      // convert the macro source
      TIter iLine(fMacro->GetListOfLines());
      TObjString* osLine = 0;
      std::stringstream ssRaw;
      while ((osLine = (TObjString*)iLine()))
         ssRaw << osLine->String() << std::endl;

      TDocParser *dparser = 0;
      if (GetDocParser()->GetCurrentClass())
         dparser = new TDocParser(*(TClassDocOutput*)GetDocOutput(), GetDocParser()->GetCurrentClass());
      else dparser = new TDocParser(*GetDocOutput());
      std::stringstream ssConverted;
      dparser->Convert(ssConverted, ssRaw, "./", kTRUE /*code*/, kFALSE /*process directives*/);
      delete dparser;

      fMacro->GetListOfLines()->Delete();
      TString line;
      while (!ssConverted.fail()) {
         if (!line.ReadLine(ssConverted, kFALSE) || ssConverted.eof())
            break;
         fMacro->AddLine(line);
      }

      TString id(gSystem->BaseName(outFileName));
      id = id(0, id.Length()-4); // remove ".gif"
      // TODO: we need an accessible version of the source, i.e. visible w/o javascript
      TString tags("</pre><div class=\"tabs\">\n"
               "<a id=\"" + id + "_A0\" class=\"tabsel\" href=\"" + gSystem->BaseName(outFileName) + "\" onclick=\"javascript:return SetDiv('" + id + "',0);\">Picture</a>\n"
               "<a id=\"" + id + "_A1\" class=\"tab\" href=\"#\" onclick=\"javascript:return SetDiv('" + id + "',1);\">Source</a>\n"
               "<br /></div><div class=\"tabcontent\">\n"
               "<div id=\"" + id + "_0\" class=\"tabvisible\">" + result + "</div>\n"
               "<div id=\"" + id + "_1\" class=\"tabhidden\"><div class=\"listing\"><pre class=\"code\">");
      iLine.Reset();
      osLine = 0;
      while ((osLine = (TObjString*) iLine()))
         if (!TString(osLine->String().Strip()).EndsWith("*HIDE*"))
            tags += osLine->String() + "\n";
      if (tags.EndsWith("\n"))
         tags.Remove(tags.Length()-1); // trailing line break
      tags += "</pre></div></div><div class=\"clear\"></div></div><pre>";
      result = tags;
      // Protect the nested comments from being stripped by a
      // TDocParser::ProcessComment() in the call stack.
      result.ReplaceAll("<span class=\"comment\">", "<span class=\"codecomment\">");
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setting fNeedGraphics if name is "GUI",
/// setting fShowSource if name is "SOURCE"

void TDocMacroDirective::AddParameter(const TString& name, const char* /*value=0*/)
{
   if (!name.CompareTo("gui", TString::kIgnoreCase))
      fNeedGraphics = kTRUE;
   else if (!name.CompareTo("source", TString::kIgnoreCase))
      fShowSource = kTRUE;
   else Warning("AddParameter", "Unknown option %s!", name.Data());
}



namespace {
   Float_t gLinePadding = 10.; //px
   Float_t gColumnPadding = 10.; //px

   class TLatexLine {
   private:
      std::vector<Float_t> fWidths;
      Float_t fHeight;
      TObjArray* fColumns; // of TObjString*

   public:
      TLatexLine(TObjArray* columns = 0):
         fHeight(0.), fColumns(columns) { if (columns) fWidths.resize(Size());}

      Float_t& Width(UInt_t col) {return fWidths[col];}
      Float_t& Height() {return fHeight;}
      TString* operator[](Int_t column) {
         if (fColumns && fColumns->GetEntriesFast() > column)
            return &(((TObjString*)fColumns->At(column))->String());
         return 0;
      }
      UInt_t Size() const { return fColumns ? fColumns->GetEntries() : 0; }
      void Delete() { delete fColumns; }
   };
}

//______________________________________________________________________________
//
// Handle a "Begin_Latex"/"End_Latex" directive.
// called as
// "Begin_Latex(fontsize=10, separator='=,', rseparator='=|,', align=lcl)"
// will create and include a TLatex-processed image, with a given fontsize
// in pixels (defaults to 16). If (r)separator is given, the formulas on the
// following lines will be grouped into columns; a new column starts with
// (regexp) match of the separator; by default there is only one column.
// separator matches any character, rseparator matches as regexp with one
// column per pattern match. Only one of separator or rseparator can be given.
// align defines the alignment for each columns; be default, all columns
// are right aligned. NOTE that the column separator counts as a column itself!
//______________________________________________________________________________


ClassImp(TDocLatexDirective);

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TDocLatexDirective::~TDocLatexDirective()
{
   gSystem->ProcessEvents();
   delete fLatex;
   delete fBBCanvas;
   gSystem->ProcessEvents();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a latex line

void TDocLatexDirective::AddLine(const TSubString& line)
{
   if (line.Length() == 0)
      return;

   if (!fLatex) {
      TString name;
      GetName(name);
      fLatex = new TMacro(name);
   }

   TString sLine(line);
   GetDocParser()->Strip(sLine);
   if (sLine.Length() == 0)
      return;

   fLatex->AddLine(sLine);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a gif file named filename from a latex expression in fLatex.
/// Called when "Begin_Latex"/"End_Latex" is processed.

void TDocLatexDirective::CreateLatex(const char* filename)
{
   if (!fLatex
      || !fLatex->GetListOfLines()
      || !fLatex->GetListOfLines()->First())
      return;

   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   TVirtualPad* oldPad = gPad;

   Bool_t wasBatch = gROOT->IsBatch();
   if (!wasBatch)
      gROOT->SetBatch();

   const Float_t canvSize = 1200.;
   if (!fBBCanvas)
      // add magic batch vs. gui canvas sizes (4, 28)
      fBBCanvas = (TVirtualPad*)gROOT->ProcessLineFast(
         Form("new TCanvas(\"R__TDocLatexDirective_BBCanvas\",\"fBBCanvas\",%g,%g);", -(canvSize + 4.), canvSize + 28.));
   if (!fBBCanvas) {
      Error("CreateLatex", "Cannot create a TCanvas via the interpreter!");
      return;
   }
   fBBCanvas->SetBorderMode(0);
   fBBCanvas->SetFillColor(kWhite);

   gSystem->ProcessEvents();

   std::list<TLatexLine> latexLines;
   std::vector<Float_t> maxWidth(20);
   UInt_t numColumns = 0;
   Float_t totalHeight = gLinePadding;

   TLatex latex;
   latex.SetTextFont(43);
   latex.SetTextSize((Float_t)fFontSize);
   latex.SetTextAlign(12);

   // calculate positions
   TIter iterLine(fLatex->GetListOfLines());
   TObjString* line = 0;
   TPRegexp regexp;
   if (fSeparator.Length()) {
      if (fSepIsRegexp)
         regexp = TPRegexp(fSeparator);
   } else fSepIsRegexp = kFALSE;

   while ((line = (TObjString*) iterLine())) {
      const TString& str = line->String();
      TObjArray* split = 0;
      if (!fSepIsRegexp) {
         split = new TObjArray();
         split->SetOwner();
      }
      if (!fSeparator.Length())
         split->Add(new TObjString(str));
      else {
         if (fSepIsRegexp)
            split = regexp.MatchS(str);
         else {
            Ssiz_t prevStart = 0;
            for (Ssiz_t pos = 0; pos < str.Length(); ++pos) {
               if (fSeparator.Index(str[pos]) != kNPOS) {
                  split->Add(new TObjString(TString(str(prevStart, pos - prevStart))));
                  split->Add(new TObjString(TString(str(pos, 1))));
                  prevStart = pos + 1;
               }
            }
            split->Add(new TObjString(TString(str(prevStart, str.Length() - prevStart))));
         }
      }

      latexLines.push_back(TLatexLine(split));
      if (numColumns < (UInt_t)split->GetEntries())
         numColumns = split->GetEntries();

      Float_t heightLine = -1.;
      for (UInt_t col = 0; col < (UInt_t)split->GetEntries(); ++col) {
         Float_t widthLatex = 0.;
         Float_t heightLatex = 0.;
         TString* strCol = latexLines.back()[col];
         if (strCol)
            GetBoundingBox(latex, *strCol, widthLatex, heightLatex);
         if (heightLine < heightLatex)   heightLine = heightLatex;
         if (maxWidth.size() < col)
            maxWidth.resize(col * 2);
         if (maxWidth[col] < widthLatex)
            maxWidth[col] = widthLatex;
         latexLines.back().Width(col) = widthLatex;
      }
      latexLines.back().Height() = heightLine;
      totalHeight += heightLine + gLinePadding;
   } // while next line

   std::vector<Float_t> posX(numColumns + 1);
   for (UInt_t col = 0; col <= numColumns; ++col) {
      if (col == 0) posX[col] = gColumnPadding;
      else          posX[col] = posX[col - 1] + maxWidth[col - 1] + gColumnPadding;
   }
   Float_t totalWidth = posX[numColumns];

   // draw
   fBBCanvas->Clear();
   fBBCanvas->cd();
   Float_t padSizeX = totalWidth;
   Float_t padSizeY = totalHeight + 8.;
   // add magic batch vs. gui canvas sizes (4, 28) + rounding
   TVirtualPad* padImg = (TVirtualPad*)gROOT->ProcessLineFast(
      Form("new TCanvas(\"R__TDocLatexDirective_padImg\",\"padImg\",-(Int_t)%g,(Int_t)%g);",
           padSizeX + 4.5, padSizeY + 28.5));
   padImg->SetBorderMode(0);
   padImg->SetFillColor(kWhite);
   padImg->cd();

   Float_t posY = 0.;
   for (std::list<TLatexLine>::iterator iLine = latexLines.begin();
      iLine != latexLines.end(); ++iLine) {
      posY += iLine->Height()/2. + gLinePadding;
      for (UInt_t iCol = 0; iCol < iLine->Size(); ++iCol) {
         TString* str = (*iLine)[iCol];
         if (!str) continue;
         char align = 'l';
         if ((UInt_t)fAlignment.Length() > iCol)
            align = fAlignment[(Int_t)iCol];
         Float_t x = posX[iCol];
         switch (align) {
            case 'l': break;
            case 'r': x += maxWidth[iCol] - iLine->Width(iCol); break;
            case 'c': x += 0.5*(maxWidth[iCol] - iLine->Width(iCol)); break;
            default:
               if (iLine == latexLines.begin())
                  Error("CreateLatex", "Invalid alignment character '%c'!", align);
         }
         latex.DrawLatex( x / padSizeX, 1. - posY / padSizeY, str->Data());
      }
      posY += iLine->Height()/2.;
   }

   padImg->Print(filename);

   // delete the latex objects
   for (std::list<TLatexLine>::iterator iLine = latexLines.begin();
      iLine != latexLines.end(); ++iLine) {
      iLine->Delete();
   }

   delete padImg;

   if (!wasBatch)
      gROOT->SetBatch(kFALSE);

   gPad = oldPad;
}

////////////////////////////////////////////////////////////////////////////////
/// Determines the bounding box for text as height and width.
/// Assumes that we are in batch mode.

void TDocLatexDirective::GetBoundingBox(TLatex& latex, const char* text, Float_t& width, Float_t& height)
{
   UInt_t uiWidth = 0;
   UInt_t uiHeight = 0;
   fBBCanvas->cd();
   latex.SetText(0.1, 0.5, text);
   latex.GetBoundingBox(uiWidth, uiHeight);

   width = uiWidth;
   height = uiHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list of lines as TObjStrings

TList* TDocLatexDirective::GetListOfLines() const
{
   return fLatex ? fLatex->GetListOfLines() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// convert fLatex to a gif by creating a TLatex, drawing it on a
/// temporary canvas, and saving that to a filename in the output
/// directory.

Bool_t TDocLatexDirective::GetResult(TString& result)
{
   TString filename;
   GetName(filename);
   filename.ReplaceAll(" ", "_");
   const TString& firstLine = ((TObjString*)fLatex->GetListOfLines()->First())->String();
   TString latexFilename(firstLine);
   for (Ssiz_t namepos = 0; namepos < latexFilename.Length(); ++namepos)
      if (!GetDocParser()->IsWord(latexFilename[namepos])) {
         latexFilename.Remove(namepos, 1);
         --namepos;
      }
   filename += "_";
   filename += latexFilename;

   GetDocOutput()->NameSpace2FileName(filename);
   filename += ".gif";

   TString altText(firstLine);
   GetDocOutput()->ReplaceSpecialChars(altText);
   altText.ReplaceAll("\"", "&quot;");
   result = "<span class=\"latex\"><img class=\"latex\" alt=\"";
   result += altText;
   result += "\" title=\"LATEX\" src=\"";
   result += filename;
   result += "\" /></span>";

   gSystem->PrependPathName(GetOutputDir(), filename);

   if (gDebug > 3)
      Info("HandleDirective_Latex", "Writing Latex \"%s\" to file %s.",
           fLatex->GetName(), filename.Data());

   CreateLatex(filename);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse fParameters, setting fFontSize, fAlignment, and fSeparator

void TDocLatexDirective::AddParameter(const TString& name, const char* value /*=0*/)
{
   if (!name.CompareTo("fontsize", TString::kIgnoreCase)) {
      if (!value || !value[0])
         Error("AddParameter", "Option \"fontsize\" needs a value!");
      else fFontSize = atol(value);
   } else if (!name.CompareTo("separator", TString::kIgnoreCase)) {
      if (!value || !value[0])
         Error("AddParameter", "Option \"separator\" needs a value!");
      else fSeparator = value;
   } else if (!name.CompareTo("align", TString::kIgnoreCase)) {
      if (!value || !value[0])
         Error("AddParameter", "Option \"align\" needs a value!");
      else fAlignment = value;
   } else
      Warning("AddParameter", "Unknown option %s!", name.Data());
}
