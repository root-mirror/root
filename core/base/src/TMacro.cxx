// @(#)root/base:$Id$
// Author: Rene Brun   16/08/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TMacro
\ingroup Base

Class supporting a collection of lines with C++ code.
A TMacro can be executed, saved to a ROOT file, edited, etc.

A macro can be built line by line by calling the AddLine function.
or it can be created directly from a file via the special constructor
when the first argument is a file name.

A macro can be executed via the Exec function.
Arguments can be specified when calling Exec.

A macro can be drawn in a pad. When the pad is updated, the macro is
automatically executed.

The code in the macro can be saved via the SaveSource function.
If the macro is in the list of primitives of a pad/canvas, the macro
will be saved in the script generated by TCanvas::SaveSource.

A macro can be written to a ROOT file via TObject::Write.

Examples:
~~~ {.cpp}
  TMacro m("Peaks.C");  //macro m with name "Peaks" is created
                        //from file  Peaks.C
  m.Exec();             //macro executed with default arguments
  m.Exec("4");          //macro executed with argument
  m.SaveSource("newPeaks.C");
  TFile f("mymacros.root","recreate");
  m.Write();   //macro saved to file with name "Peaks"
~~~
*/

#include "TEnv.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TMacro.h"
#include "TMD5.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "strlcpy.h"
#include <iostream>
#include <fstream>
#include <sstream>

ClassImp(TMacro);

////////////////////////////////////////////////////////////////////////////////
/// Create an empty macro, use AddLine() or ReadFile() to fill this macro.

TMacro::TMacro(): TNamed(), fLines(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a macro with a name and a title.
/// If name contains a '.' it is assumed to be the name of a file, and
///  - the macro is automatically filled by reading all the lines in the file,
///  - if the title is empty, it will be set to the name of the file,
///  - the name will be set to the filename without path or extension.

TMacro::TMacro(const char *name, const char *title): TNamed(name,title)
{
   fLines  = new TList();
   if (!name) return;
   Int_t nch = strlen(name);
   char *s = new char[nch+1];
   strlcpy(s,name,nch+1);
   char *slash = (char*)strrchr(s,'/');
   if (!slash) slash = s;
   else ++slash;
   char *dot   = (char*)strchr(slash,'.');
   if (dot) {
      *dot = 0;
      fName = slash;
      if (fTitle.Length() == 0) fTitle = name;
      ReadFile(name);
   }
   delete [] s;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TMacro::TMacro(const TMacro &macro): TNamed(macro)
{
   fLines = new TList();
   TIter next(macro.GetListOfLines());
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      fLines->Add(new TObjString(obj->GetName()));
   }
   fParams = macro.fParams;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete this macro.

TMacro::~TMacro()
{
   if (fLines) fLines->Delete();
   delete fLines;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TMacro& TMacro::operator=(const TMacro &macro)
{
   if(this!=&macro) {
      TNamed::operator=(macro);
      if (fLines) fLines->Delete();
      delete fLines;
      fLines = new TList();
      TIter next(macro.GetListOfLines());
      TObjString *obj;
      while ((obj = (TObjString*) next())) {
         fLines->Add(new TObjString(obj->GetName()));
      }
      fParams = macro.fParams;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Add line with text in the list of lines of this macro.

TObjString *TMacro::AddLine(const char *text)
{
   if (!fLines) fLines = new TList();
   TObjString *obj = new TObjString(text);
   fLines->Add(obj);
   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// When clicking in the browser, the following action is performed
/// on this macro, depending the content of the variable TMacro.Browse.
/// TMacro.Browse can be set in the system.rootrc or .rootrc file like:
/// ~~~ {.cpp}
///     TMacro.Browse   :  Action
/// ~~~
/// or set via gEnv->SetValue, eg
/// ~~~ {.cpp}
///     gEnv->SetValue("TMacro.Browse","Print");
/// ~~~
/// By default TMacro.Browse=""
/// -if TMacro.Browse ="" the macro is executed
/// -if TMacro.Browse ="Print" the macro is printed in stdout
/// -if TMacro.Browse is of the form "mymacro.C"
///     the macro void mymacro.C(TMacro *m) is called where m=this macro
///     An example of macro.C saving the macro into a file and viewing it
///     with emacs is shown below:
/// ~~~ {.cpp}
///        void mymacro(TMacro *m) {
///           m->SaveSource("xx.log");
///           gSystem->Exec("emacs xx.log&");
///        }
/// ~~~

void TMacro::Browse(TBrowser * /*b*/)
{
   TString opt = gEnv->GetValue("TMacro.Browse","");
   if (opt.IsNull()) {
      Exec();
      return;
   }
   if (opt == "Print") {
      Print();
      return;
   }
   if (opt.Contains(".C")) {
      const char *cmd = Form(".x %s((TMacro*)0x%lx)",opt.Data(),(ULong_t)this);
      gROOT->ProcessLine(cmd);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns checksum of the current content. The returned TMD5 object must
/// be deleted by the user. Returns 0 in case of error.

TMD5 *TMacro::Checksum()
{
   if (!fLines || fLines->GetSize() <= 0)
      return (TMD5 *)0;

   TMD5 *md5 = new TMD5;

   // Fill (same params as in TMD5::FileChecksum)
   const Int_t bufSize = 8192;
   UChar_t buf[bufSize];
   Long64_t pos = 0;
   Long64_t left = bufSize;

   TIter nxl(fLines);
   TObjString *l;
   while ((l = (TObjString *) nxl())) {
      TString line = l->GetString();
      line += '\n';
      Int_t len = line.Length();
      char *p = (char *) line.Data();
      if (left > len) {
         strlcpy((char *)&buf[pos], p, len+1);
         pos += len;
         left -= len;
      } else if (left == len) {
         strlcpy((char *)&buf[pos], p, len+1);
         md5->Update(buf, bufSize);
         pos = 0;
         left = bufSize;
      } else {
         strlcpy((char *)&buf[pos], p, left+1);
         md5->Update(buf, bufSize);
         len -= left;
         p += left;
         pos = 0;
         left = bufSize;
         strlcpy((char *)&buf[pos], p, len+1);
         pos += len;
         left -= len;
      }
   }
   md5->Update(buf, pos);

   // Finalize
   md5->Final();

   return md5;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the macro into the interpreter.
/// Return true in case the loading was successful.

Bool_t TMacro::Load() const
{
   std::stringstream ss;

   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      ss << obj->GetName() << std::endl;
   }
   return gInterpreter->LoadText(ss.str().c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// Execute this macro with params, if params is 0, default parameters
/// (set via SetParams) are used.
/// error is set to an TInterpreter::EErrorCode by TApplication::ProcessLine().
/// Returns the result of the macro (return value or value of the last
/// expression), cast to a Long_t.

Long_t TMacro::Exec(const char *params, Int_t* error)
{
   if ( !gROOT->GetGlobalFunction(GetName(), 0, kTRUE) ) {
      if (!Load()) {
         if (error) *error = 1;
         return 0;
      }
   }

   // if macro has been executed, look for global function with name
   // of macro and re-execute this global function, if not found then
   // macro is unnamed macro, which we re-execute from file
   if ( gROOT->GetGlobalFunction(GetName(), 0, kTRUE) ) {
      gROOT->SetExecutingMacro(kTRUE);
      TString exec = GetName();
      TString p = params;
      if (p == "") p = fParams;
      if (p != "")
         exec += "(" + p + ")";
      else
         exec += "()";
      Long_t ret = gROOT->ProcessLine(exec, error);
      //enable gROOT->Reset
      gROOT->SetExecutingMacro(kFALSE);
      return ret;
   }

   Error("Exec","Macro does not contains function named %s.",GetName());
   if (error) *error = 1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Search the first line containing text.

TObjString *TMacro::GetLineWith(const char *text) const
{
   if (!fLines) return 0;
   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      if (strstr(obj->GetName(),text)) return obj;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute this macro (called by TPad::Paint).

void TMacro::Paint(Option_t *option)
{
   Exec(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print contents of this macro.

void TMacro::Print(Option_t * /*option*/) const
{
   if (!fLines) return;
   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      printf("%s\n",obj->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read lines in filename in this macro.

Int_t TMacro::ReadFile(const char *filename)
{
   if (!fLines) fLines = new TList();
   std::ifstream in;
   in.open(filename);
   if (!in.good()) {
      Error("ReadFile","Cannot open file: %s",filename);
      return 0;
   }
   char *line = new char[10000];
   Int_t nlines = 0;
   while (1) {
      in.getline(line,10000);
      if (!in.good()) break;
      if (in.eof()) break;
      fLines->Add(new TObjString(line));
      nlines++;
   }
   delete [] line;
   return nlines;
}

////////////////////////////////////////////////////////////////////////////////
/// Save macro source in filename.

void TMacro::SaveSource(const char *filename)
{
   std::ofstream out;
   out.open(filename, std::ios::out);
   if (!out.good ()) {
      Printf("SaveSource cannot open file: %s",filename);
      return;
   }
   if (!fLines) {out.close(); return;}
   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      out<<obj->GetName()<<std::endl;
   }
   out.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Save macro source in file pointer fp.

void TMacro::SaveSource(FILE *fp)
{
   if (!fLines) {fclose(fp); return;}
   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      fprintf(fp, "%s\n", obj->GetName());
   }
   fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
/// Save macro source on stream out.

void TMacro::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TMacro::Class())) {
      out<<"   ";
   } else {
      out<<"   "<<ClassName()<<" *";
   }
   out<<"macro = new "<<ClassName()<<"("<<quote<<GetName()<<quote<<","<<quote<<GetTitle()<<quote<<");"<<std::endl;
   if (!fLines) return;
   TIter next(fLines);
   TObjString *obj;
   while ((obj = (TObjString*) next())) {
      TString s = obj->GetName();
      s.ReplaceAll("\"","\\\"");
      out<<"   macro->AddLine("<<quote<<s.Data()<<quote<<");"<<std::endl;
   }
   out<<"   macro->Draw("<<quote<<option<<quote<<");"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default parameters to execute this macro.

void TMacro::SetParams(const char *params)
{
   if (params) fParams = params;
}
