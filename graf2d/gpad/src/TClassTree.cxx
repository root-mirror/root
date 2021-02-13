// @(#)root/gpad:$Id$
// Author: Rene Brun   01/12/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h"

#include "TROOT.h"
#include "TBuffer.h"
#include "TClassTree.h"
#include "TClassTable.h"
#include "TClass.h"
#include "TBaseClass.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TVirtualPad.h"
#include "TPaveClass.h"
#include "TArrow.h"
#include "TSystem.h"
#include "TObjString.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <fstream>
#include <algorithm>

const Int_t kIsClassTree = BIT(7);
const Int_t kUsedByData  = BIT(11);
const Int_t kUsedByFunc  = BIT(12);
const Int_t kUsedByCode  = BIT(13);
const Int_t kUsedByClass = BIT(14);
const Int_t kUsingData   = BIT(15);
const Int_t kUsingFunc   = BIT(16);
const Int_t kUsingCode   = BIT(17);
const Int_t kUsingClass  = BIT(18);
const Int_t kUsedByCode1 = BIT(19);
const Int_t kIsaPointer  = BIT(20);
const Int_t kIsBasic     = BIT(21);

static Float_t gXsize, gYsize, gDx, gDy, gLabdx, gLabdy, gDxx, gCsize;
static Int_t *gNtsons, *gNsons;

ClassImp(TClassTree);

/** \class TClassTree
\ingroup gpad

Draw inheritance tree and their relations for a list of classes.

The following options are supported
  - Direct inheritance (default)
  - Multiple inheritance
  - Composition
  - References by data members and member functions
  - References from Code

The list of classes is specified:
  - either in the TClassTree constructor as a second argument
  - or the parameter to TClassTree::Draw

Note that the ClassTree viewer can also be started from the canvas
pull down menu "Classes".

In the list of classes, class names are separated by a ":"
wildcarding is supported.
The following formats are supported, eg in TClassTree::Draw
  1. `Draw("ClassA")`
        - Draw inheritance tree for ClassA
        - Show all classes referenced by ClassA
  2. `Draw("*ClassB")`
        - Draw inheritance tree for ClassB
          and all the classes deriving from ClassB
  3. `Draw(">ClassC")`
        - Draw inheritance tree for ClassC
        - Show classes referencing ClassC
  4. `Draw("ClassD<")`
        - Draw inheritance tree for ClassD
        - Show classes referenced by ClassD
        - Show all classes referencing ClassD
  5. `Draw("Cla*")`
        - Draw inheritance tree for all classes with name starting with "Cla"
        - Show classes referenced by these classes
  6. `Draw("ClassA:ClassB<")`
        - Draw inheritance tree for ClassA
        - Show all classes referenced by ClassA
        - Draw inheritance tree for ClassB
        - Show classes referenced by ClassB
        - Show all classes referencing ClassB

Example:  `Draw("TTree<")`
       - Draw inheritance tree for the Root class TTree
       - Show all classes referenced by TTree
       - Show all classes using TTree

By default, only direct inheritance is drawn.
Use TClassTree::ShowLinks(option) to show additional references
 - option = "H" to show links to embedded classes
 - option = "M" to show multiple inheritance
 - option = "R" to show pointers to other classes from data members
 - option = "C" to show classes used by the code(implementation) of a class

The following picture is produced directly by:
~~~ {.cpp}
      TClassTree ct("ct","*TH1")
~~~
It shows all the classes derived from the base class TH1.

\image html gpad_classtree1.png

The TClassTree class uses the services of the class TPaveClass to
show the class names. By clicking with the right mouse button in
one TPaveClass object, one can invoke the following functions of TClassTree:
  - ShowLinks(option) with by default option = "HMR"
  - Draw(classes). By default the class drawn is the one being pointed
  - ShowClassesUsedBy(classes) (by default the pointed class)
  - ShowClassesUsing(classes) (by default the pointed class)

The following picture has been generated with the following statements
~~~ {.cpp}
      TClassTree tc1("tc1","TH1");
      tc1.ShowLinks("HMR");
~~~

\image html gpad_classtree2.png

Note that in case of embedded classes or pointers to classes,
the corresponding dashed lines or arrows respectively start
in the TPaveClass object at an X position reflecting the position
in the list of data members.

 - References by data members to other classes are show with a full red line
 - Multiple inheritance is shown with a dashed blue line
 - "Has a" relation is shown with a dotted cyan line
 - References from code is shown by a full green line

Use TClassTree::SetSourceDir to specify the search path for source files.
By default the search path includes the `$ROOTSYS` directory, the current
directory and the subdirectory `src`.

The first time TClassTree::Draw is invoked, all the classes in the
current application are processed, including the parsing of the code
to find all classes referenced by the include statements.
This process may take a few seconds. The following commands will be
much faster.

A TClassTree object may be saved in a Root file.
This object can be processed later by a Root program that ignores
the original classes. This interesting possibility allows to send
the class structure of an application to a colleague who does not have
your classes.

Example:
~~~ {.cpp}
   TFile f("myClasses.root","recreate")
   TClassTree *ct = new TClassTree("ct","ATLF*")
   ct->Write();
~~~
You can send at this point the file myClass.root to a colleague who can
run the following Root basic session
~~~ {.cpp}
    TFile f("myClass.root"); //connect the file
    tt.ls();                 //to list all classes and titles
    tt.Draw("ATLFDisplay")   //show class ATLFDisplay with all its dependencies
~~~
At this point, one has still access to all the classes present
in the original session and select any combination of these classes
to be displayed.
*/

////////////////////////////////////////////////////////////////////////////////
/// TClassTree default constructor.

TClassTree::TClassTree()
{
   fShowCod  = 0;
   fShowHas  = 0;
   fShowMul  = 0;
   fShowRef  = 0;
   fNclasses = 0;
   fCstatus  = 0;
   fParents  = 0;
   fCparent  = 0;
   fCpointer = 0;
   fCnames   = 0;
   fCtitles  = 0;
   fOptions  = 0;
   fLinks    = 0;
   fDerived  = 0;
   fNdata    = 0;
   SetLabelDx();
   SetYoffset(0);
   SetSourceDir(".:src:" + TROOT::GetSourceDir());
}

////////////////////////////////////////////////////////////////////////////////
/// TClassTree constructor.

TClassTree::TClassTree(const char *name, const char *classes)
           :TNamed(name,classes)
{
   fShowCod  = 0;
   fShowHas  = 0;
   fShowMul  = 0;
   fShowRef  = 0;
   fNclasses = 0;
   fCstatus  = 0;
   fParents  = 0;
   fCparent  = 0;
   fCpointer = 0;
   fCnames   = 0;
   fCtitles  = 0;
   fOptions  = 0;
   fLinks    = 0;
   fDerived  = 0;
   fNdata    = 0;
   SetLabelDx();
   SetYoffset(0);
   SetSourceDir(".:src:" + TROOT::GetSourceDir());

   // draw list of classes (if specified)
   if (classes && strlen(classes)) {
      fClasses = classes;
      Draw();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TClassTree default destructor.

TClassTree::~TClassTree()
{
   for (Int_t i=0;i<fNclasses;i++) {
      //delete fOptions[i];
      if (fLinks[i]) fLinks[i]->Delete();
      //delete fLinks[i];
      //if (fDerived[i]) {delete [] fDerived[i]; fDerived[i] = 0;}
   }
   delete [] fCnames;
   delete [] fCtitles;
   delete [] fCstatus;
   delete [] fParents;
   delete [] fCparent;
   delete [] fCpointer;
   delete [] fOptions;
   delete [] fLinks;
   delete [] fDerived;
   delete [] fNdata;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the inheritance tree and relations for the list of classes
/// see this class header for the syntax and examples

void TClassTree::Draw(const char *classes)
{
   if (!gPad) {
      gROOT->MakeDefCanvas();
   }
   Init();
   if (classes && strlen(classes)) fClasses = classes;
   for (Int_t i=0;i<fNclasses;i++) {
      fCstatus[i]  = 0;
      fCparent[i] = -1;
   }
   Paint();
}

////////////////////////////////////////////////////////////////////////////////
///  Find class number corresponding to classname in list of local classes

Int_t TClassTree::FindClass(const char *classname)
{
   for (Int_t i=0;i<fNclasses;i++) {
      if(!fCnames[i]->CompareTo(classname)) return i;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
///  Select all classes used/referenced by the class number iclass

void TClassTree::FindClassesUsedBy(Int_t iclass)
{
   fCstatus[iclass] = 1;
   Int_t i;
   TObjString *os;
   TList *los = fLinks[iclass];
   TIter next(los);
   while ((os = (TObjString*)next())) {
      i = FindClass(os->GetName());
      if (i < 0) continue;
      if (fCstatus[i]) continue;
      Int_t udata  = os->TestBit(kUsedByData);
      Int_t ufunc  = os->TestBit(kUsedByFunc);
      Int_t ucode  = os->TestBit(kUsedByCode);
      Int_t uclass = os->TestBit(kUsedByClass);
      if (udata || ufunc || ucode || uclass) {
         fCstatus[i] = 1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Select all classes using/referencing the class number iclass

void TClassTree::FindClassesUsing(Int_t iclass)
{
   // loop on all classes
   fCstatus[iclass] = 1;
   Int_t i;
   TObjString *os;
   TList *los = fLinks[iclass];
   TIter next(los);
   while ((os = (TObjString*)next())) {
      i = FindClass(os->GetName());
      if (i < 0) continue;
      if (fCstatus[i]) continue;
      Int_t udata  = os->TestBit(kUsingData);
      Int_t ufunc  = os->TestBit(kUsingFunc);
      Int_t ucode  = os->TestBit(kUsingCode);
      Int_t uclass = os->TestBit(kUsingClass);
      if (udata || ufunc || ucode || uclass) {
         fCstatus[i] = 1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Search the TPaveClass object in the pad with label=classname
/// returns the x and y position of the center of the pave.

void TClassTree::FindClassPosition(const char *classname, Float_t &x, Float_t &y)
{
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   TPaveClass *pave;
   while((obj=next())) {
      if (obj->InheritsFrom(TPaveClass::Class())) {
         pave = (TPaveClass*)obj;
         if (!strcmp(pave->GetLabel(),classname)) {
            x = 0.5*(pave->GetX1() + pave->GetX2());
            y = 0.5*(pave->GetY1() + pave->GetY2());
            return;
         }
      }
   }
   x = y = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the data structures

void TClassTree::Init()
{
   if (fNclasses) return;

   // fill the classes structures
   gClassTable->Init();
   fNclasses   = gClassTable->Classes();   //number of classes in the application
   fCnames     = new TString*[fNclasses];  //class names
   fCtitles    = new TString*[fNclasses];  //class titles (given in ClassDef)
   fCstatus    = new Int_t[fNclasses];     //=0 if not used in current expression
   fParents    = new Int_t[fNclasses];     //parent number of classes (permanent)
   fCparent    = new Int_t[fNclasses];     //parent number of classes (local to expression)
   fNdata      = new Int_t[fNclasses];     //number of data members per class
   fCpointer   = new TClass*[fNclasses];   //pointers to the TClass
   fOptions    = new TString*[fNclasses];  //options per class
   fLinks      = new TList*[fNclasses];    //list of classes referencing/referenced
   fDerived    = new char*[fNclasses];     //derivation matrix

   Int_t i,j;
   for (i=0;i<fNclasses;i++) {
      fCnames[i]   = new TString(gClassTable->Next());
      fCpointer[i] = TClass::GetClass(fCnames[i]->Data());
      fCtitles[i]  = new TString(fCpointer[i]->GetTitle());
      fCstatus[i]  = 0;
      fOptions[i]  = new TString("ID");
      fLinks[i]    = new TList();
      fDerived[i]  = new char[fNclasses];
   }
   TBaseClass *clbase;
   TClass *cl;
   for (i=0;i<fNclasses;i++) {
      TList *lm = fCpointer[i]->GetListOfDataMembers();
      if (lm) fNdata[i] = lm->GetSize();
      else    fNdata[i] = 0;
      // build derivation matrix
      char *derived = fDerived[i];
      for (j=0;j<fNclasses;j++) {
         derived[j] = 0;
         if (fCpointer[i]->InheritsFrom(fCpointer[j])) {
            derived[j] = 1;
         }
      }
      //build list of class parent
      fParents[i] = -1;
      TList *lb = fCpointer[i]->GetListOfBases();
      if (!lb) continue;
      clbase = (TBaseClass*)lb->First();
      if (clbase == 0) continue;
      cl = (TClass*)clbase->GetClassPointer();
      for (j=0;j<fNclasses;j++) {
         if(cl == fCpointer[j]) {
            fParents[i] = j;
            break;
         }
      }
   }
   //now the real & hard stuff
   for (i=0;i<fNclasses;i++) {
      ScanClasses(i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// list classes names and titles

void TClassTree::ls(Option_t *) const
{
   char line[500];
   for (Int_t i=0;i<fNclasses;i++) {
      snprintf(line,500,"%s%s",fCnames[i]->Data(),"...........................");
      snprintf(&line[30],460,"%s",fCtitles[i]->Data());
      line[79] = 0;
      printf("%5d %s\n",i,line);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// set bit abit in class classname in list los

TObjString *TClassTree::Mark(const char *classname, TList *los, Int_t abit)
{
   if (!los) return 0;
   TObjString *os = (TObjString*)los->FindObject(classname);
   if (!os) {
      os = new TObjString(classname);
      los->Add(os);
   }
   os->SetBit(abit);
   return os;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the current class setting in fClasses and fStatus

void TClassTree::Paint(Option_t *)
{
   //delete primitives belonging to a previous paint
   if (gPad) {
      TIter next(gPad->GetListOfPrimitives());
      TObject *obj;
      while((obj=next())) {
         if (obj->TestBit(kIsClassTree)) delete obj;
      }
   }

   Int_t nch      = strlen(GetClasses());
   if (nch == 0) return;
   char *classes  = new char[nch+1];
   gNsons   = new Int_t[fNclasses];
   gNtsons  = new Int_t[fNclasses];
   strlcpy(classes,GetClasses(),nch+1);
   Int_t i,j;
   char *derived;
   char *ptr = strtok(classes,":");
    //mark referenced classes
   while (ptr) {
      nch = strlen(ptr);
      if (ptr[0] == '*') {
         j = FindClass(&ptr[1]);
         if (j >= 0) {
            for (i=0;i<fNclasses;i++) {
               derived = fDerived[i];
               if(derived[j]) fCstatus[i] = 1;
            }
         }
      } else if (ptr[0] == '>') {
         for (i=0;i<fNclasses;i++) {
            if(fCnames[i]->Contains(&ptr[1])) {
               FindClassesUsing(i);
               fCstatus[i] = 2;
               break;
            }
         }
      } else if (ptr[nch-1] == '<') {
         ptr[nch-1] = 0;
         for (i=0;i<fNclasses;i++) {
            if(fCnames[i]->Contains(ptr)) {
               FindClassesUsedBy(i);
               FindClassesUsing(i);
               fCstatus[i] = 2;
               break;
            }
         }
      } else if (ptr[nch-1] == '*') {
         ptr[nch-1] = 0;
         for (i=0;i<fNclasses;i++) {
            if(fCnames[i]->Contains(ptr)) fCstatus[i] = 1;
         }
      } else {
         for (i=0;i<fNclasses;i++) {
            if(!fCnames[i]->CompareTo(ptr)) {
               FindClassesUsedBy(i);
               fCstatus[i] = 2;
               break;
            }
         }
      }
      ptr = strtok(0,":");
   }
    //mark base classes of referenced classes
   for (i=0;i<fNclasses;i++) {
      gNsons[i] = gNtsons[i] = 0;
   }
   for (i=0;i<fNclasses;i++) {
      if (fCstatus[i] == 0) continue;
      derived = fDerived[i];
      for (j=0;j<fNclasses;j++) {
         if (j == i) continue;
         if(derived[j]) {
            fCstatus[j] = 1;
         }
      }
   }
    //find parent class number for selected classes
   for (i=0;i<fNclasses;i++) {
      if (fCstatus[i] == 0) continue;
      j = fParents[i];
      if (j >=0 ) {
         fCparent[i] = j;
         gNsons[j]++;
      }
   }
    //compute total number of sons for each node
   Int_t maxlev = 1;
   Int_t icl,ip;
   for (i=0;i<fNclasses;i++) {
      if (fCstatus[i] == 0) continue;
      if (gNsons[i] != 0) continue;
      icl = i;
      Int_t nlevel = 1;
      while (fCparent[icl] >= 0) {
         nlevel++;
         if (nlevel > maxlev) maxlev = nlevel;
         ip = fCparent[icl];
         gNtsons[ip]++;
         icl = ip;
      }
   }

    //compute levels, number and list of sons
   Int_t ndiv=0;
   Int_t nmore = 0;
   for (i=0;i<fNclasses;i++) {
      if (fCstatus[i] == 0) continue;
      if (fCparent[i] < 0) {
         ndiv += gNtsons[i]+1;
         nmore++;
      }
   }
   ndiv++;

   // We are now ready to draw the active nodes
   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t ymin = gPad->GetY1();
   Float_t ymax = gPad->GetY2();
   Float_t ytop = gYsize/20;
   gXsize = xmax - xmin;
   gYsize = ymax - ymin;
   gDy = (gYsize-ytop)/(ndiv);
   if (gDy > gYsize/10.) gDy = gYsize/10.;
   gDx = 0.9*gXsize/5;
   if (maxlev > 5) gDx = 0.97*gXsize/maxlev;
   Float_t y  = ymax -ytop;
   gLabdx = fLabelDx*gXsize;
   if (gLabdx > 0.95*gDx) gLabdx = 0.95*gDx;
   gLabdy = 0.3*gDy;
   gDxx = 0.5*gXsize/26.;
   Float_t xleft  = xmin +gDxx;
   Float_t ymore  = 0.5*nmore*gDy+fYoffset*gYsize;
   Int_t dxpixels = gPad->XtoAbsPixel(gLabdx) - gPad->XtoAbsPixel(0);
   Int_t dypixels = gPad->YtoAbsPixel(0)     - gPad->YtoAbsPixel(gLabdy);
   gCsize  = dxpixels/(10.*dypixels);
   gCsize = std::max(gCsize,Float_t(0.75));
   gCsize = std::min(gCsize,Float_t(1.1));
   // draw classes level 0
   for (i=0;i<fNclasses;i++) {
      if (fCstatus[i] == 0) continue;
      if (fCparent[i] < 0) {
         y -= gDy+0.5*gNtsons[i]*gDy;
         if (!fCnames[i]->CompareTo("TObject")) y += ymore;
         PaintClass(i,xleft,y);
         y -= 0.5*gNtsons[i]*gDy;
      }
   }

   // show all types of links corresponding to selected options
   if (fShowCod) ShowCod();
   if (fShowHas) ShowHas();
   if (fShowMul) ShowMul();
   if (fShowRef) ShowRef();

   nch      = strlen(GetClasses());
   xmax = 0.3;
   if (nch > 20) xmax = 0.5;
   if (nch > 50) xmax = 0.7;
   if (nch > 70) xmax = 0.9;
   TPaveClass *ptitle = new TPaveClass(xmin +0.1*gXsize/26.
                                      ,ymin+gYsize-0.9*gYsize/20.
                                      ,xmin+xmax*gXsize
                                      ,ymin+gYsize-0.1*gYsize/26.
                                      ,GetClasses(),this);
   ptitle->SetFillColor(42);
   ptitle->SetBit(kIsClassTree);
   ptitle->Draw();

   //cleanup
   delete [] classes;
   delete [] gNsons;
   delete [] gNtsons;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint one class level

void TClassTree::PaintClass(Int_t iclass, Float_t xleft, Float_t y)
{
   Float_t u[2],yu=0,yl=0;
   Int_t ns = gNsons[iclass];
   u[0] = xleft;
   u[1] = u[0]+gDxx;
   if(ns != 0) u[1] = u[0]+gDx;
   TLine *line = new TLine(u[0],y,u[1],y);
   line->SetBit(kIsClassTree);
   line->Draw();
   Int_t icobject = FindClass("TObject");
   TPaveClass *label = new TPaveClass(xleft+gDxx,y-gLabdy,xleft+gLabdx,y+gLabdy,fCnames[iclass]->Data(),this);
   char *derived = fDerived[iclass];
   if (icobject >= 0 && !derived[icobject]) label->SetFillColor(30);
   if (fCstatus[iclass] > 1) label->SetFillColor(kYellow);
   label->SetTextSize(gCsize);
   label->SetBit(kIsClassTree);
   label->SetToolTipText(fCtitles[iclass]->Data(),500);
   label->Draw();
   if (ns == 0) return;

   // drawing sons
   y +=  0.5*gNtsons[iclass]*gDy;
   Int_t first =0;
   for (Int_t i=0;i<fNclasses;i++) {
      if(fCparent[i] != iclass) continue;
      if (gNtsons[i] > 1) y -= 0.5*gNtsons[i]*gDy;
      else               y -= 0.5*gDy;
      if (!first) {first=1; yu = y;}
      PaintClass(i,u[1],y);
      yl = y;
      if (gNtsons[i] > 1) y -= 0.5*gNtsons[i]*gDy;
      else               y -= 0.5*gDy;
   }
   if (ns == 1) return;
   line = new TLine(u[1],yl,u[1],yu);
   line->SetBit(kIsClassTree);
   line->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// save current configuration in a Root file
/// if filename is blank, the name of the file will be the current objectname.root
/// all the current settings are preserved
/// the Root file produced can be looked at by a another Root session
/// with no access to the original classes.
/// By default a message is printed. Specify option "Q" to remove the message

void TClassTree::SaveAs(const char *filename, Option_t *option) const
{
   if (gDirectory) gDirectory->SaveObjectAs(this,filename,option);
}

////////////////////////////////////////////////////////////////////////////////
///  Select all classes used by/referenced/referencing the class number iclass
///  and build the list of these classes

void TClassTree::ScanClasses(Int_t iclass)
{
   Int_t ic, icl;
   TList *los = fLinks[iclass];
   TList *losref = 0;
   TObjString *os;

   // scan list of data members
   // =========================
   TClass *cl = fCpointer[iclass];
   TDataMember *dm;
   TList *lm = cl->GetListOfDataMembers();
   if (lm) {
      TIter      next(lm);
      Int_t imember = 0;
      while ((dm = (TDataMember *) next())) {
         imember++;
         ic = FindClass(dm->GetTypeName());
         if (ic < 0 || ic == iclass) continue;
         losref = fLinks[ic];
         os = Mark(fCnames[ic]->Data(),los,kUsedByData);
         if (os) {
            os->SetBit(kIsaPointer,dm->IsaPointer());
            os->SetBit(kIsBasic,dm->IsBasic());
            os->SetUniqueID(imember);
         }
         Mark(fCnames[iclass]->Data(),losref,kUsingData);
      }
   }

   // scan base classes
   // =================
   char *derived = fDerived[iclass];
   TBaseClass *clbase;
   Int_t numb = 0;
   TList *lb = fCpointer[iclass]->GetListOfBases();
   if (lb) {
      TIter nextb(lb);
      while ((clbase = (TBaseClass*)nextb())) {
         numb++;
         if (numb == 1) continue;
         ic = FindClass(clbase->GetName());
         derived[ic] = 2;
      }
      for (ic=0;ic<fNclasses;ic++) {
         if (ic == iclass) continue;
         if (derived[ic]) {
            losref = fLinks[ic];
            Mark(fCnames[ic]->Data(),los,kUsedByClass);
            Mark(fCnames[iclass]->Data(),losref,kUsingClass);
         }
      }
   }

   // scan member functions
   // =====================
   char *star, *cref;
   TMethod *method;
   TMethodArg *methodarg;
   TList *lf = cl->GetListOfMethods();
   if (lf) {
      TIter nextm(lf);
      TString name;
      while ((method = (TMethod*) nextm())) {
         // check return type
         name = method->GetReturnTypeName();
         star = strstr((char*)name.Data(),"*");
         if (star) *star = 0;
         cref = strstr((char*)name.Data(),"&");
         if (cref) *cref = 0;
         ic = FindClass(name);
         if (ic < 0 || ic == iclass) continue;
         losref = fLinks[ic];
         Mark(fCnames[ic]->Data(),los,kUsedByFunc);
         Mark(fCnames[iclass]->Data(),losref,kUsingFunc);

         // now loop on all method arguments
         // ================================
         TIter nexta(method->GetListOfMethodArgs());
         while ((methodarg = (TMethodArg*) nexta())) {
            name = methodarg->GetTypeName();
            star = strstr((char*)name.Data(),"*");
            if (star) *star = 0;
            cref = strstr((char*)name.Data(),"&");
            if (cref) *cref = 0;
            ic = FindClass(name);
            if (ic < 0 || ic == iclass) continue;
            losref = fLinks[ic];
            Mark(fCnames[ic]->Data(),los,kUsedByFunc);
            Mark(fCnames[iclass]->Data(),losref,kUsingFunc);
         }
      }
   }

   // Look into the source code to search the list of includes
   // here we assume that include file names are classes file names
   // we stop reading the code when
   //   - a class member function is found
   //   - any class constructor is found
   if (!cl->GetImplFileName() || !cl->GetImplFileName()[0])
      return;

   const char *source = gSystem->BaseName( gSystem->UnixPathName(cl->GetImplFileName()));
   char *sourceName = gSystem->Which( fSourceDir.Data(), source , kReadPermission );
   if (!sourceName) return;
   Int_t ncn = strlen(fCnames[iclass]->Data())+2;
   char *cname = new char[ncn+1];
   snprintf(cname,ncn,"%s::",fCnames[iclass]->Data());
       // open source file
   std::ifstream sourceFile;
   sourceFile.open( sourceName, std::ios::in );
   Int_t nlines = 0;
   if( sourceFile.good() ) {
      const Int_t kMAXLEN=1500;
      char line[kMAXLEN];
      while( !sourceFile.eof() ) {
         sourceFile.getline( line, kMAXLEN-1 );
         if( sourceFile.eof() ) break;
         Int_t nblank = strspn(line," ");
         if (!strncmp(&line[nblank],"//",2)) continue;
         char *cc = strstr(line,"::");
         if (cc) {
            *cc = 0;
            if (!strncmp(&line[nblank],cname,ncn)) break;  //reach class member function
            Int_t nl = strlen(&line[nblank]);
            if (!strncmp(&line[nblank],cc+2,nl))   break;  //reach any class constructor
         }
         nlines++; if (nlines > 1000) break;
         char *inc = strstr(line,"#include");
         if (inc) {
            char *ch = strstr(line,".h");
            if (!ch) continue;
            *ch = 0;
            char *start = strstr(line,"<");
            if (!start) start = strstr(line,"\"");
            if (!start) continue;
            start++;
            while ((start < ch) && (*start == ' ')) start++;
            icl = FindClass(start);
            if (icl < 0 || icl == iclass) continue;
            // mark this include being used by this class
            losref = fLinks[icl];
            Mark(fCnames[icl]->Data(),los,kUsedByCode1);
            Mark(fCnames[icl]->Data(),los,kUsedByCode);
            Mark(fCnames[iclass]->Data(),losref,kUsingCode);
            // and also the base classes of the class in the include
            derived = fDerived[icl];
            for (ic=0;ic<fNclasses;ic++) {
               if (ic == icl) continue;
               if (derived[ic]) {
                  losref = fLinks[ic];
                  Mark(fCnames[ic]->Data(),los,kUsedByCode);
                  Mark(fCnames[iclass]->Data(),losref,kUsingCode);
               }
            }
         }
      }
   }
   delete [] cname;
   delete [] sourceName;
   sourceFile.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the list of classes for which the hierarchy is to be drawn
/// See Paint for the syntax

void TClassTree::SetClasses(const char *classes, Option_t *)
{
   if (classes == 0) return;
   fClasses = classes;
   for (Int_t i=0;i<fNclasses;i++) {
      fCstatus[i]  = 0;
      fCparent[i] = -1;
   }
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the size along x of the TPaveLabel showing the class name

void TClassTree::SetLabelDx(Float_t labeldx)
{
   fLabelDx = labeldx;
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the offset at the top of the picture
/// The default offset is computed automatically taking into account
/// classes not inheriting from TObject.

void TClassTree::SetYoffset(Float_t offset)
{
   fYoffset = offset;
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// mark classes used by the list of classes in classes

void TClassTree::ShowClassesUsedBy(const char *classes)
{
   Int_t i,j;
   Int_t nch = strlen(classes);
   char *ptr = new char[nch+1];
   strlcpy(ptr,classes,nch+1);
   if (ptr[0] == '*') {
      i = FindClass(&ptr[1]);
      if (i >= 0) {
         char *derived = fDerived[i];
         for (j=0;j<fNclasses;j++) {
            if(derived[j]) FindClassesUsedBy(j);
         }
      }
   } else if (ptr[nch-1] == '*') {
      ptr[nch-1] = 0;
      for (j=0;j<fNclasses;j++) {
         if(fCnames[j]->Contains(ptr)) FindClassesUsedBy(j);
      }
   } else {
      for (j=0;j<fNclasses;j++) {
         if(!fCnames[j]->CompareTo(ptr)) FindClassesUsedBy(j);
      }
   }
   delete [] ptr;
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// mark classes using any class in the list of classes in classes

void TClassTree::ShowClassesUsing(const char *classes)
{
   Int_t i,j;
   Int_t nch = strlen(classes);
   char *ptr = new char[nch+1];
   strlcpy(ptr,classes,nch+1);
   if (ptr[0] == '*') {
      i = FindClass(&ptr[1]);
      if (i >= 0) {
         char *derived = fDerived[i];
         for (j=0;j<fNclasses;j++) {
            if(derived[j]) FindClassesUsing(j);
         }
      }
   } else if (ptr[nch-1] == '*') {
      ptr[nch-1] = 0;
      for (j=0;j<fNclasses;j++) {
         if(fCnames[j]->Contains(ptr)) FindClassesUsing(j);
      }
   } else {
      for (j=0;j<fNclasses;j++) {
         if(!fCnames[j]->CompareTo(ptr)) FindClassesUsing(j);
      }
   }
   delete [] ptr;
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the Code References relationships

void TClassTree::ShowCod()
{
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   TObjString *os;
   TPaveClass *pave;
   Int_t ic,icl;
   Float_t x,y,x1,y1;
   //iterate on all TPaveClass objects in the pad
   while((obj=next())) {
      if (obj->InheritsFrom(TPaveClass::Class())) {
         pave = (TPaveClass*)obj;
         icl = FindClass(pave->GetLabel());
         if (icl < 0) continue;
         char *derived = fDerived[icl];
         x  = 0.5*(pave->GetX1() + pave->GetX2());
         y  = 0.5*(pave->GetY1() + pave->GetY2());
         TIter nextos(fLinks[icl]);
         //iterate on all classes in the list of classes of this class
         while((os=(TObjString*)nextos())) {
            if (!os->TestBit(kUsedByCode1)) continue;
            ic = FindClass(os->GetName());
            if (derived[ic]) continue;
            FindClassPosition(os->GetName(),x1,y1);
            if (x1 == 0 || y1 == 0) continue; //may be pointed class was not drawn
            TArrow *arrow = new TArrow(x,y,x1,y1,0.008,"|>");
            arrow->SetLineColor(kGreen);
            arrow->SetFillColor(kGreen);
            arrow->SetBit(kIsClassTree);
            arrow->Draw();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the "Has a" relationships

void TClassTree::ShowHas()
{
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   TObjString *os;
   TPaveClass *pave;
   Int_t icl;
   Float_t y,x1,y1,dx;
   //iterate on all TPaveClass objects in the pad
   while((obj=next())) {
      if (obj->InheritsFrom(TPaveClass::Class())) {
         pave = (TPaveClass*)obj;
         icl = FindClass(pave->GetLabel());
         if (icl < 0) continue;
         y  = 0.5*(pave->GetY1() + pave->GetY2());
         Int_t nmembers = fNdata[icl];
         if (nmembers == 0) continue;
         dx = (pave->GetX2() - pave->GetX1())/nmembers;
         TIter nextos(fLinks[icl]);
         //iterate on all classes in the list of classes of this class
         while((os=(TObjString*)nextos())) {
            if (!os->TestBit(kUsedByData)) continue;
            if (os->TestBit(kIsaPointer)) continue;
            if (os->TestBit(kIsBasic)) continue;
            FindClassPosition(os->GetName(),x1,y1);
            if (x1 == 0 || y1 == 0) continue; //may be base class was not drawn
            Int_t imember = os->GetUniqueID();
            TLine *line = new TLine(pave->GetX1()+(imember+0.5)*dx,y,x1,y1);
            line->SetLineStyle(3);
            line->SetLineColor(6);
            line->SetBit(kIsClassTree);
            line->Draw();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set link options in the ClassTree object
///
///  - "C"  show References from code
///  - "H"  show Has a relations
///  - "M"  show Multiple Inheritance
///  - "R"  show References from data members

void TClassTree::ShowLinks(Option_t *option)
{
   TString opt = option;
   opt.ToUpper();
   fShowCod = fShowHas = fShowMul = fShowRef = 0;
   if (opt.Contains("C")) fShowCod = 1;
   if (opt.Contains("H")) fShowHas = 1;
   if (opt.Contains("M")) fShowMul = 1;
   if (opt.Contains("R")) fShowRef = 1;
   if (gPad) Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the Multiple inheritance relationships

void TClassTree::ShowMul()
{
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   TObjString *os;
   TPaveClass *pave;
   Int_t ic,icl;
   Float_t x,y,x1,y1;
   //iterate on all TPaveClass objects in the pad
   while((obj=next())) {
      if (obj->InheritsFrom(TPaveClass::Class())) {
         pave = (TPaveClass*)obj;
         icl = FindClass(pave->GetLabel());
         if (icl < 0) continue;
         char *derived = fDerived[icl];
         x = 0.5*(pave->GetX1() + pave->GetX2());
         y = 0.5*(pave->GetY1() + pave->GetY2());
         TIter nextos(fLinks[icl]);
         //iterate on all classes in the list of classes of this class
         while((os=(TObjString*)nextos())) {
            if (!os->TestBit(kUsedByClass)) continue;
            ic = FindClass(os->GetName());
            if (derived[ic] != 2) continue; //keep only multiple inheritance
            FindClassPosition(os->GetName(),x1,y1);
            if (x1 == 0 || y1 == 0) continue; //may be base class was not drawn
            TLine *line = new TLine(x,y,x1,y1);
            line->SetBit(kIsClassTree);
            line->SetLineStyle(2);
            line->SetLineColor(kBlue);
            line->Draw();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the References relationships (other than inheritance or composition)

void TClassTree::ShowRef()
{
   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   TObjString *os;
   TPaveClass *pave;
   Int_t ic,icl;
   Float_t y,x1,y1,dx;
   Int_t icc = FindClass("TClass");
   //iterate on all TPaveClass objects in the pad
   while((obj=next())) {
      if (obj->InheritsFrom(TPaveClass::Class())) {
         pave = (TPaveClass*)obj;
         icl = FindClass(pave->GetLabel());
         if (icl < 0) continue;
         y  = 0.5*(pave->GetY1() + pave->GetY2());
         Int_t nmembers = fNdata[icl];
         if (nmembers == 0) continue;
         dx = (pave->GetX2() - pave->GetX1())/nmembers;
         TIter nextos(fLinks[icl]);
         //iterate on all classes in the list of classes of this class
         while((os=(TObjString*)nextos())) {
            if (!os->TestBit(kUsedByData)) continue;
            ic = FindClass(os->GetName());
            if (!os->TestBit(kIsaPointer)) continue;
            if (os->TestBit(kIsBasic)) continue;
            if (ic == icc) continue; // do not show relations with TClass
            FindClassPosition(os->GetName(),x1,y1);
            if (x1 == 0 || y1 == 0) continue; //may be pointed class was not drawn
            Int_t imember = os->GetUniqueID();
            TArrow *arrow = new TArrow(pave->GetX1()+(imember+0.5)*dx,y,x1,y1,0.008,"|>");
            arrow->SetLineColor(kRed);
            arrow->SetFillColor(kRed);
            arrow->SetBit(kIsClassTree);
            arrow->Draw();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TClassTree.
/// the status of the object is saved and can be replayed in a subsequent session

void TClassTree::Streamer(TBuffer &R__b)
{
   Int_t i;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      fClasses.Streamer(R__b);
      R__b >> fYoffset;
      R__b >> fLabelDx;
      R__b >> fNclasses;
      R__b >> fShowCod;
      R__b >> fShowMul;
      R__b >> fShowHas;
      R__b >> fShowRef;
      fCnames     = new TString*[fNclasses];
      fCtitles    = new TString*[fNclasses];
      fCstatus    = new Int_t[fNclasses];
      fParents    = new Int_t[fNclasses];
      fCparent    = new Int_t[fNclasses];
      fNdata      = new Int_t[fNclasses];
      fCpointer   = new TClass*[fNclasses];
      fOptions    = new TString*[fNclasses];
      fLinks      = new TList*[fNclasses];
      fDerived    = new char*[fNclasses];
      for (i=0;i<fNclasses;i++) {
         R__b >> fCstatus[i];
         R__b >> fParents[i];
         R__b >> fNdata[i];
         fCnames[i]  = new TString();
         fCtitles[i] = new TString();
         fOptions[i] = new TString();
         fCnames[i]->Streamer(R__b);
         fCtitles[i]->Streamer(R__b);
         fOptions[i]->Streamer(R__b);
         fLinks[i] = new TList();
         fLinks[i]->Streamer(R__b);
         fDerived[i] = new char[fNclasses];
         R__b.ReadFastArray(fDerived[i],fNclasses);
      }
      fSourceDir.Streamer(R__b);
   } else {
      R__b.WriteVersion(TClassTree::IsA());
      TNamed::Streamer(R__b);
      fClasses.Streamer(R__b);
      R__b << fYoffset;
      R__b << fLabelDx;
      R__b << fNclasses;
      R__b << fShowCod;
      R__b << fShowMul;
      R__b << fShowHas;
      R__b << fShowRef;
      for (i=0;i<fNclasses;i++) {
         R__b << fCstatus[i];
         R__b << fParents[i];
         R__b << fNdata[i];
         fCnames[i]->Streamer(R__b);
         fCtitles[i]->Streamer(R__b);
         fOptions[i]->Streamer(R__b);
         fLinks[i]->Streamer(R__b);
         R__b.WriteFastArray(fDerived[i],fNclasses);
      }
      fSourceDir.Streamer(R__b);
   }
}
