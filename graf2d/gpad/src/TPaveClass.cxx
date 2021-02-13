// @(#)root/gpad:$Id$
// Author: Rene Brun   06/08/99
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "Buttons.h"
#include "TPaveClass.h"
#include "TClassTree.h"

#include <iostream>

ClassImp(TPaveClass);


/** \class TPaveClass
\ingroup gpad

A TPaveLabel specialized to process classes inside a TClassTree.
A TPaveClass object is used by the TClassTree to represent a class.
A TPaveClass has the same graphical representation as a TPaveLabel.

Using the context menu on can select additional options in the ClassTree:
  - Show classes using this class
  - Show all classes used by this class
*/

////////////////////////////////////////////////////////////////////////////////
/// PaveClass default constructor.

TPaveClass::TPaveClass(): TPaveLabel()
{
   fClassTree  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// PaveClass normal constructor.

TPaveClass::TPaveClass(Double_t x1, Double_t y1,Double_t x2, Double_t  y2, const char *label, TClassTree *classtree)
           :TPaveLabel(x1,y1,x2,y2,label,"br")
{
   fClassTree  = classtree;
   SetName(label);
   SetTextFont(61);
}

////////////////////////////////////////////////////////////////////////////////
/// PaveClass default destructor.

TPaveClass::~TPaveClass()
{
}

////////////////////////////////////////////////////////////////////////////////
/// PaveClass copy constructor.

TPaveClass::TPaveClass(const TPaveClass &PaveClass) : TPaveLabel(PaveClass)
{
   ((TPaveClass&)PaveClass).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this PaveClass to PaveClass.

void TPaveClass::Copy(TObject &obj) const
{
   TPaveLabel::Copy(obj);
   ((TPaveClass&)obj).fClassTree      = fClassTree;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw classes.

void TPaveClass::DrawClasses(const char *classes)
{
   if (!fClassTree) return;
   if (!strcmp(classes,"this")) fClassTree->Draw(GetName());
   else                         fClassTree->Draw(classes);
}

////////////////////////////////////////////////////////////////////////////////
/// Save as.

void TPaveClass::SaveAs(const char *filename, Option_t *option) const
{
   if (!fClassTree) return;
   fClassTree->SaveAs(filename,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set classes.

void TPaveClass::SetClasses(const char *classes, Option_t *option)
{
   if (!fClassTree) return;
   if (!strcmp(classes,"this")) fClassTree->SetClasses(GetName(),option);
   else                         fClassTree->SetClasses(classes,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set link options in the ClassTree object.
///
///  - "C"  show References from code
///  - "H"  show "Has a" relations
///  - "M"  show Multiple Inheritance
///  - "R"  show References from data members

void TPaveClass::ShowLinks(Option_t *option)
{
   if (!fClassTree) return;
   fClassTree->ShowLinks(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Show classes used by.

void TPaveClass::ShowClassesUsedBy(const char *classes)
{
   if (!fClassTree) return;
   if (!strcmp(classes,"this")) fClassTree->ShowClassesUsedBy(GetName());
   else                         fClassTree->ShowClassesUsedBy(classes);
}

////////////////////////////////////////////////////////////////////////////////
/// Show classes using.

void TPaveClass::ShowClassesUsing(const char *classes)
{
   if (!fClassTree) return;
   if (!strcmp(classes,"this")) fClassTree->ShowClassesUsing(GetName());
   else                         fClassTree->ShowClassesUsing(classes);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TPaveClass::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TPaveClass::Class())) {
      out<<"   ";
   } else {
      out<<"   TPaveClass *";
   }
   out<<"pclass = new TPaveClass("<<fX1<<","<<fY1<<","<<fX2<<","<<fY2
      <<","<<quote<<fLabel<<quote<<","<<quote<<fOption<<quote<<");"<<std::endl;

   SaveFillAttributes(out,"pclass",0,1001);
   SaveLineAttributes(out,"pclass",1,1,1);
   SaveTextAttributes(out,"pclass",22,0,1,62,0);

   out<<"   pclass->Draw();"<<std::endl;
}
