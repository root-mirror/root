// @(#)root/graf:$Id$
// Author: Matthew.Adam.Dobbs   06/09/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <stdio.h>

#include "TLegendEntry.h"
#include "TVirtualPad.h"
#include "TROOT.h"
#include <iostream>

ClassImp(TLegendEntry);

/** \class TLegendEntry
\ingroup BasicGraphics

Storage class for one entry of a TLegend.
*/

////////////////////////////////////////////////////////////////////////////////
/// TLegendEntry do-nothing default constructor

TLegendEntry::TLegendEntry(): TAttText(), TAttLine(), TAttFill(), TAttMarker()
{
   fObject = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TLegendEntry normal constructor for one entry in a TLegend.
///
/// obj is the object this entry will represent. If obj has
/// line/fill/marker attributes, then the TLegendEntry will display
/// these attributes.
///
/// label is the text that will describe the entry, it is displayed using
/// TLatex, so may have a complex format.
///
/// option may have values
///  - L draw line associated w/ TAttLine if obj inherits from TAttLine
///  - P draw polymarker assoc. w/ TAttMarker if obj inherits from TAttMarker
///  - F draw a box with fill associated w/ TAttFill if obj inherits TAttFill
///    default is object = "LPF"

TLegendEntry::TLegendEntry(const TObject* obj, const char* label, Option_t* option )
             :TAttText(0,0,0,0,0), TAttLine(1,1,1), TAttFill(0,0), TAttMarker(1,21,1)
{
   fObject = 0;
   if ( !label && obj ) fLabel = obj->GetTitle();
   else                 fLabel = label;
   fOption = option;
   if (obj) SetObject((TObject*)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// TLegendEntry copy constructor

TLegendEntry::TLegendEntry( const TLegendEntry &entry ) : TObject(entry), TAttText(entry), TAttLine(entry), TAttFill(entry), TAttMarker(entry)
{
   ((TLegendEntry&)entry).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// TLegendEntry default destructor

TLegendEntry::~TLegendEntry()
{
   fObject = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// copy this TLegendEntry into obj

void TLegendEntry::Copy( TObject &obj ) const
{
   TObject::Copy(obj);
   TAttText::Copy((TLegendEntry&)obj);
   TAttLine::Copy((TLegendEntry&)obj);
   TAttFill::Copy((TLegendEntry&)obj);
   TAttMarker::Copy((TLegendEntry&)obj);
   ((TLegendEntry&)obj).fObject = fObject;
   ((TLegendEntry&)obj).fLabel = fLabel;
   ((TLegendEntry&)obj).fOption = fOption;
}

////////////////////////////////////////////////////////////////////////////////
/// dump this TLegendEntry to std::cout

void TLegendEntry::Print( Option_t *) const
{
   TString output;
   std::cout << "TLegendEntry: Object ";
   if ( fObject ) output = fObject->GetName();
   else output = "NULL";
   std::cout << output << " Label ";
   if ( fLabel ) output = fLabel.Data();
   else output = "NULL";
   std::cout << output << " Option ";
   if (fOption ) output = fOption.Data();
   else output = "NULL";
   std::cout << output << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save this TLegendEntry as C++ statements on output stream out
///  to be used with the SaveAs .C option

void TLegendEntry::SaveEntry(std::ostream &out, const char* name )
{
   char quote = '"';
   if ( gROOT->ClassSaved( TLegendEntry::Class() ) ) {
      out << "   entry=";
   } else {
      out << "   TLegendEntry *entry=";
   }
   TString objname = "NULL";
   if ( fObject ) objname = fObject->GetName();
   TString tL(fLabel);
   tL.ReplaceAll("\\","\\\\");
   tL.ReplaceAll("\"","\\\"");
   out << name << "->AddEntry("<<quote<<objname<<quote<<","<<quote<<
      tL.Data()<<quote<<","<<quote<<fOption.Data()<<quote<<");"<<std::endl;
   SaveFillAttributes(out,"entry",0,0);
   SaveLineAttributes(out,"entry",0,0,0);
   SaveMarkerAttributes(out,"entry",0,0,0);
   SaveTextAttributes(out,"entry",0,0,0,0,0);
}

////////////////////////////////////////////////////////////////////////////////
/// (re)set the obj pointed to by this entry

void TLegendEntry::SetObject(TObject* obj )
{
   if ( ( fObject && fLabel == fObject->GetTitle() ) || !fLabel ) {
      if (obj) fLabel = obj->GetTitle();
   }
   fObject = obj;
}

////////////////////////////////////////////////////////////////////////////////
/// (re)set the obj pointed to by this entry

void TLegendEntry::SetObject( const char* objectName)
{
   TObject* obj = 0;
   TList* padprimitives = gPad->GetListOfPrimitives();
   if (padprimitives) obj = padprimitives->FindObject( objectName );
   if (obj) SetObject( obj );
}
