/// \file RFieldVisitor.cxx
/// \ingroup NTuple ROOT7
/// \author Simon Leisibach <simon.leisibach@gmail.com>
/// \date 2019-06-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RFieldVisitor.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


//----------------------------- RPrepareVisitor --------------------------------


void ROOT::Experimental::RPrepareVisitor::VisitField(const Detail::RFieldBase &field)
{
   auto subFields = field.GetSubFields();
   for (auto f : subFields) {
      RPrepareVisitor visitor;
      f->AcceptVisitor(visitor);
      fNumFields += visitor.fNumFields;
      fDeepestLevel = std::max(fDeepestLevel, 1 + visitor.fDeepestLevel);
   }
}


void ROOT::Experimental::RPrepareVisitor::VisitFieldZero(const RFieldZero &field)
{
   VisitField(field);
   fNumFields--;
   fDeepestLevel--;
}


//---------------------------- RPrintSchemaVisitor -----------------------------


void ROOT::Experimental::RPrintSchemaVisitor::SetDeepestLevel(int d)
{
   fDeepestLevel = d;
}

void ROOT::Experimental::RPrintSchemaVisitor::SetNumFields(int n)
{
   fNumFields = n;
   SetAvailableSpaceForStrings();
}

void ROOT::Experimental::RPrintSchemaVisitor::VisitField(const Detail::RFieldBase &field)
{
   fOutput << fFrameSymbol << ' ';

   std::string key = fTreePrefix;
   key += "Field " + fFieldNoPrefix + std::to_string(fFieldNo);
   fOutput << RNTupleFormatter::FitString(key, fAvailableSpaceKeyString);
   fOutput << " : ";

   std::string value = field.GetName() + " (" + field.GetType() + ")";
   fOutput << RNTupleFormatter::FitString(value, fAvailableSpaceValueString);
   fOutput << fFrameSymbol << std::endl;

   auto subFields = field.GetSubFields();
   auto fieldNo = 1;
   for (auto iField = subFields.begin(); iField != subFields.end(); ) {
      RPrintSchemaVisitor visitor(*this);
      visitor.fFieldNo = fieldNo++;
      visitor.fFieldNoPrefix += std::to_string(fFieldNo) + ".";

      auto f = *iField;
      ++iField;
      // TODO(jblomer): implement tree drawing
      visitor.fTreePrefix += "  ";
      f->AcceptVisitor(visitor);
   }
}


void ROOT::Experimental::RPrintSchemaVisitor::VisitFieldZero(const RFieldZero &fieldZero)
{
   auto fieldNo = 1;
   for (auto f : fieldZero.GetSubFields()) {
      RPrintSchemaVisitor visitor(*this);
      visitor.fFieldNo = fieldNo++;
      f->AcceptVisitor(visitor);
   }
}


//--------------------------- RPrintValueVisitor -------------------------------

void ROOT::Experimental::RPrintValueVisitor::PrintIndent()
{
   if (fPrintOptions.fPrintSingleLine)
      return;

   for (unsigned int i = 0; i < fLevel; ++i)
      fOutput << "  ";
}


void ROOT::Experimental::RPrintValueVisitor::PrintName(const Detail::RFieldBase &field)
{
   if (fPrintOptions.fPrintName)
      fOutput << "\"" << field.GetName() << "\": ";
}


void ROOT::Experimental::RPrintValueVisitor::PrintCollection(const Detail::RFieldBase &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << "[";
   auto elems = field.SplitValue(fValue);
   for (auto iValue = elems.begin(); iValue != elems.end(); ) {
      RPrintOptions options;
      options.fPrintSingleLine = true;
      options.fPrintName = false;
      RPrintValueVisitor elemVisitor(*iValue, fOutput, 0 /* level */, options);
      iValue->GetField()->AcceptVisitor(elemVisitor);

      if (++iValue == elems.end())
         break;
      else
         fOutput << ", ";
   }
   fOutput << "]";
}


void ROOT::Experimental::RPrintValueVisitor::VisitField(const Detail::RFieldBase &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << "\"<unsupported type: " << field.GetType() << ">\"";
}


void ROOT::Experimental::RPrintValueVisitor::VisitBoolField(const RField<bool> &field)
{
   PrintIndent();
   PrintName(field);
   if (*fValue.Get<bool>())
      fOutput << "true";
   else
      fOutput << "false";
}


void ROOT::Experimental::RPrintValueVisitor::VisitDoubleField(const RField<double> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << *fValue.Get<double>();
}


void ROOT::Experimental::RPrintValueVisitor::VisitFloatField(const RField<float> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << *fValue.Get<float>();
}


void ROOT::Experimental::RPrintValueVisitor::VisitIntField(const RField<int> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << *fValue.Get<int>();
}


void ROOT::Experimental::RPrintValueVisitor::VisitStringField(const RField<std::string> &field)
{
   PrintIndent();
   PrintName(field);
   // TODO(jblomer): escape double quotes
   fOutput << "\"" << *fValue.Get<std::string>() << "\"";
}


void ROOT::Experimental::RPrintValueVisitor::VisitUInt8Field(const RField<std::uint8_t> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << static_cast<int>(*fValue.Get<std::uint8_t>());
}


void ROOT::Experimental::RPrintValueVisitor::VisitUInt32Field(const RField<std::uint32_t> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << *fValue.Get<std::uint32_t>();
}


void ROOT::Experimental::RPrintValueVisitor::VisitUInt64Field(const RField<std::uint64_t> &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << *fValue.Get<std::uint64_t>();
}


void ROOT::Experimental::RPrintValueVisitor::VisitArrayField(const RArrayField &field)
{
   PrintCollection(field);
}


void ROOT::Experimental::RPrintValueVisitor::VisitClassField(const RClassField &field)
{
   PrintIndent();
   PrintName(field);
   fOutput << "{";
   auto elems = field.SplitValue(fValue);
   for (auto iValue = elems.begin(); iValue != elems.end(); ) {
      if (!fPrintOptions.fPrintSingleLine)
         fOutput << std::endl;

      RPrintOptions options;
      options.fPrintSingleLine = fPrintOptions.fPrintSingleLine;
      RPrintValueVisitor visitor(*iValue, fOutput, fLevel + 1, options);
      iValue->GetField()->AcceptVisitor(visitor);

      if (++iValue == elems.end()) {
         if (!fPrintOptions.fPrintSingleLine)
            fOutput << std::endl;
         break;
      } else {
         fOutput << ",";
         if (fPrintOptions.fPrintSingleLine)
           fOutput << " ";
      }
   }
   PrintIndent();
   fOutput << "}";
}


void ROOT::Experimental::RPrintValueVisitor::VisitVectorField(const RVectorField &field)
{
   PrintCollection(field);
}


void ROOT::Experimental::RPrintValueVisitor::VisitVectorBoolField(const RField<std::vector<bool>> &field)
{
   PrintCollection(field);
}


//---------------------------- RNTupleFormatter --------------------------------


std::string ROOT::Experimental::RNTupleFormatter::FitString(const std::string &str, int availableSpace)
{
   int strSize{static_cast<int>(str.size())};
   if (strSize <= availableSpace)
      return str + std::string(availableSpace - strSize, ' ');
   else if (availableSpace < 3)
      return std::string(availableSpace, '.');
   return std::string(str, 0, availableSpace - 3) + "...";
}
