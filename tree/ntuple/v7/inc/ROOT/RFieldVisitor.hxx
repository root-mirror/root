/// \file ROOT/RFieldVisitor.hxx
/// \ingroup NTuple ROOT7
/// \author Simon Leisibach <simon.satoshi.rene.leisibach@cern.ch>
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

#ifndef ROOT7_RFieldVisitor
#define ROOT7_RFieldVisitor

#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RNTupleUtil.hxx>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {
namespace Detail {


// clang-format off
/**
\class ROOT::Experimental::Detail::RFieldVisitor
\ingroup NTuple
\brief Abstract base class for classes implementing the visitor design pattern.

RFieldVisitor::VisitField() is invoked by RFieldBase::AcceptVisitor(). VisitField() is inherited for instance
by the RPrintSchemaVisitor class. The RFieldBase class and classes which inherit from it will be visited.
*/
// clang-format on
class RFieldVisitor {
public:
   virtual void VisitField(const Detail::RFieldBase &field) = 0;
   virtual void VisitFieldZero(const RFieldZero &field) { VisitField(field); }
   virtual void VisitArrayField(const RArrayField &field) { VisitField(field); }
   virtual void VisitBoolField(const RField<bool> &field) { VisitField(field); }
   virtual void VisitClassField(const RClassField &field) { VisitField(field); }
   virtual void VisitClusterSizeField(const RField<ClusterSize_t> &field) { VisitField(field); }
   virtual void VisitDoubleField(const RField<double> &field) { VisitField(field); }
   virtual void VisitFloatField(const RField<float> &field) { VisitField(field); }
   virtual void VisitIntField(const RField<int> &field) { VisitField(field); }
   virtual void VisitStringField(const RField<std::string> &field) { VisitField(field); }
   virtual void VisitUInt32Field(const RField<std::uint32_t> &field) { VisitField(field); }
   virtual void VisitUInt64Field(const RField<std::uint64_t> &field) { VisitField(field); }
   virtual void VisitUInt8Field(const RField<std::uint8_t> &field) { VisitField(field); }
   virtual void VisitVectorField(const RVectorField &field) { VisitField(field); }
   virtual void VisitVectorBoolField(const RField<std::vector<bool>> &field) { VisitField(field); }
}; // class RFieldVisitor

} // namespace Detail


// clang-format off
/**
\class ROOT::Experimental::RPrepareVisitor
\ingroup NTuple
\brief Visitor used for a pre-processing run to collect information needed by another visitor class.

 Currently used for RPrintSchemaVisitor in RNTupleReader::Print() to collect information about levels, max depth etc.
*/
// clang-format on
class RPrepareVisitor : public Detail::RFieldVisitor {
private:
   unsigned int fDeepestLevel = 1;
   unsigned int fNumFields = 1;

public:
   RPrepareVisitor() = default;
   void VisitField(const Detail::RFieldBase &field) final;
   void VisitFieldZero(const RFieldZero &field) final;

   unsigned int GetDeepestLevel() const { return fDeepestLevel; }
   unsigned int GetNumFields() const { return fNumFields; }
};


// clang-format off
/**
\class ROOT::Experimental::RPrintSchemaVisitor
\ingroup NTuple
\brief Contains settings for printing and prints a summary of an RField instance.

This visitor is used by RNTupleReader::Print()
*/
// clang-format on
class RPrintSchemaVisitor : public Detail::RFieldVisitor {
private:
   /// Where to write the printout to
   std::ostream &fOutput;
   /// To render the output, use an asterix (*) by default to draw table lines and boundaries
   char fFrameSymbol;
   /// Indicates maximal number of allowed characters per line
   int fWidth;
   int fDeepestLevel;
   int fNumFields;
   int fAvailableSpaceKeyString;
   int fAvailableSpaceValueString;
   int fFieldNo = 1;
   std::string fTreePrefix;
   std::string fFieldNoPrefix;

public:
   RPrintSchemaVisitor(std::ostream &out = std::cout, char frameSymbol = '*', int width = 80, int deepestLevel = 1,
                       int numFields = 1)
      : fOutput{out}, fFrameSymbol{frameSymbol}, fWidth{width}, fDeepestLevel{deepestLevel}, fNumFields{numFields}
   {
      SetAvailableSpaceForStrings();
   }
   /// Prints summary of Field
   void VisitField(const Detail::RFieldBase &field) final;
   void VisitFieldZero(const RFieldZero &fieldZero) final;
   void SetFrameSymbol(char s) { fFrameSymbol = s; }
   void SetWidth(int w) { fWidth = w; }
   void SetDeepestLevel(int d);
   void SetNumFields(int n);
   /// Computes how many characters should be placed between the frame symbol and ':' for left and right side of ':' for
   /// visually pleasing output.
   // E.g.
   // * Field 1       : vpx (std::vector<float>)                                     *
   // * |__Field 1.1  : vpx/vpx (float)                                              *
   // int fAvailableSpaceKeyString (num characters on left side between "* " and " : ")
   //    deepestlevel here is 2 (1.1 is deepest and has 2 numbers).
   //    For every additional level an additional "|_" and ".1" (4 characters) is added, so the number of required
   //    additional characters is 4 * fDeepestLevel. For level 1, 8 characters are required ("Field 1 "), so an
   //    additional + 4 is added. To account for cases where the total number of fields is not a single digit number and
   //    more space is required to output big numbers, fNumFields is incorporated into the calculation. To make sure
   //    that there is still enough space for the right side of " : ", an std::min comparision with fWidth - 15 is done.
   // int fAvailableSpaceValueString(num characters on right side between " : " and '*')
   //    The 6 subtracted characters are "* " (2) in the beginning,  " : " (3) and '*' (1) on the far right.
   void SetAvailableSpaceForStrings()
   {
      fAvailableSpaceKeyString =
         std::min(4 * fDeepestLevel + 4 + static_cast<int>(std::to_string(fNumFields).size()), fWidth - 15);
      fAvailableSpaceValueString = fWidth - 6 - fAvailableSpaceKeyString;
   }
};


// clang-format off
/**
\class ROOT::Experimental::RPrintValueVisitor
\ingroup NTuple
\brief Renders a JSON value corresponding to the field.
*/
// clang-format on
class RPrintValueVisitor : public Detail::RFieldVisitor {
public:
   struct RPrintOptions {
      bool fPrintSingleLine;
      bool fPrintName;

      RPrintOptions() : fPrintSingleLine(false), fPrintName(true) {}
   };

private:
   Detail::RFieldValue fValue;
   /// The output is directed to fOutput which may differ from std::cout.
   std::ostream &fOutput;
   unsigned int fLevel;
   RPrintOptions fPrintOptions;

   void PrintIndent();
   void PrintName(const Detail::RFieldBase &field);
   void PrintCollection(const Detail::RFieldBase &field);

public:
   RPrintValueVisitor(const Detail::RFieldValue &value,
                      std::ostream &output,
                      unsigned int level = 0,
                      RPrintOptions options = RPrintOptions())
      : fValue(value), fOutput{output}, fLevel(level), fPrintOptions(options) {}

   void VisitField(const Detail::RFieldBase &field) final;

   void VisitBoolField(const RField<bool> &field) final;
   void VisitDoubleField(const RField<double> &field) final;
   void VisitFloatField(const RField<float> &field) final;
   void VisitIntField(const RField<int> &field) final;
   void VisitStringField(const RField<std::string> &field) final;
   void VisitUInt8Field(const RField<std::uint8_t> &field) final;
   void VisitUInt32Field(const RField<std::uint32_t> &field) final;
   void VisitUInt64Field(const RField<std::uint64_t> &field) final;

   void VisitArrayField(const RArrayField &field) final;
   void VisitClassField(const RClassField &field) final;
   void VisitVectorField(const RVectorField &field) final;
   void VisitVectorBoolField(const RField<std::vector<bool>> &field) final;
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleFormatter
\ingroup NTuple
\brief Contains helper functions for RNTupleReader::PrintInfo() and RPrintSchemaVisitor::VisitField()

The functions in this class format strings which are displayed by RNTupleReader::PrintInfo() and RNTupleReader::Show().
*/
// clang-format on
class RNTupleFormatter {
public:
   // Can abbreviate long strings, e.g. ("ExampleString" , space= 8) => "Examp..."
   static std::string FitString(const std::string &str, int availableSpace);
};

} // namespace Experimental
} // namespace ROOT

#endif
