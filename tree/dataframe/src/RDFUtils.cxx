// Author: Enrico Guiraud, Danilo Piparo CERN  03/2017

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <stdexcept>
#include <string>
#include <typeinfo>

#include "RConfigure.h"      // R__USE_IMT
#include "ROOT/RDFNodes.hxx" // ColumnName2ColumnTypeName -> RCustomColumnBase, FindUnknownColumns -> RLoopManager
#include "ROOT/RDataSource.hxx"
#include "RtypesCore.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClassRef.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TROOT.h" // IsImplicitMTEnabled, GetImplicitMTPoolSize
#include "TTree.h"

using namespace ROOT::Detail::RDF;
using namespace ROOT::RDF;

namespace ROOT {
namespace Internal {
namespace RDF {

/// Return the type_info associated to a name. If the association fails, an
/// exception is thrown.
/// References and pointers are not supported since those cannot be stored in
/// columns.
const std::type_info &TypeName2TypeID(const std::string &name)
{
   if (auto c = TClass::GetClass(name.c_str())) {
      return *c->GetTypeInfo();
   } else if (name == "char" || name == "Char_t")
      return typeid(char);
   else if (name == "unsigned char" || name == "UChar_t")
      return typeid(unsigned char);
   else if (name == "int" || name == "Int_t")
      return typeid(int);
   else if (name == "unsigned int" || name == "UInt_t")
      return typeid(unsigned int);
   else if (name == "short" || name == "Short_t")
      return typeid(short);
   else if (name == "unsigned short" || name == "UShort_t")
      return typeid(unsigned short);
   else if (name == "long" || name == "Long_t")
      return typeid(long);
   else if (name == "unsigned long" || name == "ULong_t")
      return typeid(unsigned long);
   else if (name == "double" || name == "Double_t")
      return typeid(double);
   else if (name == "float" || name == "Float_t")
      return typeid(float);
   else if (name == "long long" || name == "long long int" || name == "Long64_t")
      return typeid(Long64_t);
   else if (name == "unsigned long long" || name == "unsigned long long int" || name == "ULong64_t")
      return typeid(ULong64_t);
   else if (name == "bool" || name == "Bool_t")
      return typeid(bool);
   else {
      std::string msg("Cannot extract type_info of type ");
      msg += name.c_str();
      msg += ".";
      throw std::runtime_error(msg);
   }
}

/// Returns the name of a type starting from its type_info
/// An empty string is returned in case of failure
/// References and pointers are not supported since those cannot be stored in
/// columns.
std::string TypeID2TypeName(const std::type_info &id)
{
   if (auto c = TClass::GetClass(id)) {
      return c->GetName();
   } else if (id == typeid(char))
      return "char";
   else if (id == typeid(unsigned char))
      return "unsigned char";
   else if (id == typeid(int))
      return "int";
   else if (id == typeid(unsigned int))
      return "unsigned int";
   else if (id == typeid(short))
      return "short";
   else if (id == typeid(unsigned short))
      return "unsigned short";
   else if (id == typeid(long))
      return "long";
   else if (id == typeid(unsigned long))
      return "unsigned long";
   else if (id == typeid(double))
      return "double";
   else if (id == typeid(float))
      return "float";
   else if (id == typeid(Long64_t))
      return "Long64_t";
   else if (id == typeid(ULong64_t))
      return "ULong64_t";
   else if (id == typeid(bool))
      return "bool";
   else
      return "";
}

std::string ComposeRVecTypeName(const std::string &valueType)
{
   return "ROOT::VecOps::RVec<" + valueType + ">";
}

std::string GetLeafTypeName(TLeaf *leaf, const std::string &colName)
{
   std::string colType = leaf->GetTypeName();
   if (colType.empty())
      throw std::runtime_error("Could not deduce type of leaf " + colName);
   if (leaf->GetLeafCount() != nullptr && leaf->GetLenStatic() == 1) {
      // this is a variable-sized array
      colType = ComposeRVecTypeName(colType);
   } else if (leaf->GetLeafCount() == nullptr && leaf->GetLenStatic() > 1) {
      // this is a fixed-sized array (we do not differentiate between variable- and fixed-sized arrays)
      colType = ComposeRVecTypeName(colType);
   } else if (leaf->GetLeafCount() != nullptr && leaf->GetLenStatic() > 1) {
      // we do not know how to deal with this branch
      throw std::runtime_error("TTree leaf " + colName +
                               " has both a leaf count and a static length. This is not supported.");
   }

   return colType;
}

/// Return the typename of object colName stored in t, if any. Return an empty string if colName is not in t.
/// Supported cases:
/// - leaves corresponding to single values, variable- and fixed-length arrays, with following syntax:
///   - "leafname", as long as TTree::GetLeaf resolves it
///   - "b1.b2...leafname", as long as TTree::GetLeaf("b1.b2....", "leafname") resolves it
/// - TBranchElements, as long as TTree::GetBranch resolves their names
std::string GetBranchOrLeafTypeName(TTree &t, const std::string &colName)
{
   // look for TLeaf either with GetLeaf(colName) or with GetLeaf(branchName, leafName) (splitting on last dot)
   auto leaf = t.GetLeaf(colName.c_str());
   if (!leaf) {
      const auto dotPos = colName.find_last_of('.');
      const auto hasDot = dotPos != std::string::npos;
      if (hasDot) {
         const auto branchName = colName.substr(0, dotPos);
         const auto leafName = colName.substr(dotPos + 1);
         leaf = t.GetLeaf(branchName.c_str(), leafName.c_str());
      }
   }
   if (leaf)
      return GetLeafTypeName(leaf, colName);

   // we could not find a leaf named colName, so we look for a TBranchElement
   auto branch = t.GetBranch(colName.c_str());
   if (branch) {
      static const TClassRef tbranchelement("TBranchElement");
      if (branch->InheritsFrom(tbranchelement)) {
         auto be = static_cast<TBranchElement *>(branch);
         if (auto currentClass = be->GetCurrentClass())
            return currentClass->GetName();
         else
            return be->GetClassName();
      }
   }

   // colName is not a leaf nor a TBranchElement
   return std::string();
}

/// Return a string containing the type of the given branch. Works both with real TTree branches and with temporary
/// column created by Define. Throws if type name deduction fails.
/// Note that for fixed- or variable-sized c-style arrays the returned type name will be RVec<T>.
/// vector2tvec specifies whether typename 'std::vector<T>' should be converted to 'RVec<T>' or returned as is
std::string ColumnName2ColumnTypeName(const std::string &colName, unsigned int namespaceID, TTree *tree,
                                      RDataSource *ds, bool isCustomColumn, bool vector2tvec)
{
   std::string colType;

   if (ds && ds->HasColumn(colName))
      colType = ds->GetTypeName(colName);

   if (colType.empty() && tree) {
      colType = GetBranchOrLeafTypeName(*tree, colName);
      if (vector2tvec && TClassEdit::IsSTLCont(colType) == ROOT::ESTLType::kSTLvector) {
         std::vector<std::string> split;
         int dummy;
         TClassEdit::GetSplit(colType.c_str(), split, dummy);
         auto &valueType = split[1];
         colType = ComposeRVecTypeName(valueType);
      }
   }

   if (colType.empty() && isCustomColumn) {
      // this must be a temporary branch, we know there is an alias for its type
      colType = "__tdf" + std::to_string(namespaceID) + "::" + colName + "_type";
   }

   if (colType.empty())
      throw std::runtime_error("Column \"" + colName +
                               "\" is not in a dataset and is not a custom column been defined.");

   return colType;
}

/// Convert type name (e.g. "Float_t") to ROOT type code (e.g. 'F') -- see TBranch documentation.
/// Return a space ' ' in case no match was found.
char TypeName2ROOTTypeName(const std::string &b)
{
   if (b == "Char_t" || b == "char")
      return 'B';
   if (b == "UChar_t" || b == "unsigned char")
      return 'b';
   if (b == "Short_t" || b == "short" || b == "short int")
      return 'S';
   if (b == "UShort_t" || b == "unsigned short" || b == "unsigned short int")
      return 's';
   if (b == "Int_t" || b == "int")
      return 'I';
   if (b == "UInt_t" || b == "unsigned" || b == "unsigned int")
      return 'i';
   if (b == "Float_t" || b == "float")
      return 'F';
   if (b == "Double_t" || b == "double")
      return 'D';
   if (b == "Long64_t" || b == "long" || b == "long int")
      return 'L';
   if (b == "ULong64_t" || b == "unsigned long" || b == "unsigned long int")
      return 'l';
   if (b == "Bool_t" || b == "bool")
      return 'O';
   return ' ';
}

unsigned int GetNSlots()
{
   unsigned int nSlots = 1;
#ifdef R__USE_IMT
   if (ROOT::IsImplicitMTEnabled())
      nSlots = ROOT::GetImplicitMTPoolSize();
#endif // R__USE_IMT
   return nSlots;
}

/// Replace occurrences of '.' with '_' in each string passed as argument.
/// An Info message is printed when this happens. Dots at the end of the string are not replaced.
/// An exception is thrown in case the resulting set of strings would contain duplicates.
std::vector<std::string> ReplaceDotWithUnderscore(const std::vector<std::string> &columnNames)
{
   auto newColNames = columnNames;
   for (auto &col : newColNames) {
      const auto dotPos = col.find('.');
      if (dotPos != std::string::npos && dotPos != col.size() - 1 && dotPos != 0u) {
         auto oldName = col;
         std::replace(col.begin(), col.end(), '.', '_');
         if (std::find(columnNames.begin(), columnNames.end(), col) != columnNames.end())
            throw std::runtime_error("Column " + oldName + " would be written as " + col +
                                     " but this column already exists. Please use Alias to select a new name for " +
                                     oldName);
         Info("Snapshot", "Column %s will be saved as %s", oldName.c_str(), col.c_str());
      }
   }

   return newColNames;
}

} // end NS RDF
} // end NS Internal
} // end NS ROOT
