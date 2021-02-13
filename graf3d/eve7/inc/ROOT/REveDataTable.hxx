// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT7_REveDataTable
#define ROOT7_REveDataTable

#include <ROOT/REveElement.hxx>

namespace ROOT {
namespace Experimental {

class REveDataCollection;

class REveDataTable : public REveElement
{
protected:
   const REveDataCollection *fCollection{nullptr};

public:
   REveDataTable(const std::string& n = "REveDataTable", const std::string& t = "");
   virtual ~REveDataTable() {}

   void SetCollection(const REveDataCollection *col) { fCollection = col; }
   const REveDataCollection *GetCollection() const { return fCollection; }

   void PrintTable();
   virtual Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset);

   void AddNewColumn(const std::string& expr, const std::string& title, int prec = 2);
};

//==============================================================================

class REveDataColumn : public REveElement
{
public:
   enum FieldType_e { FT_Double = 0, FT_Bool, FT_String };

protected:
public:
   TString fExpression;
   FieldType_e fType; // can we auto detect this?
   Int_t fPrecision{2};

   std::string fTrue{"*"};
   std::string fFalse{" "};

   std::function<double(void *)> fDoubleFoo;
   std::function<bool(void *)> fBoolFoo;
   std::function<std::string(void *)> fStringFoo;

public:
   REveDataColumn(const std::string& n = "REveDataColumn", const std::string& t = "");
   virtual ~REveDataColumn() {}

   void SetExpressionAndType(const std::string &expr, FieldType_e type);
   void SetExpressionAndType(const std::string &expr, FieldType_e type, TClass* c);
   void SetPrecision(Int_t prec);

   std::string EvalExpr(void *iptr) const;
};


} // namespace Experimental
} // namespace ROOT
#endif
