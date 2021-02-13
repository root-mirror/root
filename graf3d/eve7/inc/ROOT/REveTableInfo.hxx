// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveTableInfo
#define ROOT7_REveTableInfo

#include <ROOT/REveElement.hxx>
#include <ROOT/REveDataCollection.hxx>
#include <ROOT/REveDataTable.hxx>

namespace ROOT {
namespace Experimental {

///////////////////////////////////////////////////////////////////////////////
/// REveTableEntry
///////////////////////////////////////////////////////////////////////////////

class REveTableEntry {
public:
   std::string    fName;
   int            fPrecision;
   std::string    fExpression;
   REveDataColumn::FieldType_e fType;

   REveTableEntry() : fName("unknown"), fPrecision(2), fType(REveDataColumn::FT_Double) {}

   REveTableEntry(const std::string &name, int precision, const std::string &expression)
      : fName(name), fPrecision(precision), fExpression(expression), fType(REveDataColumn::FT_Double)
   {
   }

   void Print() const
   {
      printf("TableEntry\n");
      printf("name: %s expression: %s\n", fName.c_str(), fExpression.c_str());
   }
};

///////////////////////////////////////////////////////////////////////////////
/// REveTableHandle
///////////////////////////////////////////////////////////////////////////////

class REveTableHandle
{
   friend class REveTableViewInfo;

public:
   typedef std::vector<REveTableEntry> Entries_t;
   typedef std::map<std::string, Entries_t> Specs_t;

   // REveTableHandle() {}

   REveTableHandle&
   column(const std::string &name, int precision, const std::string &expression)
   {
      fSpecs[fClassName].emplace_back(name, precision, expression);
      return *this;
   }

   REveTableHandle &column(const std::string &label, int precision)
   {
      return column(label, precision, label);
   }

   REveTableHandle(std::string className, Specs_t &specs)
      :fClassName(className), fSpecs(specs)
   {
   }

protected:
   std::string  fClassName;
   Specs_t&  fSpecs;
};

///////////////////////////////////////////////////////////////////////////////
/// REveTableViewInfo
///////////////////////////////////////////////////////////////////////////////

class REveTableViewInfo : public REveElement
{
public:
   REveTableViewInfo(const std::string &name = "TableViewManager", const std::string &title = "");

   typedef std::function<void ()> Delegate_t;

   void SetDisplayedCollection(ElementId_t);
   ElementId_t GetDisplayedCollection() const  { return fDisplayedCollection; }

   void AddNewColumnToCurrentCollection(const std::string& expr, const std::string& title, int prec = 2);

   void AddDelegate(Delegate_t d) { fDelegates.push_back(d); }

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   // read
   REveTableHandle::Entries_t &RefTableEntries(std::string cname);

   // filling
   REveTableHandle table(std::string className)
   {
      REveTableHandle handle(className, fSpecs);
      return handle;
   }

   bool GetConfigChanged() const { return fConfigChanged; }


private:
   int fDisplayedCollection{0};
   std::vector<Delegate_t> fDelegates;
   REveTableHandle::Specs_t  fSpecs;
   bool                      fConfigChanged{false};
};


}
}

#endif
