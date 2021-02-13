/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RMenuItems
#define ROOT7_RMenuItems

#include <string>
#include <vector>
#include <memory>

#include <ROOT/RDrawableRequest.hxx>

#include "TClass.h"

namespace ROOT {
namespace Experimental {
namespace Detail {

/** \class RMenuItem
\ingroup GpadROOT7
\brief Base class for menu items, shown on JS side.
\author Sergey Linev
\date 2017-06-29
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RMenuItem {
protected:
   std::string fName;       ///< name of the menu item
   std::string fTitle;      ///< title of menu item
   std::string fExec;       ///< execute when item is activated
   std::string fClassName;  ///< class to which belongs menu item
public:
   /** Default constructor */
   RMenuItem() = default;

   /** Create menu item with the name and title
    *  name used to display item in the object context menu,
    *  title shown as hint info for that item  */
   RMenuItem(const std::string &name, const std::string &title) : fName(name), fTitle(title) {}

   /** virtual destructor need for vtable, used when vector of RMenuItem* is stored */
   virtual ~RMenuItem() = default;

   /** Set execution string with all required arguments,
    * which will be executed when menu item is selected  */
   void SetExec(const std::string &exec) { fExec = exec; }

   /** Set class name for menu item, will be used to group items together */
   void SetClassName(const std::string &clname) { fClassName = clname; }

   /** Returns menu item name */
   const std::string &GetName() const { return fName; }

   /** Returns execution string for the menu item */
   const std::string &GetExec() const { return fExec; }
};

/** \class RCheckedMenuItem
\ingroup GpadROOT7
\brief Menu item with check box
\author Sergey Linev
\date 2017-06-29
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RCheckedMenuItem : public RMenuItem {
protected:
   bool fChecked = false; ///< state of checkbox
public:
   /** Default constructor */
   RCheckedMenuItem() = default;

   /** Create checked menu item  */
   RCheckedMenuItem(const std::string &name, const std::string &title, bool checked = false)
      : RMenuItem(name, title), fChecked(checked)
   {
   }

   /** virtual destructor need for vtable, used when vector of RMenuItem* is stored */
   virtual ~RCheckedMenuItem() {}

   /** Set checked state for the item, default is none */
   void SetChecked(bool on = true) { fChecked = on; }

   bool IsChecked() const { return fChecked; }
};

/** \class RMenuArgument
\ingroup GpadROOT7
\brief Argument description for menu item which should invoke class method
\author Sergey Linev
\date 2017-06-29
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RMenuArgument {
protected:
   std::string fName;     ///<  name of call argument
   std::string fTitle;    ///<  title of call argument
   std::string fTypeName; ///<  typename
   std::string fDefault;  ///<  default value
public:
   /** Default constructor */
   RMenuArgument() = default;

   RMenuArgument(const std::string &name, const std::string &title, const std::string &typname,
                 const std::string &dflt = "")
      : fName(name), fTitle(title), fTypeName(typname), fDefault(dflt)
   {
   }

   void SetDefault(const std::string &dflt) { fDefault = dflt; }
};

/** \class RArgsMenuItem
\ingroup GpadROOT7
\brief Menu item which requires extra arguments for invoked class method
\author Sergey Linev
\date 2017-06-29
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RArgsMenuItem : public RMenuItem {
protected:
   std::vector<RMenuArgument> fArgs;

public:
   /** Default constructor */
   RArgsMenuItem() = default;

   RArgsMenuItem(const std::string &name, const std::string &title) : RMenuItem(name, title) {}

   /** virtual destructor need for vtable, used when vector of RMenuItem* is stored */
   virtual ~RArgsMenuItem() {}

   void AddArg(const RMenuArgument &arg) { fArgs.emplace_back(arg); }
};

} // namespace Detail

///////////////////////////////////////////////////////////////////////

/** \class RMenuItems
\ingroup GpadROOT7
\brief List of items for object context menu
\author Sergey Linev
\date 2017-06-29
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RMenuItems : public RDrawableReply {
protected:
   std::string fId;                                        ///< object identifier
   std::string fSpecifier;                                 ///<! extra specifier, used only on server
   std::vector<std::unique_ptr<Detail::RMenuItem>> fItems; ///< list of items in the menu
public:
   RMenuItems() = default;

   RMenuItems(const std::string &_id, const std::string &_specifier)
   {
      fId = _id;
      fSpecifier = _specifier;
   }

   virtual ~RMenuItems();

   const std::string &GetFullId() const { return fId; }
   const std::string &GetSpecifier() const { return fSpecifier; }

   auto Size() const { return fItems.size(); }

   void Add(std::unique_ptr<Detail::RMenuItem> &&item) { fItems.emplace_back(std::move(item)); }

   void AddMenuItem(const std::string &name, const std::string &title, const std::string &exec, const TClass *cl = nullptr)
   {
      auto item = std::make_unique<Detail::RMenuItem>(name, title);
      item->SetExec(exec);
      if (cl) item->SetClassName(cl->GetName());
      Add(std::move(item));
   }

   void AddChkMenuItem(const std::string &name, const std::string &title, bool checked, const std::string &toggle, const TClass *cl = nullptr)
   {
      auto item = std::make_unique<Detail::RCheckedMenuItem>(name, title, checked);
      item->SetExec(toggle);
      if (cl) item->SetClassName(cl->GetName());
      Add(std::move(item));
   }

   void PopulateObjectMenu(void *obj, TClass *cl);
};


/** \class RDrawableMenuRequest
\ingroup GpadROOT7
\brief Request menu items for the drawable object
\author Sergey Linev
\date 2020-04-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RDrawableMenuRequest : public RDrawableRequest {
   std::string menukind;
   std::string menureqid;
public:
   std::unique_ptr<RDrawableReply> Process() override;
};


} // namespace Experimental
} // namespace ROOT

#endif
