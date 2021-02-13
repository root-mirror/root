/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsabl_RSysFileItem
#define ROOT7_Browsabl_RSysFileItem

#include <ROOT/Browsable/RItem.hxx>

namespace ROOT {
namespace Experimental {
namespace Browsable {

/** \class RSysFileItem
\ingroup rbrowser
\brief Representation of single item in the file browser
*/

class RSysFileItem : public RItem {
public:
   // internal data, used for generate directory list
   int type{0};             ///<! file type
   int uid{0};              ///<! file uid
   int gid{0};              ///<! file gid
   bool islink{false};      ///<! true if symbolic link
   bool isdir{false};       ///<! true if directory
   long modtime{0};         ///<! modification time
   int64_t size{0};         ///<! file size

   // this is part for browser, visible for I/O
   std::string fsize;    ///< file size
   std::string mtime;    ///< modification time
   std::string ftype;    ///< file attributes
   std::string fuid;     ///< user id
   std::string fgid;     ///< group id

   /** Default constructor */
   RSysFileItem() = default;

   RSysFileItem(const std::string &_name, int _nchilds) : RItem(_name, _nchilds) {}

   // should be here, one needs virtual table for correct streaming of RRootBrowserReply
   virtual ~RSysFileItem() = default;

   bool IsFolder() const override { return isdir; }

   // return true for hidden files
   bool IsHidden() const override {
      auto &n = GetName();
      if ((n.length() == 0) || (n[0] != '.')) return false;
      return (n != ".") && (n != "..");
   }

   bool Compare(const RItem *b, const std::string &method) const override
   {
      if (IsFolder() != b->IsFolder())
         return IsFolder();

      if (method == "size") {
         auto fb = dynamic_cast<const RSysFileItem *> (b);
         if (fb)
            return size < fb->size;
      }

      return GetName() < b->GetName();
   }
};

} // namespace Browsable
} // namespace Experimental
} // namespace ROOT


#endif
