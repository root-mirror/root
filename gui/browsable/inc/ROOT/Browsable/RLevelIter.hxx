/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_RLevelIter
#define ROOT7_Browsable_RLevelIter

#include <memory>
#include <string>

namespace ROOT {
namespace Experimental {
namespace Browsable {

class RElement;
class RItem;

/** \class RLevelIter
\ingroup rbrowser
\brief Iterator over single level hierarchy like any array, keys list, ...
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RLevelIter {
public:
   virtual ~RLevelIter() = default;

   /** Shift to next entry */
   virtual bool Next() = 0;

   /** Returns current entry name  */
   virtual std::string GetItemName() const = 0;

   /** Returns true if current item can have childs */
   virtual bool CanItemHaveChilds() const { return false; }

   /** Create RElement for current entry - may take much time to load object or open file */
   virtual std::shared_ptr<RElement> GetElement() = 0;

   virtual std::unique_ptr<RItem> CreateItem();

   virtual bool Find(const std::string &name, int indx = -1);

};


} // namespace Browsable
} // namespace Experimental
} // namespace ROOT

#endif
