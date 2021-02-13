/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_TObjectElement
#define ROOT7_Browsable_TObjectElement

#include <ROOT/Browsable/RElement.hxx>

class TObject;
class TCollection;

namespace ROOT {
namespace Experimental {
namespace Browsable {


/** \class TObjectElement
\ingroup rbrowser
\brief Access to TObject basic properties for RBrowsable
\author Sergey Linev <S.Linev@gsi.de>
\date 2021-01-11
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


class TObjectElement : public RElement {
protected:
   std::unique_ptr<RHolder> fObject;
   TObject *fObj{nullptr};
   std::string fName;

   bool IsSame(TObject *obj) const { return obj == fObj; }

public:
   TObjectElement(TObject *obj, const std::string &name = "");

   TObjectElement(std::unique_ptr<RHolder> &obj, const std::string &name = "");

   virtual ~TObjectElement() = default;

   /** Name of TObject */
   std::string GetName() const override;

   void SetName(const std::string &name) { fName = name; }

   /** Title of TObject */
   std::string GetTitle() const override;

   bool IsFolder();

   /** Create iterator for childs elements if any */
   std::unique_ptr<RLevelIter> GetChildsIter() override;

   /** Return copy of TObject holder - if possible */
   std::unique_ptr<RHolder> GetObject() override;

   const TClass *GetClass() const;

   EActionKind GetDefaultAction() const override;

   bool IsCapable(EActionKind) const override;

   static std::unique_ptr<RLevelIter> GetCollectionIter(const TCollection *);

};

} // namespace Browsable
} // namespace Experimental
} // namespace ROOT


#endif
