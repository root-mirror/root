// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveCompound
#define ROOT7_REveCompound

#include <ROOT/REveElement.hxx>
#include <ROOT/REveProjectionBases.hxx>

namespace ROOT {
namespace Experimental {

//==============================================================================
// REveCompound
// Container for managing compounds of EveElements.
//==============================================================================

class REveCompound : public REveElement,
                     public REveProjectable
{
private:
   REveCompound(const REveCompound &) = delete;
   REveCompound &operator=(const REveCompound &) = delete;

protected:
   Short_t fCompoundOpen{0}; // If more than zero, tag new children as compound members.
   Bool_t  fDoColor{kFALSE};
   Bool_t  fDoTransparency{kFALSE};

public:
   REveCompound(const std::string& n = "REveCompound", const std::string& t = "",
                Bool_t doColor = kTRUE, Bool_t doTransparency = kFALSE);
   virtual ~REveCompound() {}

   void   OpenCompound()   { ++fCompoundOpen; }
   void   CloseCompound()  { --fCompoundOpen; }
   Bool_t IsCompoundOpen() const { return fCompoundOpen > 0; }

   void SetMainColor(Color_t color) override;
   void SetMainTransparency(Char_t t) override;

   void AddElement(REveElement *el) override;
   void RemoveElementLocal(REveElement *el) override;
   void RemoveElementsLocal() override;

   void FillImpliedSelectedSet(Set_t &impSelSet) override;

   TClass *ProjectedClass(const REveProjection *p) const override;
};

//==============================================================================
// REveCompoundProjected
// Projected EveCompund container.
//==============================================================================

class REveCompoundProjected : public REveCompound,
                              public REveProjected
{
private:
   REveCompoundProjected(const REveCompoundProjected &) = delete;
   REveCompoundProjected &operator=(const REveCompoundProjected &) = delete;

public:
   REveCompoundProjected();
   virtual ~REveCompoundProjected() {}

   virtual void SetMainColor(Color_t color);

   virtual void UpdateProjection() {}
   virtual REveElement *GetProjectedAsElement() { return this; }
};

} // namespace Experimental
} // namespace ROOT

#endif
