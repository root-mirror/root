// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveSecondarySelectable.hxx>
#include <ROOT/REveElement.hxx>

using namespace ROOT::Experimental;

/** \class REveSecondarySelectable
\ingroup REve
Semi-abstract interface for classes supporting secondary-selection.

Element class that inherits from this, should also implement the
following virtual methods from REveElement:
~~~ {.cpp}
    virtual void UnSelected();
    virtual void UnHighlighted();
~~~
and clear corresponding selection-set from there.

To support tooltips for sub-elements, implement:
~~~ {.cpp}
    virtual std::string REveElement::GetHighlightTooltip() const;
~~~
and return tooltip for the entry in the fHighlightedSet.
There should always be a single entry there.
See EveDigitSet for an example.
*/


////////////////////////////////////////////////////////////////////////////////
/// Process secondary GL selection and populate selected set accordingly.

/*

void REveSecondarySelectable::ProcessGLSelection(TGLSelectRecord& rec)
{
   if (rec.GetHighlight())
      ProcessGLSelectionInternal(rec, fHighlightedSet);
   else
      ProcessGLSelectionInternal(rec, fSelectedSet);
}

////////////////////////////////////////////////////////////////////////////////
/// Process secondary GL selection and populate given set accordingly.

void REveSecondarySelectable::ProcessGLSelectionInternal(TGLSelectRecord& rec,
                                                         SelectionSet_t& sset)
{
   Int_t id = (rec.GetN() > 1) ? (Int_t) rec.GetItem(1) : -1;

   if (sset.empty())
   {
      if (id >= 0)
      {
         sset.insert(id);
         rec.SetSecSelResult(TGLSelectRecord::kEnteringSelection);
      }
   }
   else
   {
      if (id >= 0)
      {
         if (rec.GetMultiple())
         {
            if (sset.find(id) == sset.end())
            {
               sset.insert(id);
               rec.SetSecSelResult(TGLSelectRecord::kModifyingInternalSelection);
            }
            else
            {
               sset.erase(id);
               if (sset.empty())
                  rec.SetSecSelResult(TGLSelectRecord::kLeavingSelection);
               else
                  rec.SetSecSelResult(TGLSelectRecord::kModifyingInternalSelection);
           }
         }
         else
         {
            if (sset.size() != 1 || sset.find(id) == sset.end())
            {
               sset.clear();
               sset.insert(id);
               rec.SetSecSelResult(TGLSelectRecord::kModifyingInternalSelection);
            }
         }
      }
      else
      {
         if (!rec.GetMultiple())
         {
            sset.clear();
            rec.SetSecSelResult(TGLSelectRecord::kLeavingSelection);
         }
      }
   }

   if (rec.GetSecSelResult() != TGLSelectRecord::kNone)
   {
      dynamic_cast<REveElement*>(this)->StampColorSelection();
   }
}

*/
