/// \file RHistDisplayItem.cxx
/// \ingroup Hist ROOT7
/// \author Sergey Linev <s.linev@gsi.de>
/// \date 2020-06-25
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RHistDisplayItem.hxx"

using namespace ROOT::Experimental;

RHistDisplayItem::RHistDisplayItem(const RDrawable &dr) : RIndirectDisplayItem(dr)
{
}

