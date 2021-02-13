// @(#)root/tree:$Id$
// Author: Timur Pocheptsov   30/01/2014

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TreeUtils                                                            //
//                                                                      //
// Different standalone functions to work with trees and tuples,        //
// not required to be a member of any class.                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef ROOT_TreeUtils
#define ROOT_TreeUtils

#include <iosfwd>

#include "RtypesCore.h"

namespace ROOT {
namespace TreeUtils {

//1. Function to fill tuples (TNtuple/TNtupleD) from
//a simple ASCII data file. With auto and decltype - we can
//get rid of DataType parameter :) (or with a simple typedef inside ntuple class).
//An input file consists of non-empty lines (separated by newline-characters), possibly empty lines,
//and comments (treated as empty lines). Each non-empty line should contain N numbers - entry for a tuple.
//Non-strict mode lets you to have newline-characters inside a tuple's row (as it worked
//in ROOT prior to v5.3xxx).

template<class DataType, class Tuple>
Long64_t FillNtupleFromStream(std::istream &inputStream, Tuple &tuple, char delimiter, bool strictMode);

}
}

#endif
