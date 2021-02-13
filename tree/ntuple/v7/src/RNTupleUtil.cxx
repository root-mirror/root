/// \file RNTupleUtil.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch> & Max Orok <maxwellorok@gmail.com>
/// \date 2020-07-14
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RNTupleUtil.hxx"

#include "ROOT/RLogger.hxx"
#include "ROOT/RMiniFile.hxx"

#include <iostream>

ROOT::Experimental::RLogChannel &ROOT::Experimental::NTupleLog() {
   static RLogChannel sLog("ROOT.NTuple");
   return sLog;
}


namespace ROOT {
namespace Experimental {
namespace Internal {

void PrintRNTuple(const RNTuple& ntuple, std::ostream& output) {
   output << "RNTuple {\n";
   output << "    fVersion: " << ntuple.fVersion << ",\n";
   output << "    fSize: " << ntuple.fSize << ",\n";
   output << "    fSeekHeader: " << ntuple.fSeekHeader << ",\n";
   output << "    fNBytesHeader: " << ntuple.fNBytesHeader << ",\n";
   output << "    fLenHeader: " << ntuple.fLenHeader << ",\n";
   output << "    fSeekFooter: " << ntuple.fSeekFooter << ",\n";
   output << "    fNBytesFooter: " << ntuple.fNBytesFooter << ",\n";
   output << "    fLenFooter: " << ntuple.fLenFooter << ",\n";
   output << "    fReserved: " << ntuple.fReserved << ",\n";
   output << "}";
}

} // namespace Internal
} // namespace Experimental
} // namespace ROOT
