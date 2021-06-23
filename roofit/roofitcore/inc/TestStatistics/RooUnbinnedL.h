/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooUnbinnedL.h,v 1.10 2007/07/21 21:32:52 wouter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *                                                                           *
 * Copyright (c) 2000-2020, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOT_ROOFIT_TESTSTATISTICS_RooUnbinnedL
#define ROOT_ROOFIT_TESTSTATISTICS_RooUnbinnedL

#include <TestStatistics/RooAbsL.h>

// forward declarations
class RooAbsPdf;
class RooAbsData;
class RooArgSet;

namespace RooFit {
namespace TestStatistics {

class RooUnbinnedL :
   public RooAbsL {
public:
   RooUnbinnedL(RooAbsPdf* pdf, RooAbsData* data, RooAbsL::Extended extended = RooAbsL::Extended::Auto);
   RooUnbinnedL(const RooUnbinnedL &other);
   bool set_apply_weight_squared(bool flag);

   double evaluate_partition(Section events, std::size_t components_begin,
                             std::size_t components_end) override;

private:
   bool processEmptyDataSets() const;
   bool apply_weight_squared = false;            // Apply weights squared?
   mutable bool _first = true;       //!
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_RooUnbinnedL
