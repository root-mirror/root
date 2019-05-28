// \file ROOT/RFitPanel6.hxx
/// \ingroup WebGui ROOT7
/// \author Sergey Linev <S.Linev@gsi.de>
/// \author Iliana Betsou <Iliana.Betsou@cern.ch>
/// \date 2019-04-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RFitPanel6
#define ROOT_RFitPanel6

#include <ROOT/RWebWindow.hxx>

#include <ROOT/RFitPanel6Model.hxx>

#include "TH1.h"

#include <vector>

namespace ROOT {
namespace Experimental {

class RFitPanel6 {

   std::string fTitle;  ///<! title
   unsigned fConnId{0}; ///<! connection id
   TH1 *fHist{nullptr};

   std::shared_ptr<RWebWindow> fWindow; ///!< configured display

   /// process data from UI
   void ProcessData(unsigned connid, const std::string &arg);

   void DoFit(const std::string &model);

   void DrawContour(const std::string &model);

   void DrawScan(const std::string &model);

public:
   /// normal constructor
   RFitPanel6(const std::string &title = "Fit panel") : fTitle(title) {}

   /// destructor
   virtual ~RFitPanel6() {}

   // method required when any panel want to be inserted into the RCanvas
   std::shared_ptr<RWebWindow> GetWindow();

   void AssignHistogram(TH1 *hist) { fHist = hist; }

   /// show FitPanel in specified place
   void Show(const std::string &where = "");

   /// hide FitPanel
   void Hide();
};

} // namespace Experimental
} // namespace ROOT

#endif
