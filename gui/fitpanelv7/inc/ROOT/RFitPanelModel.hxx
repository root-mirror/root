// \file ROOT/RFitPanelModel.hxx
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

#ifndef ROOT_RFitPanelModel
#define ROOT_RFitPanelModel

#include "Foption.h"
#include "Fit/DataRange.h"
#include "Math/MinimizerOptions.h"

#include <vector>
#include <string>

class TH1;
class TF1;

namespace ROOT {
namespace Experimental {

/** Data structure for the fit panel */

struct RFitPanelModel {

   /// Generic item for ui5 ComboBox
   struct RComboBoxItem {
      std::string key;
      std::string value;
      RComboBoxItem() = default;
      RComboBoxItem(const std::string &_key, const std::string &_value) : key(_key), value(_value) {}
   };

   /// Entry in minimizer algorithm combo
   struct RMinimezerAlgorithm {
      int lib{0};        // to which library belongs
      int id{0};         // identifier
      std::string text;  // text shown in combobox
      RMinimezerAlgorithm() = default;
      RMinimezerAlgorithm(int _lib, int _id, const std::string &_text) : lib(_lib), id(_id), text(_text) {}
   };

   /// Basic function info, used in combo boxes
   struct RFitFuncInfo {
      std::string group;
      std::string name;
      std::string id;

      RFitFuncInfo() = default;
      RFitFuncInfo(const std::string &_name) : group("Prefefined"), name(_name) { id = "dflt::"; id.append(_name); }
   };

   /// Function parameter info, used in edit parameters dialog

   struct RFitFuncParameter {
      int ipar{0};
      std::string name;
      std::string value;
      bool fixed{false};
      std::string error;
      std::string min;
      std::string max;
      RFitFuncParameter() = default;
      RFitFuncParameter(int _ipar, const std::string &_name) : ipar(_ipar), name(_name) {}
   };

   /// Class used to transfer functions parameters list from/to client
   struct RFitFuncParsList {
      bool haspars{false};
      std::string name;
      std::vector<RFitFuncParameter> pars;
      void GetParameters(TF1 *f1);
      void SetParameters(TF1 *f1);
      void Clear();
   };

   std::string fTitle;                      ///< title of the fit panel

   std::vector<RComboBoxItem> fDataSet;     ///< list of available data sources
   std::string fSelectedData;               ///< selected data

   int fDim{0};                             ///< number of dimensions in selected data object

   std::vector<RFitFuncInfo> fFuncList;     ///< all available fit functions
   std::string fSelectedFunc;               ///< id of selected fit function like dflt::gaus


   // General tab

   // Method
   std::vector<RComboBoxItem> fFitMethods;   ///< all supported fit methods
   std::string fFitMethod;                   ///< selected fit method

   bool fLinearFit{false};
   bool fRobust{false};
   float fRobustLevel{0.95};

   // Fit Options
   bool fIntegral{false};
   bool fUseRange{false};
   bool fBestErrors{false};
   bool fImproveFitResults{false};
   bool fAllWeights1{false};
   bool fAddToList{false};
   bool fEmptyBins1{false};
   bool fUseGradient{false};

   // Draw Options
   bool fSame{false};
   bool fNoDrawing{false};
   bool fNoStoreDraw{false};


   // all combo items for all methods

   // Minimization Tab
   int fLibrary{0};                                 ///< selected minimization library
   bool fHasGenetics{false};                        ///< is genetics available
   std::vector<RMinimezerAlgorithm> fMethodMinAll;  ///< all items for all methods
   int fSelectMethodMin{0};

   float fErrorDef{1.00};
   float fMaxTolerance{0.01};
   int fMaxIterations{0};
   int fPrint{0};


   // range selection, shown depending on fDim
   float fMinRangeX{0};
   float fMaxRangeX{1};
   float fStepX{0.01};
   float fRangeX[2] = {0,1};

   float fMinRangeY{0};
   float fMaxRangeY{1};
   float fStepY{0.01};
   float fRangeY[2] = {0,1};


   /// Parameters

   RFitFuncParsList fFuncPars;

   /////////Advanced Options

   bool fHasAdvanced{false};
   std::vector<RComboBoxItem> fContour1;
   std::string fContourPar1Id;
   std::vector<RComboBoxItem> fContour2;
   std::string fContourPar2Id;
   std::vector<RComboBoxItem> fScan;
   std::string fScanId;

   // Contour Tab
   int fContourPar1{0};
   int fContourPar2{0};
   int fContourPoints{0};
   float fConfLevel{0.};
   bool fContourImpose{false};
   std::string fColorContour[3];

   // Scan Tab
   int fScanPoints{0};
   int fScanPar{0};
   int fScanMin{0};
   int fScanMax{0};

   RFitPanelModel() { Initialize(); }

   void Initialize();

   bool SelectHistogram(const std::string &hname, TH1 *hist);

   void SelectedFunc(const std::string &name, TF1 *func);

   void UpdateRange(TH1 *hist);

   void UpdateAdvanced(TF1 *func);

   void UpdateFuncList();

   bool IsSelectedHistogram() const { return !fSelectedData.empty(); }

   TH1* GetSelectedHistogram(TH1 *hist = nullptr);

   void GetRanges(ROOT::Fit::DataRange &drange);
   void GetFitOptions(Foption_t &fitOpts);
   void GetMinimizerOptions(ROOT::Math::MinimizerOptions &minOpts);

   std::string GetDrawOption();

   std::string GetFitOption();
};

} // namespace Experimental
} // namespace ROOT

#endif
