
#ifndef HIST_FACTORY_MODEL_UTILS_H
#define HIST_FACTORY_MODEL_UTILS_H

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooStats/HistFactory/ParamHistFunc.h"

#include <vector>
#include <map>
#include <string>

namespace RooStats {
namespace HistFactory {
  ///\ingroup HistFactory
  std::string channelNameFromPdf( RooAbsPdf* channelPdf );

//   void getChannelsFromModel( RooAbsPdf* model, RooArgSet* channels, 
// 			     RooArgSet* channelsWithConstraints );

  ///\ingroup HistFactory
  void FactorizeHistFactoryPdf(const RooArgSet&, RooAbsPdf&, RooArgList&, RooArgList&);
  ///\ingroup HistFactory
  bool getStatUncertaintyFromChannel( RooAbsPdf* channel, ParamHistFunc*& paramfunc, 
				      RooArgList* gammaList );

  ///\ingroup HistFactory
  RooAbsPdf* getSumPdfFromChannel( RooAbsPdf* channel );

  ///\ingroup HistFactory
  void getDataValuesForObservables( std::map< std::string, std::vector<double> >& ChannelBinDataMap, 
				    RooAbsData* data, RooAbsPdf* simPdf );

  ///\ingroup HistFactory
  int getStatUncertaintyConstraintTerm( RooArgList* constraints, RooRealVar* gamma_stat, 
					RooAbsReal*& pois_mean, RooRealVar*& tau );

}
}



#endif
