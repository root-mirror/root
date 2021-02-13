/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooSimultaneous.cxx
\class RooSimultaneous
\ingroup Roofitcore

RooSimultaneous facilitates simultaneous fitting of multiple PDFs
to subsets of a given dataset.
The class takes an index category, which is used as a selector
for PDFs, and a list of PDFs, each associated
with a state of the index category. RooSimultaneous always returns
the value of the PDF that is associated with the current value
of the index category.

Extended likelihood fitting is supported if all components support
extended likelihood mode. The expected number of events by a RooSimultaneous
is that of the component p.d.f. selected by the index category.

The index category can be accessed using indexCategory().

###Generating events
When generating events from a RooSimultaneous, the index category has to be added to
the dataset. Further, the PDF needs to know the relative probabilities of each category, i.e.,
how many events are in which category. This can be achieved in two ways:
- Generating with proto data that have category entries: An event from the same category as
in the proto data is created for each event in the proto data.
See RooAbsPdf::generate(const RooArgSet&,const RooDataSet&,Int_t,Bool_t,Bool_t,Bool_t) const.
- No proto data: A category is chosen randomly.
\note This requires that the PDFs building the simultaneous are extended. In this way,
the relative probability of each category can be calculated from the number of events
in each category.
**/

#include "RooFit.h"

#include "RooSimultaneous.h"
#include "RooAbsCategoryLValue.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
#include "Roo1DTable.h"
#include "RooSimGenContext.h"
#include "RooSimSplitGenContext.h"
#include "RooDataSet.h"
#include "RooCmdConfig.h"
#include "RooNameReg.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooDataHist.h"
#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooHelpers.h"

#include <iostream>

using namespace std;

ClassImp(RooSimultaneous);



////////////////////////////////////////////////////////////////////////////////
/// Constructor with index category. PDFs associated with indexCat
/// states can be added after construction with the addPdf() function.
///
/// RooSimultaneous can function without having a PDF associated
/// with every single state. The normalization in such cases is taken
/// from the number of registered PDFs, but getVal() will assert if
/// when called for an unregistered index state.

RooSimultaneous::RooSimultaneous(const char *name, const char *title,
				 RooAbsCategoryLValue& inIndexCat) :
  RooAbsPdf(name,title),
  _plotCoefNormSet("!plotCoefNormSet","plotCoefNormSet",this,kFALSE,kFALSE),
  _plotCoefNormRange(0),
  _partIntMgr(this,10),
  _indexCat("indexCat","Index category",this,inIndexCat),
  _numPdf(0)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from index category and full list of PDFs.
/// In this constructor form, a PDF must be supplied for each indexCat state
/// to avoid ambiguities. The PDFs are associated with the states of the
/// index category as they appear when iterating through the category states
/// with RooAbsCategory::begin() and RooAbsCategory::end(). This usually means
/// they are associated by ascending index numbers.
///
/// PDFs may not overlap (i.e. share any variables) with the index category (function)

RooSimultaneous::RooSimultaneous(const char *name, const char *title,
				 const RooArgList& inPdfList, RooAbsCategoryLValue& inIndexCat) :
  RooAbsPdf(name,title),
  _plotCoefNormSet("!plotCoefNormSet","plotCoefNormSet",this,kFALSE,kFALSE),
  _plotCoefNormRange(0),
  _partIntMgr(this,10),
  _indexCat("indexCat","Index category",this,inIndexCat),
  _numPdf(0)
{
  if (inPdfList.size() != inIndexCat.size()) {
    coutE(InputArguments) << "RooSimultaneous::ctor(" << GetName()
			  << " ERROR: Number PDF list entries must match number of index category states, no PDFs added" << endl ;
    return ;
  }

  map<string,RooAbsPdf*> pdfMap ;
  auto indexCatIt = inIndexCat.begin();
  for (unsigned int i=0; i < inPdfList.size(); ++i) {
    auto pdf = static_cast<RooAbsPdf*>(&inPdfList[i]);
    const auto& nameIdx = (*indexCatIt++);
    pdfMap[nameIdx.first] = pdf;
  }

  initialize(inIndexCat,pdfMap) ;
}


////////////////////////////////////////////////////////////////////////////////

RooSimultaneous::RooSimultaneous(const char *name, const char *title,
				 map<string,RooAbsPdf*> pdfMap, RooAbsCategoryLValue& inIndexCat) :
  RooAbsPdf(name,title),
  _plotCoefNormSet("!plotCoefNormSet","plotCoefNormSet",this,kFALSE,kFALSE),
  _plotCoefNormRange(0),
  _partIntMgr(this,10),
  _indexCat("indexCat","Index category",this,inIndexCat),
  _numPdf(0)
{
  initialize(inIndexCat,pdfMap) ;
}




// This class cannot be locally defined in initialize as it cannot be
// used as a template argument in that case
namespace RooSimultaneousAux {
  struct CompInfo {
    RooAbsPdf* pdf ;
    RooSimultaneous* simPdf ;
    const RooAbsCategoryLValue* subIndex ;
    RooArgSet* subIndexComps ;
  } ;
}

void RooSimultaneous::initialize(RooAbsCategoryLValue& inIndexCat, std::map<std::string,RooAbsPdf*> pdfMap)
{
  // First see if there are any RooSimultaneous input components
  Bool_t simComps(kFALSE) ;
  for (map<string,RooAbsPdf*>::iterator iter=pdfMap.begin() ; iter!=pdfMap.end() ; ++iter) {
    if (dynamic_cast<RooSimultaneous*>(iter->second)) {
      simComps = kTRUE ;
      break ;
    }
  }

  // If there are no simultaneous component p.d.f. do simple processing through addPdf()
  if (!simComps) {
    bool failure = false;
    for (map<string,RooAbsPdf*>::iterator iter=pdfMap.begin() ; iter!=pdfMap.end() ; ++iter) {
      failure |= addPdf(*iter->second,iter->first.c_str()) ;
    }

    if (failure) {
      throw std::invalid_argument(std::string("At least one of the PDFs of the RooSimultaneous ")
      + GetName() + " is invalid.");
    }
    return ;
  }

  // Issue info message that we are about to do some rearraning
  coutI(InputArguments) << "RooSimultaneous::initialize(" << GetName() << ") INFO: one or more input component of simultaneous p.d.f.s are"
			<< " simultaneous p.d.f.s themselves, rewriting composite expressions as one-level simultaneous p.d.f. in terms of"
			<< " final constituents and extended index category" << endl ;


  RooArgSet allAuxCats ;
  map<string,RooSimultaneousAux::CompInfo> compMap ;
  for (map<string,RooAbsPdf*>::iterator iter=pdfMap.begin() ; iter!=pdfMap.end() ; ++iter) {
    RooSimultaneousAux::CompInfo ci ;
    ci.pdf = iter->second ;
    RooSimultaneous* simComp = dynamic_cast<RooSimultaneous*>(iter->second) ;
    if (simComp) {
      ci.simPdf = simComp ;
      ci.subIndex = &simComp->indexCat() ;
      ci.subIndexComps = simComp->indexCat().isFundamental() ? new RooArgSet(simComp->indexCat()) : simComp->indexCat().getVariables() ;
      allAuxCats.add(*(ci.subIndexComps),kTRUE) ;
    } else {
      ci.simPdf = 0 ;
      ci.subIndex = 0 ;
      ci.subIndexComps = 0 ;
    }
    compMap[iter->first] = ci ;
  }

  // Construct the 'superIndex' from the nominal index category and all auxiliary components
  RooArgSet allCats(inIndexCat) ;
  allCats.add(allAuxCats) ;
  string siname = Form("%s_index",GetName()) ;
  RooSuperCategory* superIndex = new RooSuperCategory(siname.c_str(),siname.c_str(),allCats) ;
  bool failure = false;

  // Now process each of original pdf/state map entries
  for (map<string,RooSimultaneousAux::CompInfo>::iterator citer = compMap.begin() ; citer != compMap.end() ; ++citer) {

    RooArgSet repliCats(allAuxCats) ;
    if (citer->second.subIndexComps) {
      repliCats.remove(*citer->second.subIndexComps) ;
      delete citer->second.subIndexComps ;
    }
    inIndexCat.setLabel(citer->first.c_str()) ;

    if (!citer->second.simPdf) {

      // Entry is a plain p.d.f. assign it to every state permutation of the repliCats set
      RooSuperCategory repliSuperCat("tmp","tmp",repliCats) ;

      // Iterator over all states of repliSuperCat
      for (const auto& nameIdx : repliSuperCat) {
        // Set value
        repliSuperCat.setLabel(nameIdx.first) ;
        // Retrieve corresponding label of superIndex
        string superLabel = superIndex->getCurrentLabel() ;
        failure |= addPdf(*citer->second.pdf,superLabel.c_str()) ;
        cxcoutD(InputArguments) << "RooSimultaneous::initialize(" << GetName()
				    << ") assigning pdf " << citer->second.pdf->GetName() << " to super label " << superLabel << endl ;
      }
    } else {

      // Entry is a simultaneous p.d.f

      if (repliCats.getSize()==0) {

        // Case 1 -- No replication of components of RooSim component are required

        for (const auto& type : *citer->second.subIndex) {
          const_cast<RooAbsCategoryLValue*>(citer->second.subIndex)->setLabel(type.first.c_str());
          string superLabel = superIndex->getCurrentLabel() ;
          RooAbsPdf* compPdf = citer->second.simPdf->getPdf(type.first.c_str());
          if (compPdf) {
            failure |= addPdf(*compPdf,superLabel.c_str()) ;
            cxcoutD(InputArguments) << "RooSimultaneous::initialize(" << GetName()
				        << ") assigning pdf " << compPdf->GetName() << "(member of " << citer->second.pdf->GetName()
				        << ") to super label " << superLabel << endl ;
          } else {
            coutW(InputArguments) << "RooSimultaneous::initialize(" << GetName() << ") WARNING: No p.d.f. associated with label "
                << type.second << " for component RooSimultaneous p.d.f " << citer->second.pdf->GetName()
                << "which is associated with master index label " << citer->first << endl ;
          }
        }

      } else {

        // Case 2 -- Replication of components of RooSim component are required

        // Make replication supercat
        RooSuperCategory repliSuperCat("tmp","tmp",repliCats) ;

        for (const auto& stype : *citer->second.subIndex) {
          const_cast<RooAbsCategoryLValue*>(citer->second.subIndex)->setLabel(stype.first.c_str());

          for (const auto& nameIdx : repliSuperCat) {
            repliSuperCat.setLabel(nameIdx.first) ;
            const string superLabel = superIndex->getCurrentLabel() ;
            RooAbsPdf* compPdf = citer->second.simPdf->getPdf(stype.first.c_str());
            if (compPdf) {
              failure |= addPdf(*compPdf,superLabel.c_str()) ;
              cxcoutD(InputArguments) << "RooSimultaneous::initialize(" << GetName()
				          << ") assigning pdf " << compPdf->GetName() << "(member of " << citer->second.pdf->GetName()
				          << ") to super label " << superLabel << endl ;
            } else {
              coutW(InputArguments) << "RooSimultaneous::initialize(" << GetName() << ") WARNING: No p.d.f. associated with label "
                  << stype.second << " for component RooSimultaneous p.d.f " << citer->second.pdf->GetName()
                  << "which is associated with master index label " << citer->first << endl ;
            }
          }
        }
      }
    }
  }

  if (failure) {
    throw std::invalid_argument(std::string("Failed to initialise RooSimultaneous ") + GetName());
  }

  // Change original master index to super index and take ownership of it
  _indexCat.setArg(*superIndex) ;
  addOwnedComponents(*superIndex) ;

}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooSimultaneous::RooSimultaneous(const RooSimultaneous& other, const char* name) :
  RooAbsPdf(other,name),
  _plotCoefNormSet("!plotCoefNormSet",this,other._plotCoefNormSet),
  _plotCoefNormRange(other._plotCoefNormRange),
  _partIntMgr(other._partIntMgr,this),
  _indexCat("indexCat",this,other._indexCat),
  _numPdf(other._numPdf)
{
  // Copy proxy list
  TIterator* pIter = other._pdfProxyList.MakeIterator() ;
  RooRealProxy* proxy ;
  while ((proxy=(RooRealProxy*)pIter->Next())) {
    _pdfProxyList.Add(new RooRealProxy(proxy->GetName(),this,*proxy)) ;
  }
  delete pIter ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSimultaneous::~RooSimultaneous()
{
  _pdfProxyList.Delete() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the p.d.f associated with the given index category name

RooAbsPdf* RooSimultaneous::getPdf(const char* catName) const
{
  RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.FindObject(catName) ;
  return proxy ? ((RooAbsPdf*)proxy->absArg()) : 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Associate given PDF with index category state label 'catLabel'.
/// The name state must be already defined in the index category.
///
/// RooSimultaneous can function without having a PDF associated
/// with every single state. The normalization in such cases is taken
/// from the number of registered PDFs, but getVal() will fail if
/// called for an unregistered index state.
///
/// PDFs may not overlap (i.e. share any variables) with the index category (function).
/// \param[in] pdf PDF to be added.
/// \param[in] catLabel Name of the category state to be associated to the PDF.
/// \return `true` in case of failure.

Bool_t RooSimultaneous::addPdf(const RooAbsPdf& pdf, const char* catLabel)
{
  // PDFs cannot overlap with the index category
  if (pdf.dependsOn(_indexCat.arg())) {
    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName() << "): PDF '" << pdf.GetName()
			  << "' overlaps with index category '" << _indexCat.arg().GetName() << "'."<< endl ;
    return kTRUE ;
  }

  // Each index state can only have one PDF associated with it
  if (_pdfProxyList.FindObject(catLabel)) {
    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName() << "): index state '"
			  << catLabel << "' has already an associated PDF." << endl ;
    return kTRUE ;
  }

  const RooSimultaneous* simPdf = dynamic_cast<const RooSimultaneous*>(&pdf) ;
  if (simPdf) {

    coutE(InputArguments) << "RooSimultaneous::addPdf(" << GetName()
			  << ") ERROR: you cannot add a RooSimultaneous component to a RooSimultaneous using addPdf()."
			  << " Use the constructor with RooArgList if input p.d.f.s or the map<string,RooAbsPdf&> instead." << endl ;
    return kTRUE ;

  } else {

    // Create a proxy named after the associated index state
    TObject* proxy = new RooRealProxy(catLabel,catLabel,this,(RooAbsPdf&)pdf) ;
    _pdfProxyList.Add(proxy) ;
    _numPdf += 1 ;
  }

  return kFALSE ;
}





////////////////////////////////////////////////////////////////////////////////
/// Examine the pdf components and check if one of them can be extended or must be extended
/// It is enough to have one component that can be exteded or must be extended to return the flag in
///  the total simultaneous pdf

RooAbsPdf::ExtendMode RooSimultaneous::extendMode() const
{
  Bool_t anyCanExtend(kFALSE) ;
  Bool_t anyMustExtend(kFALSE) ;

  for (Int_t i=0 ; i<_numPdf ; i++) {
    RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.At(i);
    if (proxy) {
      RooAbsPdf* pdf = (RooAbsPdf*) proxy->absArg() ;
      //cout << " now processing pdf " << pdf->GetName() << endl;
      if (pdf->canBeExtended()) {
         //cout << "RooSim::extendedMode(" << GetName() << ") component " << pdf->GetName() << " can be extended"
         //     << endl;
         anyCanExtend = kTRUE;
      }
      if (pdf->mustBeExtended()) {
         //cout << "RooSim::extendedMode(" << GetName() << ") component " << pdf->GetName() << " MUST be extended" << endl;
         anyMustExtend = kTRUE;
      }
    }
  }
  if (anyMustExtend) {
    //cout << "RooSim::extendedMode(" << GetName() << ") returning MustBeExtended" << endl ;
    return MustBeExtended ;
  }
  if (anyCanExtend) {
    //cout << "RooSim::extendedMode(" << GetName() << ") returning CanBeExtended" << endl ;
    return CanBeExtended ;
  }
  //cout << "RooSim::extendedMode(" << GetName() << ") returning CanNotBeExtended" << endl ;
  return CanNotBeExtended ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return the current value:
/// the value of the PDF associated with the current index category state

Double_t RooSimultaneous::evaluate() const
{
  // Retrieve the proxy by index name
  RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.FindObject(_indexCat.label()) ;

  //assert(proxy!=0) ;
  if (proxy==0) return 0 ;

  // Calculate relative weighting factor for sim-pdfs of all extendable components
  Double_t catFrac(1) ;
  if (canBeExtended()) {
    Double_t nEvtCat = ((RooAbsPdf*)(proxy->absArg()))->expectedEvents(_normSet) ;

    Double_t nEvtTot(0) ;
    TIterator* iter = _pdfProxyList.MakeIterator() ;
    RooRealProxy* proxy2 ;
    while((proxy2=(RooRealProxy*)iter->Next())) {
      nEvtTot += ((RooAbsPdf*)(proxy2->absArg()))->expectedEvents(_normSet) ;
    }
    delete iter ;
    catFrac=nEvtCat/nEvtTot ;
  }

  // Return the selected PDF value, normalized by the number of index states
  return ((RooAbsPdf*)(proxy->absArg()))->getVal(_normSet)*catFrac ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the number of expected events: If the index is in nset,
/// then return the sum of the expected events of all components,
/// otherwise return the number of expected events of the PDF
/// associated with the current index category state

Double_t RooSimultaneous::expectedEvents(const RooArgSet* nset) const
{
  if (nset->contains(_indexCat.arg())) {

    Double_t sum(0) ;

    TIterator* iter = _pdfProxyList.MakeIterator() ;
    RooRealProxy* proxy ;
    while((proxy=(RooRealProxy*)iter->Next())) {
      sum += ((RooAbsPdf*)(proxy->absArg()))->expectedEvents(nset) ;
    }
    delete iter ;

    return sum ;

  } else {

    // Retrieve the proxy by index name
    RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.FindObject(_indexCat.label()) ;

    //assert(proxy!=0) ;
    if (proxy==0) return 0 ;

    // Return the selected PDF value, normalized by the number of index states
    return ((RooAbsPdf*)(proxy->absArg()))->expectedEvents(nset);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward determination of analytical integration capabilities to component p.d.f.s
/// A unique code is assigned to the combined integration capabilities of all associated
/// p.d.f.s

Int_t RooSimultaneous::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars,
					       const RooArgSet* normSet, const char* rangeName) const
{
  // Declare that we can analytically integrate all requested observables
  analVars.add(allVars) ;

  // Retrieve (or create) the required partial integral list
  Int_t code ;

  // Check if this configuration was created before
  CacheElem* cache = (CacheElem*) _partIntMgr.getObj(normSet,&analVars,0,RooNameReg::ptr(rangeName)) ;
  if (cache) {
    code = _partIntMgr.lastIndex() ;
    return code+1 ;
  }
  cache = new CacheElem ;

  // Create the partial integral set for this request
  TIterator* iter = _pdfProxyList.MakeIterator() ;
  RooRealProxy* proxy ;
  while((proxy=(RooRealProxy*)iter->Next())) {
    RooAbsReal* pdfInt = proxy->arg().createIntegral(analVars,normSet,0,rangeName) ;
    cache->_partIntList.addOwned(*pdfInt) ;
  }
  delete iter ;

  // Store the partial integral list and return the assigned code ;
  code = _partIntMgr.setObj(normSet,&analVars,cache,RooNameReg::ptr(rangeName)) ;

  return code+1 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return analytical integration defined by given code

Double_t RooSimultaneous::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* /*rangeName*/) const
{
  // No integration scenario
  if (code==0) {
    return getVal(normSet) ;
  }

  // Partial integration scenarios, rangeName already encoded in 'code'
  CacheElem* cache = (CacheElem*) _partIntMgr.getObjByIndex(code-1) ;

  RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.FindObject(_indexCat.label()) ;
  Int_t idx = _pdfProxyList.IndexOf(proxy) ;
  return ((RooAbsReal*)cache->_partIntList.at(idx))->getVal(normSet) ;
}






////////////////////////////////////////////////////////////////////////////////
/// Back-end for plotOn() implementation on RooSimultaneous which
/// needs special handling because a RooSimultaneous PDF cannot
/// project out its index category via integration. plotOn() will
/// abort if this is requested without providing a projection dataset.

RooPlot* RooSimultaneous::plotOn(RooPlot *frame, RooLinkedList& cmdList) const
{
  // Sanity checks
  if (plotSanityChecks(frame)) return frame ;

  // Extract projection configuration from command list
  RooCmdConfig pc(Form("RooSimultaneous::plotOn(%s)",GetName())) ;
  pc.defineString("sliceCatState","SliceCat",0,"",kTRUE) ;
  pc.defineDouble("scaleFactor","Normalization",0,1.0) ;
  pc.defineInt("scaleType","Normalization",0,RooAbsPdf::Relative) ;
  pc.defineObject("sliceCatList","SliceCat",0,0,kTRUE) ;
  pc.defineObject("projSet","Project",0) ;
  pc.defineObject("sliceSet","SliceVars",0) ;
  pc.defineObject("projDataSet","ProjData",0) ;
  pc.defineObject("projData","ProjData",1) ;
  pc.defineMutex("Project","SliceVars") ;
  pc.allowUndefined() ; // there may be commands we don't handle here

  // Process and check varargs
  pc.process(cmdList) ;
  if (!pc.ok(kTRUE)) {
    return frame ;
  }

  const RooAbsData* projData = (const RooAbsData*) pc.getObject("projData") ;
  const RooArgSet* projDataSet = (const RooArgSet*) pc.getObject("projDataSet") ;
  const RooArgSet* sliceSetTmp = (const RooArgSet*) pc.getObject("sliceSet") ;
  std::unique_ptr<RooArgSet> sliceSet( sliceSetTmp ? ((RooArgSet*) sliceSetTmp->Clone()) : nullptr );
  const RooArgSet* projSet = (const RooArgSet*) pc.getObject("projSet") ;
  Double_t scaleFactor = pc.getDouble("scaleFactor") ;
  ScaleType stype = (ScaleType) pc.getInt("scaleType") ;


  // Look for category slice arguments and add them to the master slice list if found
  const char* sliceCatState = pc.getString("sliceCatState",0,kTRUE) ;
  const RooLinkedList& sliceCatList = pc.getObjectList("sliceCatList") ;
  if (sliceCatState) {

    // Make the master slice set if it doesnt exist
    if (!sliceSet) {
      sliceSet.reset(new RooArgSet);
    }

    // Prepare comma separated label list for parsing
    auto catTokens = RooHelpers::tokenise(sliceCatState, ",");

    // Loop over all categories provided by (multiple) Slice() arguments
    TIterator* iter = sliceCatList.MakeIterator() ;
    RooCategory* scat ;
    unsigned int tokenIndex = 0;
    while((scat=(RooCategory*)iter->Next())) {
      const char* slabel = tokenIndex >= catTokens.size() ? nullptr : catTokens[tokenIndex++].c_str();

      if (slabel) {
        // Set the slice position to the value indicated by slabel
        scat->setLabel(slabel) ;
        // Add the slice category to the master slice set
        sliceSet->add(*scat,kFALSE) ;
      }
    }
    delete iter ;
  }

  // Check if we have a projection dataset
  if (!projData) {
    coutE(InputArguments) << "RooSimultaneous::plotOn(" << GetName() << ") ERROR: must have a projection dataset for index category" << endl ;
    return frame ;
  }

  // Make list of variables to be projected
  RooArgSet projectedVars ;
  if (sliceSet) {
    makeProjectionSet(frame->getPlotVar(),frame->getNormVars(),projectedVars,kTRUE) ;

    // Take out the sliced variables
    for (const auto sliceArg : *sliceSet) {
      RooAbsArg* arg = projectedVars.find(sliceArg->GetName()) ;
      if (arg) {
        projectedVars.remove(*arg) ;
      } else {
        coutI(Plotting) << "RooAbsReal::plotOn(" << GetName() << ") slice variable "
            << sliceArg->GetName() << " was not projected anyway" << endl ;
      }
    }
  } else if (projSet) {
    makeProjectionSet(frame->getPlotVar(),projSet,projectedVars,kFALSE) ;
  } else {
    makeProjectionSet(frame->getPlotVar(),frame->getNormVars(),projectedVars,kTRUE) ;
  }

  Bool_t projIndex(kFALSE) ;

  if (!_indexCat.arg().isDerived()) {
    // *** Error checking for a fundamental index category ***
    //cout << "RooSim::plotOn: index is fundamental" << endl ;

    // Check that the provided projection dataset contains our index variable
    if (!projData->get()->find(_indexCat.arg().GetName())) {
      coutE(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") ERROR: Projection over index category "
		      << "requested, but projection data set doesn't contain index category" << endl ;
      return frame ;
    }

    if (projectedVars.find(_indexCat.arg().GetName())) {
      projIndex=kTRUE ;
    }

  } else {
    // *** Error checking for a composite index category ***

    // Determine if any servers of the index category are in the projectedVars
    RooArgSet projIdxServers ;
    Bool_t anyServers(kFALSE) ;
    for (const auto server : _indexCat->servers()) {
      if (projectedVars.find(server->GetName())) {
        anyServers=kTRUE ;
        projIdxServers.add(*server) ;
      }
    }

    // Check that the projection dataset contains all the
    // index category components we're projecting over

    // Determine if all projected servers of the index category are in the projection dataset
    Bool_t allServers(kTRUE) ;
    std::string missing;
    for (const auto server : projIdxServers) {
      if (!projData->get()->find(server->GetName())) {
        allServers=kFALSE ;
        missing = server->GetName();
      }
    }

    if (!allServers) {
      coutE(Plotting) << "RooSimultaneous::plotOn(" << GetName()
	       << ") ERROR: Projection dataset doesn't contain complete set of index categories to do projection."
	       << "\n\tcategory " << missing << " is missing." << endl ;
      return frame ;
    }

    if (anyServers) {
      projIndex = kTRUE ;
    }
  }

  // Calculate relative weight fractions of components
  std::unique_ptr<Roo1DTable> wTable( projData->table(_indexCat.arg()) );

  // Clone the index category to be able to cycle through the category states for plotting without
  // affecting the category state of our instance
  std::unique_ptr<RooArgSet> idxCloneSet( RooArgSet(*_indexCat).snapshot(true) );
  auto idxCatClone = static_cast<RooAbsCategoryLValue*>( idxCloneSet->find(_indexCat->GetName()) );
  assert(idxCatClone);

  // Make list of category columns to exclude from projection data
  std::unique_ptr<RooArgSet> idxCompSliceSet( idxCatClone->getObservables(frame->getNormVars()) );

  // If we don't project over the index, just do the regular plotOn
  if (!projIndex) {

    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
		    << " represents a slice in the index category ("  << _indexCat.arg().GetName() << ")" << endl ;

    // Reduce projData: take out fitCat (component) columns and entries that don't match selected slice
    // Construct cut string to only select projection data event that match the current slice

    // Make cut string to exclude rows from projection data
    TString cutString ;
    Bool_t first(kTRUE) ;
    for (const auto arg : *idxCompSliceSet) {
      auto idxComp = static_cast<RooCategory*>(arg);
      RooAbsArg* slicedComponent = nullptr;
      if (sliceSet && (slicedComponent = sliceSet->find(*idxComp)) != nullptr) {
        auto theCat = static_cast<const RooAbsCategory*>(slicedComponent);
        idxComp->setIndex(theCat->getCurrentIndex(), false);
      }

      if (!first) {
        cutString.Append("&&") ;
      } else {
        first=kFALSE ;
      }
      cutString.Append(Form("%s==%d",idxComp->GetName(),idxComp->getCurrentIndex())) ;
    }

    // Make temporary projData without RooSim index category components
    RooArgSet projDataVars(*projData->get()) ;
    projDataVars.remove(*idxCompSliceSet,kTRUE,kTRUE) ;

    std::unique_ptr<RooAbsData> projDataTmp( const_cast<RooAbsData*>(projData)->reduce(projDataVars,cutString) );

    // Override normalization and projection dataset
    RooCmdArg tmp1 = RooFit::Normalization(scaleFactor*wTable->getFrac(idxCatClone->getCurrentLabel()),stype) ;
    RooCmdArg tmp2 = RooFit::ProjWData(*projDataSet,*projDataTmp) ;

    // WVE -- do not adjust normalization for asymmetry plots
    RooLinkedList cmdList2(cmdList) ;
    if (!cmdList.find("Asymmetry")) {
      cmdList2.Add(&tmp1) ;
    }
    cmdList2.Add(&tmp2) ;

    // Plot single component
    RooPlot* retFrame = getPdf(idxCatClone->getCurrentLabel())->plotOn(frame,cmdList2);
    return retFrame ;
  }

  // If we project over the index, plot using a temporary RooAddPdf
  // using the weights from the data as coefficients

  // Build the list of indexCat components that are sliced
  idxCompSliceSet->remove(projectedVars,kTRUE,kTRUE) ;

  // Make a new expression that is the weighted sum of requested components
  RooArgList pdfCompList ;
  RooArgList wgtCompList ;
//RooAbsPdf* pdf ;
  RooRealProxy* proxy ;
  TIterator* pIter = _pdfProxyList.MakeIterator() ;
  Double_t sumWeight(0) ;
  while((proxy=(RooRealProxy*)pIter->Next())) {

    idxCatClone->setLabel(proxy->name()) ;

    // Determine if this component is the current slice (if we slice)
    Bool_t skip(kFALSE) ;
    for (const auto idxSliceCompArg : *idxCompSliceSet) {
      const auto idxSliceComp = static_cast<RooAbsCategory*>(idxSliceCompArg);
      RooAbsCategory* idxComp = (RooAbsCategory*) idxCloneSet->find(idxSliceComp->GetName()) ;
      if (idxComp->getCurrentIndex()!=idxSliceComp->getCurrentIndex()) {
        skip=kTRUE ;
        break ;
      }
    }
    if (skip) continue ;

    // Instantiate a RRV holding this pdfs weight fraction
    RooRealVar *wgtVar = new RooRealVar(proxy->name(),"coef",wTable->getFrac(proxy->name())) ;
    wgtCompList.addOwned(*wgtVar) ;
    sumWeight += wTable->getFrac(proxy->name()) ;

    // Add the PDF to list list
    pdfCompList.add(proxy->arg()) ;
  }

  TString plotVarName(GetName()) ;
  RooAddPdf *plotVar = new RooAddPdf(plotVarName,"weighted sum of RS components",pdfCompList,wgtCompList) ;

  // Fix appropriate coefficient normalization in plot function
  if (_plotCoefNormSet.getSize()>0) {
    plotVar->fixAddCoefNormalization(_plotCoefNormSet) ;
  }

  std::unique_ptr<RooAbsData> projDataTmp;
  RooArgSet projSetTmp ;
  if (projData) {

    // Construct cut string to only select projection data event that match the current slice
    TString cutString ;
    if (idxCompSliceSet->getSize()>0) {
      Bool_t first(kTRUE) ;
      for (const auto idxSliceCompArg : *idxCompSliceSet) {
        const auto idxSliceComp = static_cast<RooAbsCategory*>(idxSliceCompArg);
        if (!first) {
          cutString.Append("&&") ;
        } else {
          first=kFALSE ;
        }
        cutString.Append(Form("%s==%d",idxSliceComp->GetName(),idxSliceComp->getCurrentIndex())) ;
      }
    }

    // Make temporary projData without RooSim index category components
    RooArgSet projDataVars(*projData->get()) ;
    RooArgSet* idxCatServers = _indexCat.arg().getObservables(frame->getNormVars()) ;

    projDataVars.remove(*idxCatServers,kTRUE,kTRUE) ;

    if (idxCompSliceSet->getSize()>0) {
      projDataTmp.reset( const_cast<RooAbsData*>(projData)->reduce(projDataVars,cutString) );
    } else {
      projDataTmp.reset( const_cast<RooAbsData*>(projData)->reduce(projDataVars) );
    }



    if (projSet) {
      projSetTmp.add(*projSet) ;
      projSetTmp.remove(*idxCatServers,kTRUE,kTRUE);
    }


    delete idxCatServers ;
  }


  if (_indexCat.arg().isDerived() && idxCompSliceSet->getSize()>0) {
    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
		    << " represents a slice in index category components " << *idxCompSliceSet << endl ;

    RooArgSet* idxCompProjSet = _indexCat.arg().getObservables(frame->getNormVars()) ;
    idxCompProjSet->remove(*idxCompSliceSet,kTRUE,kTRUE) ;
    if (idxCompProjSet->getSize()>0) {
      coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
		      << " averages with data index category components " << *idxCompProjSet << endl ;
    }
    delete idxCompProjSet ;
  } else {
    coutI(Plotting) << "RooSimultaneous::plotOn(" << GetName() << ") plot on " << frame->getPlotVar()->GetName()
		    << " averages with data index category (" << _indexCat.arg().GetName() << ")" << endl ;
  }


  // Override normalization and projection dataset
  RooLinkedList cmdList2(cmdList) ;

  RooCmdArg tmp1 = RooFit::Normalization(scaleFactor*sumWeight,stype) ;
  RooCmdArg tmp2 = RooFit::ProjWData(*projDataSet,*projDataTmp) ;
  // WVE -- do not adjust normalization for asymmetry plots
  if (!cmdList.find("Asymmetry")) {
    cmdList2.Add(&tmp1) ;
  }
  cmdList2.Add(&tmp2) ;

  RooPlot* frame2 ;
  if (projSetTmp.getSize()>0) {
    // Plot temporary function
    RooCmdArg tmp3 = RooFit::Project(projSetTmp) ;
    cmdList2.Add(&tmp3) ;
    frame2 = plotVar->plotOn(frame,cmdList2) ;
  } else {
    // Plot temporary function
    frame2 = plotVar->plotOn(frame,cmdList2) ;
  }

  // Cleanup
  delete plotVar ;

  return frame2 ;
}



////////////////////////////////////////////////////////////////////////////////
/// OBSOLETE -- Retained for backward compatibility

RooPlot* RooSimultaneous::plotOn(RooPlot *frame, Option_t* drawOptions, Double_t scaleFactor,
				 ScaleType stype, const RooAbsData* projData, const RooArgSet* projSet,
				 Double_t /*precision*/, Bool_t /*shiftToZero*/, const RooArgSet* /*projDataSet*/,
				 Double_t /*rangeLo*/, Double_t /*rangeHi*/, RooCurve::WingMode /*wmode*/) const
{
  // Make command list
  RooLinkedList cmdList ;
  cmdList.Add(new RooCmdArg(RooFit::DrawOption(drawOptions))) ;
  cmdList.Add(new RooCmdArg(RooFit::Normalization(scaleFactor,stype))) ;
  if (projData) cmdList.Add(new RooCmdArg(RooFit::ProjWData(*projData))) ;
  if (projSet) cmdList.Add(new RooCmdArg(RooFit::Project(*projSet))) ;

  // Call new method
  RooPlot* ret = plotOn(frame,cmdList) ;

  // Cleanup
  cmdList.Delete() ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Interface function used by test statistics to freeze choice of observables
/// for interpretation of fraction coefficients. Needed here because a RooSimultaneous
/// works like a RooAddPdf when plotted

void RooSimultaneous::selectNormalization(const RooArgSet* normSet, Bool_t /*force*/)
{
  _plotCoefNormSet.removeAll() ;
  if (normSet) _plotCoefNormSet.add(*normSet) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Interface function used by test statistics to freeze choice of range
/// for interpretation of fraction coefficients. Needed here because a RooSimultaneous
/// works like a RooAddPdf when plotted

void RooSimultaneous::selectNormalizationRange(const char* normRange2, Bool_t /*force*/)
{
  _plotCoefNormRange = RooNameReg::ptr(normRange2) ;
}




////////////////////////////////////////////////////////////////////////////////

RooAbsGenContext* RooSimultaneous::autoGenContext(const RooArgSet &vars, const RooDataSet* prototype,
						  const RooArgSet* auxProto, Bool_t verbose, Bool_t autoBinned, const char* binnedTag) const
{
  const char* idxCatName = _indexCat.arg().GetName() ;

  if (vars.find(idxCatName) && prototype==0
      && (auxProto==0 || auxProto->getSize()==0)
      && (autoBinned || (binnedTag && strlen(binnedTag)))) {

    // Return special generator config that can also do binned generation for selected states
    return new RooSimSplitGenContext(*this,vars,verbose,autoBinned,binnedTag) ;

  } else {

    // Return regular generator config ;
    return genContext(vars,prototype,auxProto,verbose) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Return specialized generator context for simultaneous p.d.f.s

RooAbsGenContext* RooSimultaneous::genContext(const RooArgSet &vars, const RooDataSet *prototype,
					      const RooArgSet* auxProto, Bool_t verbose) const
{
  const char* idxCatName = _indexCat.arg().GetName() ;
  const RooArgSet* protoVars = prototype ? prototype->get() : 0 ;

  if (vars.find(idxCatName) || (protoVars && protoVars->find(idxCatName))) {

    // Generating index category: return special sim-context
    return new RooSimGenContext(*this,vars,prototype,auxProto,verbose) ;

  } else if (_indexCat.arg().isDerived()) {
    // Generating dependents of a derived index category

    // Determine if we none,any or all servers
    Bool_t anyServer(kFALSE), allServers(kTRUE) ;
    if (prototype) {
      TIterator* sIter = _indexCat.arg().serverIterator() ;
      RooAbsArg* server ;
      while((server=(RooAbsArg*)sIter->Next())) {
	if (prototype->get()->find(server->GetName())) {
	  anyServer=kTRUE ;
	} else {
	  allServers=kFALSE ;
	}
      }
      delete sIter ;
    } else {
      allServers=kTRUE ;
    }

    if (allServers) {
      // Use simcontext if we have all servers

      return new RooSimGenContext(*this,vars,prototype,auxProto,verbose) ;
    } else if (!allServers && anyServer) {
      // Abort if we have only part of the servers
      coutE(Plotting) << "RooSimultaneous::genContext: ERROR: prototype must include either all "
		      << " components of the RooSimultaneous index category or none " << endl ;
      return 0 ;
    }
    // Otherwise make single gencontext for current state
  }

  // Not generating index cat: return context for pdf associated with present index state
  RooRealProxy* proxy = (RooRealProxy*) _pdfProxyList.FindObject(_indexCat.arg().getCurrentLabel()) ;
  if (!proxy) {
    coutE(InputArguments) << "RooSimultaneous::genContext(" << GetName()
			  << ") ERROR: no PDF associated with current state ("
			  << _indexCat.arg().GetName() << "=" << _indexCat.arg().getCurrentLabel() << ")" << endl ;
    return 0 ;
  }
  return ((RooAbsPdf*)proxy->absArg())->genContext(vars,prototype,auxProto,verbose) ;
}




////////////////////////////////////////////////////////////////////////////////

RooDataHist* RooSimultaneous::fillDataHist(RooDataHist *hist,
                                           const RooArgSet* nset,
                                           Double_t scaleFactor,
                                           Bool_t correctForBinVolume,
                                           Bool_t showProgress) const
{
  if (RooAbsReal::fillDataHist (hist, nset, scaleFactor,
                                correctForBinVolume, showProgress) == 0)
    return 0;

  Double_t sum = 0;
  for (int i=0 ; i<hist->numEntries() ; i++) {
    hist->get(i) ;
    sum += hist->weight();
  }
  if (sum != 0) {
    for (int i=0 ; i<hist->numEntries() ; i++) {
      hist->get(i) ;
      hist->set (hist->weight() / sum);
    }
  }

  return hist;
}




////////////////////////////////////////////////////////////////////////////////
/// Special generator interface for generation of 'global observables' -- for RooStats tools

RooDataSet* RooSimultaneous::generateSimGlobal(const RooArgSet& whatVars, Int_t nEvents)
{
  // Make set with clone of variables (placeholder for output)
  RooArgSet* globClone = (RooArgSet*) whatVars.snapshot() ;

  RooDataSet* data = new RooDataSet("gensimglobal","gensimglobal",whatVars) ;

  for (Int_t i=0 ; i<nEvents ; i++) {
    for (const auto& nameIdx : indexCat()) {

      // Get pdf associated with state from simpdf
      RooAbsPdf* pdftmp = getPdf(nameIdx.first.c_str());

      // Generate only global variables defined by the pdf associated with this state
      RooArgSet* globtmp = pdftmp->getObservables(whatVars) ;
      RooDataSet* tmp = pdftmp->generate(*globtmp,1) ;

      // Transfer values to output placeholder
      *globClone = *tmp->get(0) ;

      // Cleanup
      delete globtmp ;
      delete tmp ;
    }
    data->add(*globClone) ;
  }

  delete globClone ;
  return data ;
}
