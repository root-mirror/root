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

//////////////////////////////////////////////////////////////////////////////
///
/// \class RooSimPdfBuilder
///
/// <b>This tool has now been superseded by RooSimWSTool</b>
///
///  <p>
///    <tt>RooSimPdfBuilder</tt> is a powerful tool to build <tt>RooSimultaneous</tt>
///    PDFs that are defined in terms component PDFs that are identical in
///    structure, but have different parameters.
///  </p>
///
///  <h2>Example</h2>
///
///  <p>
///    The following example demonstrates the essence of <tt>RooSimPdfBuilder</tt>:
///    Given a dataset D with a <tt>RooRealVar X</tt> and a <tt>RooCategory C</tt> that has
///    state C1 and C2.
///    <ul>
///    <li> We want to fit the distribution of <tt>X</tt> with a Gaussian+ArgusBG PDF,
///    <li> We want to fit the data subsets <tt>D(C==C1)</tt> and <tt>D(C==C2)</tt> separately and simultaneously.
///    <li> The PDFs to fit data subsets D_C1 and D_C2 are identical except for
///         <ul>
///         <li> the kappa parameter of the ArgusBG PDF and
///         <li> the sigma parameter of the gaussian PDF
///         </ul>
///         where each PDF will have its own copy of the parameter
///    </ul>
///  </p>
///  <p>
///    Coding this example directly with RooFit classes gives
///    (we assume dataset D and variables C and X have been declared previously)
///  </p>
///  <pre>
/// RooRealVar m("m","mean of gaussian",-10,10) ;
/// RooRealVar s_C1("s_C1","sigma of gaussian C1",0,20) ;
/// RooRealVar s_C2("s_C2","sigma of gaussian C2",0,20) ;
/// RooGaussian gauss_C1("gauss_C1","gaussian C1",X,m,s_C1) ;
/// RooGaussian gauss_C2("gauss_C2","gaussian C2",X,m,s_C2) ;
///
/// RooRealVar k_C1("k_C1","ArgusBG kappa parameter C1",-50,0) ;
/// RooRealVar k_C2("k_C2","ArgusBG kappa parameter C2",-50,0) ;
/// RooRealVar xm("xm","ArgusBG cutoff point",5.29) ;
/// RooArgusBG argus_C1("argus_C1","argus background C1",X,k_C1,xm) ;
/// RooArgusBG argus_C2("argus_C2","argus background C2",X,k_C2,xm) ;
///
/// RooRealVar gfrac("gfrac","fraction of gaussian",0.,1.) ;
/// RooAddPdf pdf_C1("pdf_C1","gauss+argus_C1",RooArgList(gauss_C1,argus_C1),gfrac) ;
/// RooAddPdf pdf_C2("pdf_C2","gauss+argus_C2",RooArgList(gauss_C2,argus_C2),gfrac) ;
///
/// RooSimultaneous simPdf("simPdf","simPdf",C) ;
/// simPdf.addPdf(pdf_C1,"C1") ;
/// simPdf.addPdf(pdf_C2,"C2") ;
///  </pre>
///  <p>
///    Coding this example with RooSimPdfBuilder gives
///  </p>
///  <pre>
/// RooRealVar m("m","mean of gaussian",-10,10) ;
/// RooRealVar s("s","sigma of gaussian",0,20) ;
/// RooGaussian gauss("gauss","gaussian",X,m,s) ;
///
/// RooRealVar k("k","ArgusBG kappa parameter",-50,0) ;
/// RooRealVar xm("xm","ArgusBG cutoff point",5.29) ;
/// RooArgusBG argus("argus","argus background",X,k,xm) ;
///
/// RooRealVar gfrac("gfrac","fraction of gaussian",0.,1.) ;
/// RooAddPdf pdf("pdf","gauss+argus",RooArgList(gauss,argus),gfrac) ;
///
/// RooSimPdfBuilder builder(pdf) ;
/// RooArgSet* config = builder.createProtoBuildConfig() ;
/// (*config)["physModels"] = "pdf" ;      /// Name of the PDF we are going to work with
/// (*config)["splitCats"]  = "C" ;        /// Category used to differentiate sub-datasets
/// (*config)["pdf"]        = "C : k,s" ;  /// Prescription to taylor PDF parameters k and s
///                                        /// for each data subset designated by C states
/// RooSimultaneous* simPdf = builder.buildPdf(*config,&D) ;
///  </pre>
///  <p>
///    The above snippet of code demonstrates the concept of <tt>RooSimPdfBuilder</tt>:
///    the user defines a single <i>'prototype' PDF</i> that defines the structure of all
///    PDF components of the <tt>RooSimultaneous</tt> PDF to be built. <tt>RooSimPdfBuilder</tt>
///    then takes this prototype and replicates it as a component
///    PDF for each state of the C index category.
///  </p>
///  <p>
///    In the above example <tt>RooSimPdfBuilder</tt>
///    will first replicate <tt>k</tt> and <tt>s</tt> into
///    <tt>k_C1,k_C2</tt> and <tt>s_C1,s_C2</tt>, as prescribed in the
///    configuration. Then it will recursively replicate all PDF nodes that depend on
///    the 'split' parameter nodes: <tt>gauss</tt> into <tt>gauss_C1,C2</tt>, <tt>argus</tt>
///    into <tt>argus_C1,C2</tt> and finally <tt>pdf</tt> into <tt>pdf_C1,pdf_C2</tt>.
///     When PDFs for all states of C have been replicated
///    they are assembled into a <tt>RooSimultaneous</tt> PDF, which is returned by the <tt>buildPdf()</tt>
///    method.
///  </p>
///  <p>
///    Although in this very simple example the use of <tt>RooSimPdfBuilder</tt> doesn't
///    reduce the amount of code much, it is already easier to read and maintain
///    because there is no duplicate code. As the complexity of the <tt>RooSimultaneous</tt>
///    to be built increases, the advantages of <tt>RooSimPdfBuilder</tt> will become more and
///    more apparent.
///  </p>
///
///
///  <h2>Builder configuration rules for a single prototype PDF</h2>
///  <p>
///    Each builder configuration needs at minumum two lines, <tt>physModels</tt> and <tt>splitCats</tt>, which identify
///    the ingredients of the build. In this section we only explain the building rules for
///    builds from a single prototype PDF. In that case the <tt>physModels</tt> line always reads
///  </p>
///  <pre>
///  physModels = {pdfName}
///  </pre>
///  <p>
///    The second line, <tt>splitCats</tt>, indicates which categories are going to be used to
///    differentiate the various subsets of the 'master' input data set. You can enter
///    a single category here, or multiple if necessary:
///  </p>
///  <pre>
/// splitCats = {catName} [{catName} ...]
///  </pre>
///  <p>
///    All listed splitcats must be <tt>RooCategories</tt> that appear in the dataset provided to
///    <tt>RooSimPdfBuilder::buildPdf()</tt>
///  </p>
///  <p>
///    The parameter splitting prescriptions, the essence of each build configuration
///    can be supplied in a third line carrying the name of the pdf listed in <tt>physModels</tt>
///  </p>
///  <pre>
/// pdfName = {splitCat} : {parameter} [,{parameter},....]
///  </pre>
///  <p>
///    Each pdf can have only one line with splitting rules, but multiple rules can be
///    supplied in each line, e.g.
///  </p>
///  <pre>
/// pdfName = {splitCat} : {parameter} [,{parameter},....]
///           {splitCat} : {parameter} [,{parameter},....]
///  </pre>
///  <p>
///    Conversely, each parameter can only have one splitting prescription, but it may be split
///    by multiple categories, e.g.
///  </p>
///  <pre>
/// pdfName = {splitCat1},{splitCat2} : {parameter}
///  </pre>
///  <p>
///    instructs <tt>RooSimPdfBuilder</tt> to build a <tt>RooSuperCategory</tt>
///    of <tt>{splitCat1}</tt> and <tt>{splitCat2}</tt>
///    and split <tt>{parameter}</tt> with that <tt>RooSuperCategory</tt>
///  </p>
///  <p>
///    Here is an example of a builder configuration that uses several of the options discussed
///    above:
///  </p>
///  <pre>
///   physModels = pdf
///   splitCats  = tagCat runBlock
///   pdf        = tagCat          : signalRes,bkgRes
///                runBlock        : fudgeFactor
///                tagCat,runBlock : kludgeParam
///  </pre>
///
///  <h2>How to enter configuration data</h2>
///
///  <p>
///    The prototype builder configuration returned by
///    <tt>RooSimPdfBuilder::createProtoBuildConfig()</tt> is a pointer to a <tt>RooArgSet</tt> filled with
///    initially blank <tt>RooStringVars</tt> named <tt>physModels,splitCats</tt> and one additional for each
///    PDF supplied to the <tt>RooSimPdfBuilders</tt> constructor (with the same name)
///  </p>
///  <p>
///    In macro code, the easiest way to assign new values to these <tt>RooStringVars</tt>
///    is to use <tt>RooArgSet</tt>s array operator and the <tt>RooStringVar</tt>s assignment operator, e.g.
///  </p>
///  <pre>
/// (*config)["physModels"] = "Blah" ;
///  </pre>
///  <p>
///    To enter multiple splitting rules simply separate consecutive rules by whitespace
///    (not newlines), e.g.
///  </p>
///  <pre>
/// (*config)["physModels"] = "Blah " /// << note trailing space here
///                          "Blah 2" ;
///  </pre>
///  <p>
///    In this example, the C++ compiler will concatenate the two string literals (without inserting
///    any whitespace), so the extra space after 'Blah' is important here.
///  </p>
///  <p>
///    Alternatively, you can read the configuration from an ASCII file, as you can
///    for any <tt>RooArgSet</tt> using <tt>RooArgSet::readFromFile()</tt>. In that case the ASCII file
///    can follow the syntax of the examples above and the '<tt>\\</tt>' line continuation
///    sequence can be used to fold a long splitting rule over multiple lines.
///  </p>
///  <pre>
/// RooArgSet* config = builder.createProtoBuildConfig() ;
/// config->readFromFile("config.txt") ;
///
/// --- config.txt ----------------
/// physModels = pdf
/// splitCats  = tagCat
/// pdf        = tagCat : bogusPar
/// -------------------------------
///  </pre>
///
///
///  <h2>Working with multiple prototype PDFs</h2>
///  <p>
///    It is also possible to build a <tt>RooSimultaneous</tt> PDF from multiple PDF prototypes.
///    This is appropriate for cases where the input prototype PDF would otherwise be
///    a <tt>RooSimultaneous</tt> PDF by itself. In such cases we don't feed a single
///    <tt>RooSimultaneous</tt> PDF into <tt>RooSimPdfBuilder</tt>, instead we feed it its ingredients and
///    add a prescription to the builder configuration that corresponds to the
///    PDF-category state mapping of the prototype <tt>RooSimultaneous</tt>.
///  </p>
///  <p>
///    The constructor of the <tt>RooSimPdfBuilder</tt> will look as follows:
///  </p>
///  <pre>
///  RooSimPdfBuilder builder(RooArgSet(pdfA,pdfB,...)) ;
///  </pre>
///  <p>
///    The <tt>physModels</tt> line is now expanded to carry the pdf->state mapping information
///    that the prototype <tt>RooSimultaneous</tt> would have. I.e.
///  </p>
///  <pre>
/// physModels = mode : pdfA=modeA  pdfB=modeB
///  </pre>
///  <p>
///    is equivalent to a prototype <tt>RooSimultaneous</tt> constructed as
///  </p>
///  <pre>
/// RooSimultanous simPdf("simPdf","simPdf",mode);
/// simPdf.addPdf(pdfA,"modeA") ;
/// simPdf.addPdf(pdfB,"modeB") ;
///  </pre>
///  <p>
///    The rest of the builder configuration works the same, except that
///    each prototype PDF now has its own set of splitting rules, e.g.
///  </p>
///  <pre>
/// physModels = mode : pdfA=modeA  pdfB=modeB
/// splitCats  = tagCat
/// pdfA       = tagCat : bogusPar
/// pdfB       = tagCat : fudgeFactor
///  </pre>
///  <p>
///    Please note that
///    <ul>
///    <li> The master index category ('mode' above) doesn't have to be listed in
///         <tt>splitCats</tt>, this is implicit.
///
///    <li> The number of splitting prescriptions goes by the
///          number of prototype PDFs and not by the number of states of the
///          master index category (mode in the above and below example).
///  </ul>
///
///  In the following case:
///</p>
///  <pre>
///    physModels = mode : pdfA=modeA  pdfB=modeB  pdfA=modeC  pdfB=modeD
///  </pre>
///  <p>
///    there are still only 2 sets of splitting rules: one for <tt>pdfA</tt> and one
///    for <tt>pdfB</tt>. However, you <i>can</i> differentiate between <tt>modeA</tt> and <tt>modeC</tt> in
///    the above example. The technique is to use <tt>mode</tt> as splitting category, e.g.
///  </p>
///  <pre>
///    physModels = mode : pdfA=modeA  pdfB=modeB  pdfA=modeC  pdfB=modeD
///    splitCats = tagCat
///    pdfA      = tagCat : bogusPar
///                mode   : funnyPar
///    pdfB      = mode   : kludgeFactor
///  </pre>
///  <p>
///    will result in an individual set of <tt>funnyPar</tt> parameters for <tt>modeA</tt> and <tt>modeC</tt>
///    labeled <tt>funnyPar_modeA</tt> and <tt>funnyPar_modeB</tt> and an individual set of
///    kludgeFactor parameters for <tt>pdfB</tt>, <tt>kludgeFactor_modeB</tt> and <tt>kludgeFactor_modeD</tt>.
///    Please note that for splits in the master index category (mode) only the
///    applicable states are built (A,C for <tt>pdfA</tt>, B,D for <tt>pdfB</tt>)
///  </p>
///
///
///  <h2>Advanced options</h2>
///
///  <h4>Partial splits</h4>
///  <p>
///    You can request to limit the list of states of each splitCat that
///    will be considered in the build. This limitation is requested in the
///    each build as follows:
///  </p>
///  <pre>
/// splitCats = tagCat(Lep,Kao) RunBlock(Run1)
///  </pre>
///  <p>
///    In this example the splitting of <tt>tagCat</tt> is limited to states <tt>Lep,Kao</tt>
///    and the splitting of <tt>runBlock</tt> is limited to <tt>Run1</tt>. The splits apply
///    globally to each build, i.e. every parameter split requested in this
///    build will be limited according to these specifications.
///  </p>
///  <p>
///    NB: Partial builds have no pdf associated with the unbuilt states of the
///    limited splits. Running such a pdf on a dataset that contains data with
///    unbuilt states will result in this data being ignored completely.
///  </p>
///
///
///  <h4>Non-trivial splits</h4>
///  <p>
///    It is possible to make non-trivial parameter splits with <tt>RooSimPdfBuilder</tt>.
///    Trivial splits are considered simple splits in one (fundamental) category
///    in the dataset or a split in a <tt>RooSuperCategory</tt> 'product' of multiple
///    fundamental categories in the dataset. Non-trivial splits can be performed
///    using an intermediate 'category function' (<tt>RooMappedCategory,
///    RooGenericCategory,RooThresholdCategory</tt> etc), i.e. any <tt>RooAbsCategory</tt>
///    derived objects that calculates its output as function of one or more
///    input <tt>RooRealVars</tt> and/or <tt>RooCategories</tt>.
///  </p>
///  <p>
///    Such 'function categories' objects must be constructed by the user prior
///    to building the PDF. In the <tt>RooSimPdfBuilder::buildPdf()</tt> function these
///    objects can be passed in an optional <tt>RooArgSet</tt> called 'auxiliary categories':
///  </p>
///  <pre>
///   const <tt>RooSimultaneous</tt>* buildPdf(const RooArgSet& buildConfig, const RooAbsData* dataSet,
///                                   const RooArgSet& auxSplitCats, Bool_t verbose=kFALSE) {
///                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
///  </pre>
///  <p>
///    Objects passed in this argset can subsequently be used in the build configuration, e.g.
///  </p>
///  <pre>
/// RooMappedCategory tagMap("tagMap","Mapped tagging category",tagCat,"CutBased") ;
/// tagMap.map("Lep","CutBased") ;
/// tagMap.map("Kao","CutBased") ;
/// tagMap.map("NT*","NeuralNet") ;
/// ...
/// builder.buildPdf(config,D,tagMap) ;
///                          ^^^^^^
///<Contents of config>
///   physModels = pdf
///   splitCats  = tagCat runBlock
///   pdf        = tagCat          : signalRes
///                tagMap          : fudgeFactor
///                ^^^^^^
///  </pre>
///  <p>
///    In the above example <tt>signalRes</tt> will be split in <tt>signalRes_Kao,signalRes_Lep,
///    signalRes_NT1,signalRes_NT2</tt>, while <tt>fudgeFactor</tt> will be split in <tt>fudgeFactor_CutBased</tt>
///    and <tt>fudgeFactor_NeuralNet</tt>.
///  </p>
///  <p>
///    Category functions passed in the auxSplitCats <tt>RooArgSet</tt> can be used regularly
///    in the splitting configuration. They should not be listed in <tt>splitCats</tt>,
///    but must be able to be expressed <i>completely</i> in terms of the <tt>splitCats</tt> that
///    are listed.
///  </p>
///
///
///  <h4>Multiple connected builds</h4>
///  <p>
///    Sometimes you want to build multiple PDFs for independent consecutive fits
///    that share some of their parameters. For example, we have two prototype PDFs
///    <tt>pdfA(x;p,q)</tt> and <tt>pdfB(x;p,r)</tt> that have a common parameter <tt>p</tt>.
///    We want to build a <tt>RooSimultaneous</tt> for both <tt>pdfA</tt> and <tt>B</tt>,
///    which involves a split of parameter <tt>p</tt> and we would like to build the
///    simultaneous pdfs <tt>simA</tt> and <tt>simB</tt> such that still share their (now split) parameters
///    <tt>p_XXX</tt>. This is accomplished by letting a single instance of <tt>RooSimPdfBuilder</tt> handle
///    the builds of both <tt>pdfA</tt> and <tt>pdfB</tt>, as illustrated in this example:
///  </p>
///  <pre>
/// RooSimPdfBuilder builder(RooArgSet(pdfA,pdfB)) ;
///
/// RooArgSet* configA = builder.createProtoBuildConfig() ;
/// (*configA)["physModels"] = "pdfA" ;
/// (*configA)["splitCats"]  = "C" ;
/// (*configA)["pdf"]        = "C : p" ;
/// RooSimultaneous* simA = builder.buildPdf(*configA,&D) ;
///
/// RooArgSet* configB = builder.createProtoBuildConfig() ;
/// (*configA)["physModels"] = "pdfB" ;
/// (*configA)["splitCats"]  = "C" ;
/// (*configA)["pdf"]        = "C : p" ;
/// RooSimultaneous* simB = builder.buildPdf(*configB,&D) ;
///  </pre>
///
///  <h2>Ownership of constructed PDFs</h2>
///  <p>
///    The <tt>RooSimPdfBuilder</tt> instance owns all the objects it creates, including the top-level
///    <tt>RooSimultaneous</tt> returned by <tt>buildPdf()</tt>. Therefore the builder instance should
///    exist as long as the constructed PDFs needs to exist.
///  </p>
///
///
///


#include "RooFit.h"

#include <cstring>
#include "strtok.h"

#ifndef _WIN32
#include <strings.h>
#endif

#include "Riostream.h"
#include "RooSimPdfBuilder.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsCategory.h"
#include "RooCategory.h"
#include "RooStringVar.h"
#include "RooMappedCategory.h"
#include "RooRealIntegral.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooLinearVar.h"
#include "RooTruthModel.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooCustomizer.h"
#include "RooThresholdCategory.h"
#include "RooMultiCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooTrace.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooMsgService.h"
#include "RooHelpers.h"

using namespace std ;

ClassImp(RooSimPdfBuilder);
;



////////////////////////////////////////////////////////////////////////////////

RooSimPdfBuilder::RooSimPdfBuilder(const RooArgSet& protoPdfSet) :
  _protoPdfSet(protoPdfSet)
{
  _compSplitCatSet.useHashMapForFind(true);
  _splitNodeList.useHashMapForFind(true);
  _splitNodeListOwned.useHashMapForFind(true);
}




////////////////////////////////////////////////////////////////////////////////
/// Make RooArgSet of configuration objects

RooArgSet* RooSimPdfBuilder::createProtoBuildConfig()
{
  RooArgSet* buildConfig = new RooArgSet ;
  buildConfig->addOwned(* new RooStringVar("physModels","List and mapping of physics models to include in build","",4096)) ;
  buildConfig->addOwned(* new RooStringVar("splitCats","List of categories used for splitting","",1024)) ;

  TIterator* iter = _protoPdfSet.createIterator() ;
  RooAbsPdf* proto ;
  while ((proto=(RooAbsPdf*)iter->Next())) {
    buildConfig->addOwned(* new RooStringVar(proto->GetName(),proto->GetName(),"",4096)) ;
  }
  delete iter ;

  return buildConfig ;
}



////////////////////////////////////////////////////////////////////////////////

void RooSimPdfBuilder::addSpecializations(const RooArgSet& specSet) 
{
  _splitNodeList.add(specSet) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Initialize needed components

RooSimultaneous* RooSimPdfBuilder::buildPdf(const RooArgSet& buildConfig, const RooArgSet& dependents,
    const RooArgSet* auxSplitCats, Bool_t verbose)
{
  const char* spaceChars = " \t" ;

  // Retrieve physics index category
  Int_t buflen = strlen(((RooStringVar*)buildConfig.find("physModels"))->getVal())+1 ;
  char *buf = new char[buflen] ;

  strlcpy(buf,((RooStringVar*)buildConfig.find("physModels"))->getVal(),buflen) ;
  RooAbsCategoryLValue* physCat(0) ;
  if (strstr(buf," : ")) {
    const char* physCatName = strtok(buf,spaceChars) ;
    physCat = dynamic_cast<RooAbsCategoryLValue*>(dependents.find(physCatName)) ;
    if (!physCat) {
      coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR physics index category " << physCatName 
          << " not found in dataset variables" << endl ;
      delete[] buf ;
      return 0 ;      
    }
    coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: category indexing physics model: " << physCatName << endl ;
  }

  // Create list of physics models to be built
  char *physName ;
  RooArgSet physModelSet ;
  if (physCat) {
    // Absorb colon token
    strtok(0,spaceChars) ;
    physName = strtok(0,spaceChars) ;
  } else {
    physName = strtok(buf,spaceChars) ;
  }

  if (!physName) {
    coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR: No models specified, nothing to do!" << endl ;
    delete[] buf ;
    return 0 ;
  }

  Bool_t first(kTRUE) ;
  RooArgSet stateMap ;
  while(physName) {

    char *stateName(0) ;

    // physName may be <state>=<pdfName> or just <pdfName> is state and pdf have identical names
    if (strchr(physName,'=')) {
      // Must have a physics category for mapping to make sense
      if (!physCat) {
        coutW(ObjectHandling) << "RooSimPdfBuilder::buildPdf: WARNING: without physCat specification "
            << "<physCatState>=<pdfProtoName> association is meaningless" << endl ;
      }
      stateName = physName ;
      physName = strchr(stateName,'=') ;
      if (physName) {
        *(physName++) = 0 ;
      }
    } else {
      stateName = physName ;
    }

    RooAbsPdf* physModel = (RooAbsPdf*) (physName ? _protoPdfSet.find(physName) : 0 );
    if (!physModel) {
      coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR requested physics model " 
          << (physName?physName:"(null)") << " is not defined" << endl ;
      delete[] buf ;
      return 0 ;
    }    

    // Check if state mapping has already been defined
    if (stateMap.find(stateName)) {
      coutW(InputArguments) << "RooSimPdfBuilder::buildPdf: WARNING: multiple PDFs specified for state " 
          << stateName << ", only first will be used" << endl ;
      continue ;
    }

    // Add pdf to list of models to be processed
    physModelSet.add(*physModel,kTRUE) ; // silence duplicate insertion warnings

    // Store state->pdf mapping    
    stateMap.addOwned(* new RooStringVar(stateName,stateName,physName)) ;

    // Continue with next mapping
    physName = strtok(0,spaceChars) ;
    if (first) {
      first = kFALSE ;
    } else if (physCat==0) {
      coutW(InputArguments) << "RooSimPdfBuilder::buildPdf: WARNING: without physCat specification, only the first model will be used" << endl ;
      break ;
    }
  }
  coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: list of physics models " << physModelSet << endl ;



  // Create list of dataset categories to be used in splitting
  TList splitStateList ;
  RooArgSet splitCatSet ;

  delete[] buf ; 
  buflen = strlen(((RooStringVar*)buildConfig.find("splitCats"))->getVal())+1 ;
  buf = new char[buflen] ;
  strlcpy(buf,((RooStringVar*)buildConfig.find("splitCats"))->getVal(),buflen) ;

  char *catName = strtok(buf,spaceChars) ;
  char *stateList(0) ;
  while(catName) {

    // Chop off optional list of selected states
    char* tokenPtr(0) ;
    if (strchr(catName,'(')) {

      catName = R__STRTOK_R(catName,"(",&tokenPtr) ;
      stateList = R__STRTOK_R(0,")",&tokenPtr) ;

    } else {
      stateList = 0 ;
    }

    RooCategory* splitCat = catName ? dynamic_cast<RooCategory*>(dependents.find(catName)) : 0 ;
    if (!splitCat) {
      coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR requested split category " << (catName?catName:"(null)") 
			        << " is not a RooCategory in the dataset" << endl ;
      delete[] buf ;
      return 0 ;
    }
    splitCatSet.add(*splitCat) ;

    // Process optional state list
    if (stateList) {
      coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: splitting of category " << catName 
          << " restricted to states (" << stateList << ")" << endl ;

      // Create list named after this splitCat holding its selected states
      TList* slist = new TList ;
      slist->SetName(catName) ;
      splitStateList.Add(slist) ;

      char* stateLabel = R__STRTOK_R(stateList,",",&tokenPtr) ;

      while(stateLabel) {
        // Lookup state label and require it exists
        const RooCatType* type = splitCat->lookupType(stateLabel) ;
        if (!type) {
          coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR splitCat " << splitCat->GetName()
				    << " doesn't have a state named " << stateLabel << endl ;
          splitStateList.Delete() ;
          delete[] buf ;
          return 0 ;
        }
        slist->Add((TObject*)type) ;

        stateLabel = R__STRTOK_R(0,",",&tokenPtr) ;
      }
    }

    catName = strtok(0,spaceChars) ;
  }
  if (physCat) splitCatSet.add(*physCat) ;
  RooSuperCategory masterSplitCat("masterSplitCat","Master splitting category",splitCatSet) ;

  coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: list of splitting categories " << splitCatSet << endl ;

  // Clone auxiliary split cats and attach to splitCatSet
  RooArgSet auxSplitSet ;
  RooArgSet* auxSplitCloneSet(0) ;
  if (auxSplitCats) {
    // Deep clone auxililary split cats
    auxSplitCloneSet = (RooArgSet*) auxSplitCats->snapshot(kTRUE) ;
    if (!auxSplitCloneSet) {
      coutE(InputArguments) << "RooSimPdfBuilder::buildPdf(" << GetName() << ") Couldn't deep-clone set auxiliary splitcats, abort." << endl ;
      delete[] buf ;
      return 0 ;
    }

    for (const auto arg : *auxSplitCats) {
      // Find counterpart in cloned set
      RooAbsArg* aux = auxSplitCats->find(arg->GetName()) ;

      // Check that there is no fundamental splitCat in the dataset with the bane of the auxiliary split
      if (splitCatSet.find(aux->GetName())) {
        coutW(InputArguments) << "RooSimPdfBuilder::buildPdf: WARNING: dataset contains a fundamental splitting category " << endl
            << " with the same name as an auxiliary split function (" << aux->GetName() << "). " << endl
            << " Auxiliary split function will be ignored" << endl ;
        continue ;
      }

      // Check that all servers of this aux cat are contained in splitCatSet
      RooArgSet* parSet = aux->getParameters(splitCatSet) ;
      if (parSet->getSize()>0) {
        coutW(InputArguments) << "RooSimPdfBuilder::buildPdf: WARNING: ignoring auxiliary category " << aux->GetName()
			          << " because it has servers that are not listed in splitCatSet: " << *parSet << endl ;
        delete parSet ;
        continue ;
      }

      // Redirect servers to splitCatSet
      aux->recursiveRedirectServers(splitCatSet) ;

      // Add top level nodes to auxSplitSet
      auxSplitSet.add(*aux) ;
    }

    coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: list of auxiliary splitting categories " << auxSplitSet << endl ;
  }

  delete[] buf;

  TList* customizerList = new TList ;

  // Loop over requested physics models and build components
  for (const auto arg : physModelSet) {
    auto physModel = static_cast<const RooAbsPdf*>(arg);

    coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: processing physics model " << physModel->GetName() << endl ;

    RooCustomizer* physCustomizer = new RooCustomizer(*physModel,masterSplitCat,_splitNodeList) ;
    customizerList->Add(physCustomizer) ;

    // Parse the splitting rules for this physics model
    RooStringVar* ruleStr = (RooStringVar*) buildConfig.find(physModel->GetName()) ;
    if (ruleStr) {

      enum Mode { SplitCat, Colon, ParamList } ;
      Mode mode(SplitCat) ;

      const char* splitCatName ;
      RooAbsCategory* splitCat(0) ;

      for (const auto& token : RooHelpers::tokenise(ruleStr->getVal(), spaceChars)) {

        switch (mode) {
        case SplitCat:
        {
          splitCatName = token.data();

          if (token.find(',') != std::string::npos) {
            // Composite splitting category

            // Check if already instantiated
            splitCat = (RooAbsCategory*) _compSplitCatSet.find(splitCatName) ;
            TString origCompCatName(splitCatName) ;
            if (!splitCat) {
              // Build now

              RooArgSet compCatSet ;
              for (const auto& catName2 : RooHelpers::tokenise(token, ",")) {
                RooAbsArg* cat = splitCatSet.find(catName2.data()) ;

                // If not, check if it is an auxiliary splitcat
                if (!cat) {
                  cat = (RooAbsCategory*) auxSplitSet.find(catName2.data()) ;
                }

                if (!cat) {
                  coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR " << catName2
                      << " not found in the primary or auxilary splitcat list" << endl ;
                  customizerList->Delete() ;
                  delete customizerList ;

                  splitStateList.Delete() ;
                  return 0 ;
                }
                compCatSet.add(*cat) ;
              }


              // check that any auxSplitCats in compCatSet do not depend on any other
              // fundamental or auxiliary splitting categories in the composite set.
              for (const auto theArg : compCatSet) {
                RooArgSet tmp(compCatSet) ;
                tmp.remove(*theArg) ;
                if (theArg->dependsOnValue(tmp)) {
                  coutE(InputArguments) << "RooSimPdfBuilder::buildPDF: ERROR: Ill defined split: auxiliary splitting category " << theArg->GetName()
					            << " used in composite split " << compCatSet << " depends on one or more of the other splitting categories in the composite split" << endl ;

                  // Cleanup and axit
                  customizerList->Delete() ;
                  delete customizerList ;
                  splitStateList.Delete() ;
                  return 0 ;
                }
              }

              splitCat = new RooMultiCategory(origCompCatName,origCompCatName,compCatSet) ;
              _compSplitCatSet.addOwned(*splitCat) ;
              //cout << "composite splitcat: " << splitCat->GetName() ;
            }
          } else {
            // Simple splitting category

            // First see if it is a simple splitting category
            splitCat = (RooAbsCategory*) splitCatSet.find(splitCatName) ;

            // If not, check if it is an auxiliary splitcat
            if (!splitCat) {
              splitCat = (RooAbsCategory*) auxSplitSet.find(splitCatName) ;
            }

            if (!splitCat) {
              coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR splitting category "
                  << splitCatName << " not found in the primary or auxiliary splitcat list" << endl ;
              customizerList->Delete() ;
              delete customizerList ;
              splitStateList.Delete() ;
              return 0 ;
            }
          }

          mode = Colon ;
          break ;
        }
        case Colon:
        {
          if (token != ":") {
            coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR in parsing, expected ':' after "
                << splitCat << ", found " << token << endl ;
            customizerList->Delete() ;
            delete customizerList ;
            splitStateList.Delete() ;
            return 0 ;
          }
          mode = ParamList ;
          break ;
        }
        case ParamList:
        {
          // Verify the validity of the parameter list and build the corresponding argset
          RooArgSet splitParamList ;
          RooArgSet* paramList = physModel->getParameters(dependents) ;

          // wve -- add nodes to parameter list
          RooArgSet* compList = physModel->getComponents() ;
          paramList->add(*compList) ;
          delete compList ;

          const bool lastCharIsComma = (token[token.size()-1]==',') ;

          for (const auto& paramName : RooHelpers::tokenise(token, ",")) {
            // Check for fractional split option 'param_name[remainder_state]'
            std::string remainderState;
            if (const auto pos = paramName.find('[') != std::string::npos) {
              const auto posEnd = paramName.find(']');
              remainderState = paramName.substr(pos+1, posEnd - pos - 1);
            }


            // If fractional split is specified, check that remainder state is a valid state of this split cat
            if (!remainderState.empty()) {
              if (!splitCat->hasLabel(remainderState)) {
                coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR fraction split of parameter "
                    << paramName << " has invalid remainder state name: " << remainderState << endl ;
                delete paramList ;
                customizerList->Delete() ;
                delete customizerList ;
                splitStateList.Delete() ;
                return 0 ;
              }
            }

            RooAbsArg* param = paramList->find(paramName.data()) ;
            if (!param) {
              coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR " << paramName
                  << " is not a parameter of physics model " << physModel->GetName() << endl ;
              delete paramList ;
              customizerList->Delete() ;
              delete customizerList ;
              splitStateList.Delete() ;
              return 0 ;
            }
            splitParamList.add(*param) ;

            // Build split leaf of fraction splits here
            if (!remainderState.empty()) {

              // Check if we are splitting a real-valued parameter
              if (!dynamic_cast<RooAbsReal*>(param)) {
                coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR fraction split requested of non-real valued parameter "
                    << param->GetName() << endl ;
                delete paramList ;
                customizerList->Delete() ;
                delete customizerList ;
                splitStateList.Delete() ;
                return 0 ;
              }

              // Check if we are doing a restricted build
              TList* remStateSplitList = static_cast<TList*>(splitStateList.FindObject(splitCat->GetName())) ;

              // If so, check if remainder state is actually being built.
              if (remStateSplitList && !remStateSplitList->FindObject(remainderState.data())) {
                coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR " << paramName
                    << " remainder state " << remainderState << " in parameter split "
                    << param->GetName() << " is not actually being built" << endl ;
                delete paramList ;
                customizerList->Delete() ;
                delete customizerList ;
                splitStateList.Delete() ;
                return 0 ;
              }

              TIterator* iter = splitCat->typeIterator() ;
              RooCatType* type ;
              RooArgList fracLeafList ;
              TString formExpr("1") ;
              Int_t i(0) ;

              while((type=(RooCatType*)iter->Next())) {

                // Skip remainder state
                if (!TString(type->GetName()).CompareTo(remainderState)) continue ;

                // If restricted build is requested, skip states of splitcat that are not built
                if (remStateSplitList && !remStateSplitList->FindObject(type->GetName())) {
                  continue ;
                }

                // Construct name of split leaf
                TString splitLeafName(param->GetName()) ;
                splitLeafName.Append("_") ;
                splitLeafName.Append(type->GetName()) ;

                // Check if split leaf already exists
                RooAbsArg* splitLeaf = _splitNodeList.find(splitLeafName) ;
                if (!splitLeaf) {
                  // If not create it now
                  splitLeaf = (RooAbsArg*) param->clone(splitLeafName) ;
                  _splitNodeList.add(*splitLeaf) ;
                  _splitNodeListOwned.addOwned(*splitLeaf) ;
                }
                fracLeafList.add(*splitLeaf) ;
                formExpr.Append(Form("-@%d",i++)) ;
              }
              delete iter ;

              // Construct RooFormulaVar expresssing remainder of fraction
              TString remLeafName(param->GetName()) ;
              remLeafName.Append("_") ;
              remLeafName.Append(remainderState) ;

              // Check if no specialization was already specified for remainder state
              if (!_splitNodeList.find(remLeafName)) {
                RooAbsArg* remLeaf = new RooFormulaVar(remLeafName,formExpr,fracLeafList) ;
                _splitNodeList.add(*remLeaf) ;
                _splitNodeListOwned.addOwned(*remLeaf) ;
                coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: creating remainder fraction formula for " << remainderState
                    << " specialization of split parameter " << param->GetName() << " " << formExpr << endl ;
              }
            }
          }

          // Add the rule to the appropriate customizer ;
          physCustomizer->splitArgs(splitParamList,*splitCat) ;

          delete paramList ;

          if (!lastCharIsComma) mode = SplitCat ;
          break ;
        }
        }
      }
      if (mode!=SplitCat) {
        coutE(InputArguments) << "RooSimPdfBuilder::buildPdf: ERROR in parsing, expected "
            << (mode==Colon?":":"parameter list") << " in " << ruleStr->getVal() << endl ;
      }

      //RooArgSet* paramSet = physModel->getParameters(dependents) ;
    } else {
      coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: no splitting rules for " << physModel->GetName() << endl ;
    }    
  }

  coutI(ObjectHandling)  << "RooSimPdfBuilder::buildPdf: configured customizers for all physics models" << endl ;
  if (oodologI((TObject*)0,ObjectHandling)) {
    customizerList->Print() ;
  }

  // Create fit category from physCat and splitCatList ;
  RooArgSet fitCatList ;
  if (physCat) fitCatList.add(*physCat) ;
  fitCatList.add(splitCatSet) ;
  RooSuperCategory *fitCat = new RooSuperCategory("fitCat","fitCat",fitCatList) ;

  // Create master PDF 
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",*fitCat) ;

  // Add component PDFs to master PDF
  TIterator* fcIter = fitCat->typeIterator() ;

  RooCatType* fcState ;  
  while((fcState=(RooCatType*)fcIter->Next())) {
    // Select fitCat state
    fitCat->setLabel(fcState->GetName()) ;

    // Check if this fitCat state is selected
    Bool_t select(kTRUE) ;
    for (const auto arg : fitCatList) {
      auto splitCat = static_cast<const RooAbsCategory*>(arg);

      // Find selected state list 
      TList* slist = (TList*) splitStateList.FindObject(splitCat->GetName()) ;
      if (!slist) continue ;
      RooCatType* type = (RooCatType*) slist->FindObject(splitCat->getCurrentLabel()) ;
      if (!type) {
        select = kFALSE ;
      }
    }
    if (!select) continue ;


    // Select appropriate PDF for this physCat state
    RooCustomizer* physCustomizer ;
    if (physCat) {      
      RooStringVar* physNameVar = (RooStringVar*) stateMap.find(physCat->getCurrentLabel()) ;
      if (!physNameVar) continue ;
      physCustomizer = (RooCustomizer*) customizerList->FindObject(physNameVar->getVal());  
    } else {
      physCustomizer = (RooCustomizer*) customizerList->First() ;
    }

    coutI(ObjectHandling) << "RooSimPdfBuilder::buildPdf: Customizing physics model " << physCustomizer->GetName() 
       << " for mode " << fcState->GetName() << endl ;

    // Customizer PDF for current state and add to master simPdf
    RooAbsPdf* fcPdf = (RooAbsPdf*) physCustomizer->build(masterSplitCat.getCurrentLabel(),verbose) ;
    simPdf->addPdf(*fcPdf,fcState->GetName()) ;
  }
  delete fcIter ;

  // Move customizers (owning the cloned branch node components) to the attic
  _retiredCustomizerList.AddAll(customizerList) ;
  delete customizerList ;

  splitStateList.Delete() ;

  if (auxSplitCloneSet) delete auxSplitCloneSet ;

  _simPdfList.push_back(simPdf) ;
  _fitCatList.push_back(fitCat) ;
  return simPdf ;
}





////////////////////////////////////////////////////////////////////////////////

RooSimPdfBuilder::~RooSimPdfBuilder() 
{
  _retiredCustomizerList.Delete() ;

  std::list<RooSimultaneous*>::iterator iter = _simPdfList.begin() ;
  while(iter != _simPdfList.end()) {
    delete *iter ;
    ++iter ;
  }

  std::list<RooSuperCategory*>::iterator iter2 = _fitCatList.begin() ;
  while(iter2 != _fitCatList.end()) {
    delete *iter2 ;
    ++iter2 ;
  }
  
}
 

