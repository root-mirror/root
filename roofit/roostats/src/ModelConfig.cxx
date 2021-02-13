// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Sven Kreiss
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::ModelConfig
    \ingroup Roostats

ModelConfig is a simple class that holds configuration information specifying how a model
should be used in the context of various RooStats tools.  A single model can be used
in different ways, and this class should carry all that is needed to specify how it should be used.
ModelConfig requires a workspace to be set.

A ModelConfig holds sets of parameters of the likelihood function that have different interpretations:
- **Parameter of interest** Parameters that are measured (*i.e.* fitted).
- **Nuisance parameters** Parameters that are fitted, but their post-fit value is not interesting. Often,
they might be constrained because external knowledge about them exists, *e.g.* from external measurements.
- **Constraint parameters** No direct use in RooFit/RooStats. Can be used by the user for bookkeeping.
- **Observables** Parameters that have been measured externally, *i.e.* they exist in a dataset. These are not fitted,
but read during fitting from the entries of a dataset.
- **Conditional observables** Observables that are not integrated when the normalisation of the PDF is calculated.
See *e.g.* `rf306_condpereventerrors` in the RooFit tutorials.
- **Global observables** Observables that to the fit look like "constant" values, *i.e.* they are not being
fitted and they are not loaded from a dataset, but some knowledge exists that allows to set them to a
specific value. Examples:
-- A signal efficiency measured in a Monte Carlo study.
-- When constraining a parameter \f$ b \f$, the target value (\f$ b_0 \f$) that this parameter is constrained to:
\f[
  \mathrm{Constraint}_b = \mathrm{Gauss}(b_0 \, | \, b, 0.2)
\f]
*/

#include "RooStats/ModelConfig.h"

#include "RooMsgService.h"

#include "RooStats/RooStatsUtils.h"

#include <sstream>


ClassImp(RooStats::ModelConfig);

using namespace std;

namespace RooStats {

////////////////////////////////////////////////////////////////////////////////
/// Makes sensible guesses of observables, parameters of interest
/// and nuisance parameters if one or multiple have been set by the creator of this ModelConfig.
///
/// Defaults:
/// - Observables: determined from data,
/// - Global observables: explicit obs  -  obs from data  -  constant observables
/// - Parameters of interest: empty,
/// - Nuisance parameters: all parameters except parameters of interest
///
/// We use NULL to mean not set, so we don't want to fill
/// with empty RooArgSets.

void ModelConfig::GuessObsAndNuisance(const RooAbsData& data, bool printModelConfig) {

   // observables
  if (!GetObservables()) {
     const RooArgSet * obs = GetPdf()->getObservables(data);
     SetObservables(*obs);
     delete obs;
   }
  // global observables
   if (!GetGlobalObservables()) {
      RooArgSet co(*GetObservables());
      const RooArgSet * obs = GetPdf()->getObservables(data);
      co.remove(*obs);
      RemoveConstantParameters(&co);
      if(co.getSize()>0)
   SetGlobalObservables(co);

      // TODO BUG This does not work as observables with the same name are already in the workspace.
      /*
      RooArgSet o(*GetObservables());
      o.remove(co);
      SetObservables(o);
      */
      delete obs;
   }

   // parameters
   //   if (!GetParametersOfInterest()) {
   //      SetParametersOfInterest(RooArgSet());
   //   }
   if (!GetNuisanceParameters()) {
      const RooArgSet * params = GetPdf()->getParameters(data);
      RooArgSet p(*params);
      p.remove(*GetParametersOfInterest());
      RemoveConstantParameters(&p);
      if(p.getSize()>0)
   SetNuisanceParameters(p);
      delete params;
   }

   // print Modelconfig as an info message

   if (printModelConfig) {
     std::ostream& oldstream = RooPrintable::defaultPrintStream(&ccoutI(InputArguments));
     Print();
     RooPrintable::defaultPrintStream(&oldstream);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print contents of Model on the default print stream
/// It can be changed using RooPrintable

void ModelConfig::Print(Option_t*) const {
   ostream& os = RooPrintable::defaultPrintStream();

   os << endl << "=== Using the following for " << GetName() << " ===" << endl;


   // args
   if(GetObservables()){
      os << "Observables:             ";
      GetObservables()->Print("");
   }
   if(GetParametersOfInterest()) {
      os << "Parameters of Interest:  ";
      GetParametersOfInterest()->Print("");
   }
   if(GetNuisanceParameters()){
      os << "Nuisance Parameters:     ";
      GetNuisanceParameters()->Print("");
   }
   if(GetGlobalObservables()){
      os << "Global Observables:      ";
      GetGlobalObservables()->Print("");
   }
   if(GetConstraintParameters()){
      os << "Constraint Parameters:   ";
      GetConstraintParameters()->Print("");
   }
   if(GetConditionalObservables()){
      os << "Conditional Observables: ";
      GetConditionalObservables()->Print("");
   }
   if(GetProtoData()){
      os << "Proto Data:              ";
      GetProtoData()->Print("");
   }

   // pdfs
   if(GetPdf()) {
      os << "PDF:                     ";
      GetPdf()->Print("");
   }
   if(GetPriorPdf()) {
      os << "Prior PDF:               ";
      GetPriorPdf()->Print("");
   }

   // snapshot
   const RooArgSet * snapshot = GetSnapshot();
   if(snapshot) {
      os << "Snapshot:                " << endl;
      snapshot->Print("v");
      delete snapshot;
   }

   os << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// If a workspace already exists in this ModelConfig, RooWorkspace::merge(ws) will be called
/// on the existing workspace.

void ModelConfig::SetWS(RooWorkspace & ws) {
   if( !fRefWS.GetObject() ) {
      fRefWS = &ws;
      fWSName = ws.GetName();
   }
   else{
      RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
      GetWS()->merge(ws);
      RooMsgService::instance().setGlobalKillBelow(level) ;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// get from TRef

RooWorkspace * ModelConfig::GetWS() const {
   RooWorkspace *ws = dynamic_cast<RooWorkspace *>(fRefWS.GetObject() );
   if(!ws) {
      coutE(ObjectHandling) << "workspace not set" << endl;
      return NULL;
   }
   return ws;
}

////////////////////////////////////////////////////////////////////////////////
/// save snapshot in the workspace
/// and use values passed with the set

void ModelConfig::SetSnapshot(const RooArgSet& set) {
   if ( !GetWS() ) return;

   fSnapshotName = GetName();
   if (fSnapshotName.size()  > 0) fSnapshotName += "_";
   fSnapshotName += set.GetName();
   if (fSnapshotName.size()  > 0) fSnapshotName += "_";
   fSnapshotName += "snapshot";
   GetWS()->saveSnapshot(fSnapshotName.c_str(), set, true);  // import also the given parameter values
   DefineSetInWS(fSnapshotName.c_str(), set);
}

////////////////////////////////////////////////////////////////////////////////
/// Load the snapshot from ws and return the corresponding set with the snapshot values.
/// User must delete returned RooArgSet.

const RooArgSet * ModelConfig::GetSnapshot() const{
   if ( !GetWS() ) return 0;
   if (!fSnapshotName.length()) return 0;
   // calling loadSnapshot will also copy the current parameter values in the workspaces
   // since we do not want to change the model parameters - we restore the previous ones
   if (! GetWS()->set(fSnapshotName.c_str() ) )return 0;
   RooArgSet snapshotVars(*GetWS()->set(fSnapshotName.c_str() ) );
   if (snapshotVars.getSize() == 0) return 0;
   // make my snapshot which will contain a copy of the snapshot variables
   RooArgSet tempSnapshot;
   snapshotVars.snapshot(tempSnapshot);
   // load snapshot value from the workspace
   if (!(GetWS()->loadSnapshot(fSnapshotName.c_str())) ) return 0;
   // by doing this snapshotVars will have the snapshot values - make the snapshot to return
   const RooArgSet * modelSnapshot = dynamic_cast<const RooArgSet*>( snapshotVars.snapshot());
   // restore now the variables of snapshot in ws to their original values
   // need to const cast since assign is not const (but in reality in just assign values and does not change the set)
   // and anyway the set is const
   snapshotVars.assignFast(tempSnapshot);
   return modelSnapshot;
}

////////////////////////////////////////////////////////////////////////////////
/// load the snapshot from ws if it exists

void ModelConfig::LoadSnapshot() const{
   if ( !GetWS() ) return;
   GetWS()->loadSnapshot(fSnapshotName.c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// helper functions to avoid code duplication

void ModelConfig::DefineSetInWS(const char* name, const RooArgSet& set) {
   if ( !GetWS() ) return;

   const RooArgSet * prevSet = GetWS()->set(name);
   if (  prevSet ) {
      //be careful not to remove passed set in case it is the same updated
      if (prevSet != &set)
         GetWS()->removeSet(name);
   }

   // suppress warning when we re-define a previously defined set (when set == prevSet )
   // and set is not removed in that case
   RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;


   GetWS()->defineSet(name, set,true);

   RooMsgService::instance().setGlobalKillBelow(level) ;

}

////////////////////////////////////////////////////////////////////////////////
/// internal function to import Pdf in WS

void ModelConfig::ImportPdfInWS(const RooAbsPdf & pdf) {
   if ( !GetWS() ) return;

   if (! GetWS()->pdf( pdf.GetName() ) ){
      RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
      GetWS()->import(pdf, RooFit::RecycleConflictNodes());
      RooMsgService::instance().setGlobalKillBelow(level) ;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// internal function to import data in WS

void ModelConfig::ImportDataInWS(RooAbsData & data) {
   if ( !GetWS() ) return;

   if (! GetWS()->data( data.GetName() ) ){
      RooFit::MsgLevel level = RooMsgService::instance().globalKillBelow();
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
      GetWS()->import(data);
      RooMsgService::instance().setGlobalKillBelow(level) ;
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t ModelConfig::SetHasOnlyParameters(const RooArgSet& set, const char* errorMsgPrefix) {

    RooArgSet nonparams ;
    RooFIter iter = set.fwdIterator() ;
    RooAbsArg* arg ;
    while ((arg=iter.next())) {
      if (!arg->isFundamental()) {
         nonparams.add(*arg) ;
      }
    }

    if (errorMsgPrefix && nonparams.getSize()>0) {
      cout << errorMsgPrefix << " ERROR: specified set contains non-parameters: " << nonparams << endl ;
    }
    return (nonparams.getSize()==0) ;
  }

} // end namespace RooStats
