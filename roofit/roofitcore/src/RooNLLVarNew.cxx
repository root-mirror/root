#include "RooNLLVarNew.h"
#include "RooBatchCompute.h"

#include <numeric>
#include <stdexcept>
#include <vector>

ClassImp(RooNLLVarNew);

RooNLLVarNew::RooNLLVarNew(const char *name, const char *title,
                           RooAbsPdf &pdf, RooArgSet const& observables,
                           RooAbsReal* weight, RooAbsReal* constraints, bool isExtended)
   : RooAbsReal(name, title),
    _pdf{"pdf", "pdf", this, pdf},
    _observables{&observables},
    _isExtended{isExtended}
{
   if(weight) {
       _weight = std::make_unique<RooTemplateProxy<RooAbsReal>>("_weight", "_weight", this, *weight);
    }
   if(constraints) {
       _constraints = constraints;
   }
}

RooNLLVarNew::RooNLLVarNew(const RooNLLVarNew &other, const char *name)
   : RooAbsReal(other, name),
    _pdf{"pdf", this, other._pdf},
    _observables{other._observables}
{
    if(other._weight)
        _weight = std::make_unique<RooTemplateProxy<RooAbsReal>>("_weight", this, *other._weight);
    if(other._constraints)
        _constraints = other._constraints;
}

void RooNLLVarNew::computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const 
{
  rbc::VarVector vars = {&*_pdf};
  if (_weight) {
      vars.push_back(&**_weight);
  }
  rbc::ArgVector args = {static_cast<double>(vars.size()-1)};
  rbc::dispatch->compute(rbc::NegativeLogarithms, output, nEvents, dataMap, vars, args);

  if (_isExtended && _sumWeight == 0.0) {
    if(!_weight) {
      _sumWeight = nEvents;
    } else {
      auto weightSpan = dataMap[&**_weight];
      _sumWeight = weightSpan.size() == 1 ? weightSpan[0] * nEvents
                                          : rbc::dispatch->sumReduce(dataMap[&**_weight].data(), nEvents);
    }
  }
}

double RooNLLVarNew::reduce(const double* input, size_t nEvents) const
{
  double nll = rbc::dispatch->sumReduce(input, nEvents);
  if (_constraints) {
    nll += _constraints->getVal();
  }
  if (_isExtended) {
    assert(_sumWeight != 0.0);
    nll += _pdf->extendedTerm(_sumWeight, _observables);
  }
  return nll;
}

double RooNLLVarNew::getValV(const RooArgSet *) const
{
   throw std::runtime_error("RooNLLVarNew::getValV was called directly which should not happen!");
}

double RooNLLVarNew::evaluate() const
{
   throw std::runtime_error("RooNLLVarNew::evaluate was called directly which should not happen!");
}

RooSpan<double> RooNLLVarNew::evaluateSpan(RooBatchCompute::RunContext &, const RooArgSet *) const
{
   throw std::runtime_error("RooNLLVarNew::evaluatSpan was called directly which should not happen!");
}

RooSpan<const double> RooNLLVarNew::getValues(RooBatchCompute::RunContext &, const RooArgSet *) const
{
   throw std::runtime_error("RooNLLVarNew::getValues was called directly which should not happen!");
}

bool RooNLLVarNew::getParameters(const RooArgSet *depList, RooArgSet& outSet, bool stripDisconnected) const
{
   return _pdf->getParameters(depList, outSet, stripDisconnected);
}
