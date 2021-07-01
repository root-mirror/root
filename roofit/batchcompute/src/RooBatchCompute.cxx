// RooBatchCompute library created September 2020 by Emmanouil Michalainas

#include "RooBatchCompute.h"
#include "Batches.h"

#include "ROOT/TExecutor.hxx"

#include <cstdlib>

namespace RooBatchCompute {
namespace RF_ARCH {

std::vector<void(*)(Batches)> getFunctions();

class RooBatchComputeClass : public RooBatchComputeInterface {
  private:
    const std::vector<void(*)(Batches)> computeFunctions;
  public:
    RooBatchComputeClass()
      : computeFunctions(getFunctions())
    {
      // Set the dispatch pointer to this instance of the library upon loading
      dispatch_cpu = this;
    }

    void compute(Computer computer, RestrictArr output, size_t nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& extraArgs) override  
    {
      double buffer[maxParams][bufferSize];
      ROOT::Internal::TExecutor ex;
      unsigned int nThreads = ROOT::IsImplicitMTEnabled() ? ex.GetPoolSize() : 1u;

      // Fill a std::vector<Batches> with the same object and with ~nEvents/nThreads
      // Then advance every object but the first to split the work between threads
      std::vector<Batches> batches(nThreads, Batches(output, nEvents/nThreads +(nEvents%nThreads>0), varData, vars, extraArgs, buffer));
      for (unsigned int i=1; i<nThreads; i++)
        batches[i].advance( batches[0].getNEvents()*i );

      // Set the number of events of the last Btches object as the remaining events
      if (nThreads>1)
        batches.back().setNEvents( nEvents -(nThreads-1)*batches[0].getNEvents() );

      auto task = [this, computer](Batches _batches) -> int
      {
          int _events = _batches.getNEvents();
          _batches.setNEvents(bufferSize);
          while (_events > bufferSize)
          {
            computeFunctions[computer](_batches);
            _batches.advance(bufferSize);
            _events -= bufferSize;
          }
          _batches.setNEvents(_events);
          computeFunctions[computer](_batches);
          return 0;
      };
      ex.Map(task, batches);
    }

    double sumReduce(InputArr input, size_t n) override
    {
      long double sum=0.0;
      for (size_t i=0; i<n; i++) sum += input[i];
      return sum;
    }

    void* malloc(size_t size) override
    {
      return std::malloc(size);
    }

    void free(void* ptr) override
    {
      std::free(ptr);
    }
}; // End class RooBatchComputeClass

/// Static object to trigger the constructor which overwrites the dispatch pointer.
static RooBatchComputeClass computeObj;

Batches::Batches(RestrictArr _output, size_t _nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& _extraArgs, double stackArr[maxParams][bufferSize])
  : nEvents(_nEvents), nBatches(vars.size()), nExtraArgs(_extraArgs.size()), output(_output)
{
  for (size_t i=0; i<vars.size(); i++)
  {
    const RooSpan<const double>& span = varData.at(vars[i]);
    if (span.size()>1) arrays[i].set(span.data()[0], span.data(), true);
    else
    {
      std::fill_n(stackArr[i], bufferSize, span.data()[0]);
      arrays[i].set(span.data()[0], stackArr[i], false);
    }
  }
  std::copy(_extraArgs.cbegin(), _extraArgs.cend(), extraArgs);
}

} // End namespace RF_ARCH
} //End namespace RooBatchCompute
