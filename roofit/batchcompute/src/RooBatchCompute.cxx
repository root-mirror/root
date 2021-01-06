// RooBatchCompute library created September 2020 by Emmanouil Michalainas

#include "RooBatchCompute.h"
#include "Batches.h"

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
      Batches batches(output, nEvents, varData, vars, extraArgs, buffer);
      batches.setNEvents(bufferSize);
      while (nEvents > bufferSize)
      {
        computeFunctions[computer](batches);
        batches.advance();
        nEvents -= bufferSize;
      }
      batches.setNEvents(nEvents);
      computeFunctions[computer](batches);
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
