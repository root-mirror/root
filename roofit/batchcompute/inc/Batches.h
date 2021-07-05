// Author: Emmanouil Michalainas, CERN 3 March 2021

#ifndef ROOFIT_BATCHCOMPUTE_BATCHES_H
#define ROOFIT_BATCHCOMPUTE_BATCHES_H

#include "rbc.h"

#include <stdint.h>

namespace RooBatchCompute {

constexpr uint8_t maxParams=8;
constexpr uint8_t maxExtraArgs=16;
constexpr uint16_t bufferSize=64;

namespace RF_ARCH {

class Batch {
  private:
    double _scalar=0;
    const double* __restrict _array=nullptr;
    bool _isVector=false;

  public:
    Batch() = default;
    inline Batch(InputArr array, bool isVector)
      : _scalar{array[0]}, _array{array}, _isVector{isVector}
    {}
    
    __device__ constexpr bool isItVector() const { return _isVector; }
    inline void set(double scalar, InputArr array, bool isVector)
    {
      _scalar = scalar;
      _array = array;
      _isVector = isVector;
    }
    inline void advance(size_t _nEvents) { _array += _isVector*_nEvents; }
#ifdef __CUDACC__
    __device__ constexpr double operator[](size_t i) const noexcept { return _isVector ? _array[i] : _scalar; }
#else
    constexpr double operator[](size_t i) const noexcept { return _array[i]; }
#endif // #ifdef __CUDACC__
}; //end class Batch

/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Batches {
  private:
    Batch _arrays[maxParams];
    size_t _nEvents=0;
    double _extraArgs[maxExtraArgs];
    uint8_t _nBatches=0;
    uint8_t _nExtraArgs;

  public:
    RestrictArr _output=nullptr;
    
    Batches(RestrictArr output, size_t nEvents, const DataMap& varData, const VarVector& vars, const ArgVector& extraArgs={}, double stackArr[maxParams][bufferSize]=nullptr);
    __device__ constexpr size_t getNEvents() const { return _nEvents; }
    __device__ constexpr uint8_t getNExtraArgs() const { return _nExtraArgs; }
    __device__ constexpr double extraArg(uint8_t i) const { return _extraArgs[i]; }
    __device__ constexpr Batch operator[] (int batchIdx) const { return _arrays[batchIdx]; }
    inline void setNEvents(size_t n=bufferSize) { _nEvents = n; }
    inline void advance(size_t nEvents)
    {
      for (int i=0; i<_nBatches; i++) _arrays[i].advance(nEvents);
      _output += nEvents;
    }
}; //end class Batches
} // End namespace RF_ARCH
} //end namespace RooBatchCompute
#endif // #ifdef ROOFIT_BATCHCOMPUTE_BATCHES_H 
