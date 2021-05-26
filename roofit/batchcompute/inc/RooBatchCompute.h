#ifndef ROOFIT_BATCHCOMPUTE_ROOBATCHCOMPUTE_H
#define ROOFIT_BATCHCOMPUTE_ROOBATCHCOMPUTE_H

#include "rbc.h"
#include "Batches.h"
#include "RunContext.h"
#include "RooVDTHeaders.h"

#include "DllImport.h" //for R__EXTERN, needed for windows

#include <functional>

class RooAbsReal;
class RooListProxy;

/**
 * Namespace for dispatching RooFit computations to various backends.
 *
 * This namespace contains an interface for providing high-performance computation functions for use in RooAbsReal::evaluateSpan(),
 * see RooBatchComputeInterface.
 *
 * Furthermore, several implementations of this interface can be created, which reside in RooBatchCompute::RF_ARCH, where
 * RF_ARCH may be replaced by the architecture that this implementation targets, e.g. SSE, AVX, etc.
 *
 * Using the pointer RooBatchCompute::dispatch, a computation request can be dispatched to the fastest backend that is available
 * on a specific platform.
 */
namespace rbc = RooBatchCompute;
namespace RooBatchCompute {

enum Computer{AddPdf, ArgusBG, Bernstein, BifurGauss, BreitWigner, Bukin, CBShape, Chebychev, ChiSquare, 
              DstD0BG, Exponential, Gaussian, Gamma, Johnson, Lognormal, NegativeLogarithms, Novosibirsk, 
              Polynomial, ProdPdf};

struct RunContext;
/**
 * \brief The interface which should be implemented to provide optimised computation functions for implementations of RooAbsReal::evaluateSpan().
 *
 * This interface contains the signatures of the compute functions of every PDF that has an optimised implementation available.
 * These are the functions that perform the actual computations in batches.
 *
 * Several implementations of this interface may be provided, e.g. SSE, AVX, AVX2 etc. At run time, the fastest implementation of this interface
 * is selected, and using a virtual call, the computation is dispatched to the best backend.
 *
 * \see RooBatchCompute::dispatch, RooBatchComputeClass, RF_ARCH
 */ 
class RooBatchComputeInterface {
  public:
    virtual ~RooBatchComputeInterface() = default;
    virtual void   init() { throw std::bad_function_call(); }
    virtual void   compute(Computer, RestrictArr, size_t, const DataMap&, const VarVector&, const ArgVector& ={}) = 0;
    virtual double sumReduce(InputArr, size_t) = 0;
    virtual void*  malloc(size_t) = 0;
    virtual void   free(void*) = 0;
    virtual void   memcpyToGPU(void* dest, const void* src, size_t) { (void)dest; (void)src; throw std::bad_function_call(); }
    virtual void   memcpyToCPU(void* dest, const void* src, size_t) { (void)dest; (void)src; throw std::bad_function_call(); }

                            
    virtual RooSpan<double> computeLandau(const RooAbsReal*, RunContext&, RooSpan<const double> , RooSpan<const double> , RooSpan<const double> ) {return{};}
    virtual RooSpan<double> computeLognormal(const RooAbsReal*, RunContext&, RooSpan<const double> , RooSpan<const double> , RooSpan<const double> ) {return{};}
    virtual RooSpan<double> computeNovosibirsk(const RooAbsReal*, RunContext&, RooSpan<const double> , RooSpan<const double> , RooSpan<const double> , RooSpan<const double> ) {return{};}
    virtual RooSpan<double> computePoisson(const RooAbsReal*, RunContext&, RooSpan<const double> , RooSpan<const double> , bool , bool ) {return{};}
    //~  virtual void computePolynomial(size_t , double * __restrict , const double * __restrict const , int , std::vector<BracketAdapterWithMask>&) {return;}
    virtual RooSpan<double> computeVoigtian(const RooAbsReal*, RunContext&, RooSpan<const double> , RooSpan<const double> , RooSpan<const double> , RooSpan<const double> ) {return{};}    
};

/**
 * This dispatch pointer points to an implementation of the compute library, provided one has been loaded.
 * Using a virtual call, computation requests are dispatched to backends with architecture-specific functions
 * such as SSE, AVX, AVX2, etc.
 *
 * \see RooBatchComputeInterface, RooBatchComputeClass, RF_ARCH
 */
R__EXTERN RooBatchComputeInterface *dispatch, *dispatch_cpu, *dispatch_gpu;
} // End namespace RooBatchCompute

#endif
