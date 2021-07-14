#ifndef ROO_FIT_DRIVER_H
#define ROO_FIT_DRIVER_H

#include "rbc.h"

#include <queue>
#include <unordered_map>

class RooAbsData;
class RooAbsArg;
class RooAbsReal;
class RooArgSet;
class RooNLLVarNew;
struct cudaStream_t;
struct cudaEvent_t;

class RooFitDriver {
  public:
     RooFitDriver(const RooAbsData& data, const RooNLLVarNew& topNode, int batchMode);
     ~RooFitDriver();
     double getVal();
     RooArgSet* getParameters() const;
     
    struct NodeInfo {
      int nClients = 0;
      int computeStage = 0; // 0=not_processed, 1=computing, 2=copying, 3=finished
      bool computeInScalarMode = false;
      bool computeInGPU = false;
      bool copyAfterEvaluation = false;
      cudaStream_t* stream = nullptr;
      cudaEvent_t* eventAfterComputation = nullptr;
    };

  private:
    double* getAvailableCPUBuffer();
    double* getAvailableGPUBuffer();
    double* getAvailablePinnedBuffer();
    cudaStream_t* getAvailableCudaStream();
    void updateMyServers(const RooAbsReal* node, std::unordered_map<const RooAbsReal*, NodeInfo>&);
    void checkMyClients(const RooAbsReal* node, std::unordered_map<const RooAbsReal*, NodeInfo>&);

    const int _batchMode = 0;
    double* _cudaMemDataset = nullptr;

    // used for preserving static info about the computation graph
    rbc::DataMap _dataMapCPU;
    rbc::DataMap _dataMapGPU;
    const RooNLLVarNew& _topNode;
    const RooAbsData* const _data = nullptr;
    const size_t _nEvents;
    std::vector<const RooAbsReal*> _initialQueue;
    std::unordered_map<const RooAbsReal*, NodeInfo> _nodeInfos;

    // used for dynamically scheduling each step's computations
    std::queue<double*> _cpuBuffers;
    std::queue<double*> _gpuBuffers;
    std::queue<double*> _pinnedBuffers;
    std::queue<cudaStream_t*> _cudaStreamBuffers;
};

#endif //ROO_FIT_DRIVER_H
