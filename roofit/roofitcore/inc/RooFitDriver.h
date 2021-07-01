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

class RooFitDriver {
  public:
     RooFitDriver(const RooAbsData& data, const RooNLLVarNew& topNode, int batchMode);
     ~RooFitDriver();
     double getVal();
     RooArgSet* getParameters() const;
     
    struct NodeInfo {
        int nServers = 0;
        int nClients = 0;
        bool dependsOnObservables = true;
        bool computeInScalarMode = false;
    };

  private:
    double* getAvailableBuffer();

    const int _batchMode=0;
    double* _cudaMemDataset=nullptr;

    // used for preserving static info about the computation graph
    rbc::DataMap _dataMap;
    const RooNLLVarNew& _topNode;
    size_t _nEvents;
    RooAbsData const* _data = nullptr;
    std::queue<const RooAbsReal*> _initialQueue;
    std::unordered_map<const RooAbsArg*,NodeInfo> _nodeInfos;

    // used for dynamically scheduling each step's computations
    std::queue<const RooAbsReal*> _computeQueue;
    std::queue<double*> _vectorBuffers;
};

#endif //ROO_FIT_DRIVER_H
