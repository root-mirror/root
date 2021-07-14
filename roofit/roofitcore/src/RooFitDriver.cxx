#include "RooFitDriver.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooBatchCompute.h"
#include "RooNLLVarNew.h"
#include "RooRealVar.h"
#include "RunContext.h"

RooFitDriver::RooFitDriver(const RooAbsData& data, const RooNLLVarNew& topNode, int batchMode)
  : _batchMode{batchMode}, _topNode{topNode}, _data{&data},
    _nEvents{static_cast<size_t>( data.numEntries() )}
{
  // fill the RunContext with the observable data and map the observables
  // by namePtr in order to replace their memory addresses later, with
  // the ones from the variables that are actually in the computation graph. 
  rbc::RunContext evalData;
  data.getBatches(evalData, 0, _nEvents);
  _dataMapCPU = evalData.spans;
  std::unordered_map<const TNamed*,const RooAbsReal*> nameResolver;
  for (auto& it:_dataMapCPU) nameResolver[it.first->namePtr()]=it.first;

  // Check if there is a batch for weights and if it's already in the dataMap.
  // If not, we need to put the batch and give as a key a RooRealVar* that has
  // the same name as RooNLLVarNew's _weight proxy, so that it gets renamed like
  // every other observable.
  RooSpan<const double> weights = data.getWeightBatch(0, _nEvents);
  std::string weightVarName = data.getWeightVarName()!="" ? data.getWeightVarName() : "_weight";
  RooRealVar dummy (weightVarName.c_str(), "dummy", 0.0);
  const TNamed* pTNamed = dummy.namePtr();
  if (!weights.empty() && nameResolver.count(pTNamed)==0)
  {
    _dataMapCPU[&dummy] = weights;
    nameResolver[pTNamed] = &dummy;
  }

  // Get a serial list of the nodes in the computation graph.
  // treeNodeServelList() is recursive and adds the top node before the children,
  // so reversing the list gives us a topological ordering of the graph.
  RooArgList list;
  _topNode.treeNodeServerList(&list);
  for (int i=list.size()-1; i>=0; i--)
  {
    auto pAbsReal = dynamic_cast<RooAbsReal*>(&list[i]);
    if (!pAbsReal) continue;
    const bool alreadyExists = nameResolver.count(pAbsReal->namePtr());
    const RooAbsReal* pClone = nameResolver[pAbsReal->namePtr()];
    if (alreadyExists && !pClone) continue; // node included multiple times in the list
      
    if (pClone) //this node is an observable, update the RunContext and don't add it in `_initialQueue`.
    {
      auto it = _dataMapCPU.find(pClone);
      _dataMapCPU[pAbsReal]=it->second;
      _dataMapCPU.erase(it);

      // set nameResolver to nullptr to be able to detect future duplicates
      nameResolver[pAbsReal->namePtr()] = nullptr;
    }
    else //this node needs evaluation, mark it's clients
    {
      // If the node doesn't depend on any observables, there is no need to
      // loop over events and we don't need to use the batched evaluation.
      RooArgSet observablesForNode;
      pAbsReal->getObservables(_data->get(), observablesForNode);
      _nodeInfos[pAbsReal].computeInScalarMode = observablesForNode.empty() || !pAbsReal->isDerived();

      if (_nodeInfos[pAbsReal].nServers==0) _initialQueue.push(pAbsReal);

      auto clients = pAbsReal->valueClients();
      for (auto* client:clients)
        if(list.find(*client))
        {
          auto pClient = static_cast<const RooAbsReal*>(client);
          ++_nodeInfos[pClient].nServers;
          ++_nodeInfos[pAbsReal].nClients;
          // client needs to wait for this node to get evaluated first, so it shouldn't
          // be in the initial queue
          _initialQueue.erase(pClient);

          // If this node's client is not computed in the same place (CPU/GPU)
          // then we need to make a copy after it gets evaluated
          if (_batchMode==-1 && pAbsReal->canComputeBatchWithCuda() != pClient->canComputeBatchWithCuda())
            _nodeInfos[pAbsReal].copyAfterEvaluation = true;
        }
    }
  }
  rbc::dispatch = rbc::dispatch_cpu;
  // If cuda mode is on, copy all observable data to device memory
  if (_batchMode == -1)
  {
    if (!rbc::dispatch_gpu) 
      throw std::runtime_error(std::string("In: ")+__func__+"(), "+__FILE__+":"+__LINE__+": Cuda implementation of the computing library is not available\n");
    rbc::dispatch = rbc::dispatch_gpu;
    _cudaMemDataset = static_cast<double*>(rbc::dispatch->malloc( _nEvents*_dataMapCPU.size()*sizeof(double) ));
    size_t idx=0;
    for (auto& record:_dataMapCPU)
    {
      _dataMapGPU[record.first] = RooSpan<double>(_cudaMemDataset+idx, _nEvents);
      rbc::dispatch->memcpyToGPU(_cudaMemDataset+idx, record.second.data(), _nEvents*sizeof(double));
      idx += _nEvents;
    }
  }
}

RooFitDriver::~RooFitDriver()
{
  !!!!!!!!!!!!!!!!
  while (!_vectorBuffers.empty())
  {
    rbc::dispatch->free( _vectorBuffers.front() );
    _vectorBuffers.pop();
  }
  rbc::dispatch->free(_cudaMemDataset);
}

double RooFitDriver::getVal()
{
  std::queue<const RooAbsReal*> cpuQueue;
  std::list<cudaStream_t*> activeStreams;
  std::unordered_map<cudaStream_t*, const RooAbsReal*> streamToNode;
  std::unordered_map<const RooAbsReal*, NodeInfo> remaining = _nodeInfos;
  std::vector<double> nonDerivedValues;
  nonDerivedValues.reserve(_nodeInfos.size()); // to avoid reallocation

  for (auto pAbsReal:_initialQueue)
    if (!_nodeInfos.at(pAbsReal).computeInGPU) cpuQueue.push(pAbsReal);
    else assignToGPU(pAbsReal, remaining, activeStreams, streamToNode);

  while (!remaining.empty())
  {
    // STEP 1: Compute next cpu node (if exists)
    if (!cpuQueue.empty())
    {
      auto node = cpuQueue.front();
      cpuQueue.pop();
      if (remaining.at(node).computeInScalarMode)
      {
        nonDerivedValues.push_back(node->getVal(_data->get()));
        _dataMapCPU[node] = _dataMapGPU[node] = RooSpan<const double>(&nonDerivedValues.back(),1);
      }
      else
      {
        !!!!!!!!!!!!!!!!!!!!
        bool copyAfterEvaluation=false;
        for (auto* client : node->valueClients())
          if (static_cast<const RooAbsReal*>(client)->canComputeBatchWithCUDA())
            copyAfterEvaluation=true;
        
        // get an available buffer for storing the comptation results
        double* buffer = copyAfterEvaluation ? getAvailablePinnedBuffer() : getAvailableCPUBuffer();

        // TODO: Put integrals seperately in the computation queue
        // For now, we just assume they are scalar and assign them some temporary memory
        double normVal=1.0;
        if (auto pAbsPdf = dynamic_cast<const RooAbsPdf*>(node))
        {
          auto integral = pAbsPdf->getIntegral(*_data->get());
          normVal = integral->getVal();
          _dataMapCPU[integral] = _dataMapGPU[integral] = RooSpan<const double>(&normVal,1);
        }

        // compute this node and register the result in the dataMap
        node->computeBatch(buffer, _nEvents, _dataMapCPU);
        _dataMapCPU[node] = RooSpan<const double>(buffer, _nEvents);
        if (copyAfterEvaluation)
        {
          cudaStream_t* stream = getAvailableCudaStream();
          remaining.at(node).computeStage = 2;
          remaining.at(node).stream = stream;
          streamToNode[stream]=node;
          double* gpuBuffer = getAvailableGPUBuffer();
          _dataMapGPU[node] = RooSpan<const double>(gpuBuffer, _nEvents);
          rbc::dispatch->memcpyToGPU(gpuBuffer, buffer, _nEvents*sizeof(double), stream);
        }
        else remaining.at(node0).computeStage = 3;
      } // else compute in batch mode

        // STEP 2: handle the node that was just computed by the CPU
        updateMyServers(node, remaining);
        checkMyClients(node, remaining);
    } // if (!cpuQueue.empty())

    // STEP 3: handle nodes computed by the GPU while the CPU was computing a node
    if (_batchMode == -1)
      for (auto it=activeStreams.begin(); it!=activeStreams.end(); )
        if (!rbc::dispatch->cudaStreamHasFinished(*it)) it++;
        else
        {
          remaining.at(streamToNode.at(*it)).stream = nullptr;
          remaining.at(streamToNode.at(*it)).computeStage = 3;
          updateMyServers(streamToNode.at(*it), remaining);
          checkMyClients(streamToNode.at(*it), remaining);
          _cudaStreamBuffers.push(*it);
          it = activeStreams.erase(it);
        }
  } // while (!remaining.empty())
  
  // recycle the top node's buffer and return the final value
  _vectorBuffers.push(const_cast<double*>( _dataMapCPU[&_topNode].data() ));
  return _topNode.reduce(_dataMapCPU[&_topNode].data(), _nEvents);
}

void RooFitDriver::updateMyServers(const RooAbsReal* node, std::unordered_map<const RooAbsReal*, NodeInfo>& nodeInfos)
{
  for (auto* server : node->servers())
  {
    auto pServer = static_cast<const RooAbsReal*>(server);
    auto& info = nodeInfos[pServer];
    if (--info.nClients>0) continue;
    if (!info.computeInScalarMode)
      if (info.copyAfterEvaluation)
      {
        _gpuBuffers.push(const_cast<double*>( _dataMapGPU[pServer].data() );
        _pinnnedBuffers.push(const_cast<double*>( _dataMapCPU[pServer].data() ));
      }
      else if (info.computeInGPU)
        _gpuBuffers.push(const_cast<double*>( _dataMapGPU[pServer].data() );
      else
        _cpuBuffers.push(const_cast<double*>( _dataMapCPU[pServer].data() ));
    nodeInfos.erase(pServer);
  }
}

void checkMyClients(const RooAbsReal* node, std::unordered_map<const RooAbsReal*, NodeInfo>& nodeInfos)
{
  for (auto* client : node->valueClients())
  {
    auto pClient = static_cast<const RooAbsReal*>(client); 
    if (nodeInfos.count(pClient)==0) continue;
    if (nodeInfos.at(pClient).computeStage>0) continue;
    for (auto* depend : pClient->servers())
    {
      pDependency = static_cast<const RooAbsReal*>(depend);
      if (nodeInfos.at(pDependency).computeStage==0
      computeQueue.insert(pClient);
    }
  }
}


double* RooFitDriver::getAvailableBuffer()
{
  if (_vectorBuffers.empty())
    return static_cast<double*>(rbc::dispatch->malloc( _nEvents*sizeof(double) ));
  else
  {
    double* buffer = _vectorBuffers.front();
    _vectorBuffers.pop();
    return buffer;
  }
}

RooArgSet* RooFitDriver::getParameters() const
{
  return _topNode.getParameters(_data->get(), true);
}
