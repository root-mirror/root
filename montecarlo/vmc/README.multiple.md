# The Virtual Monte Carlo (VMC) Package

## Overview

The VMC package provides an abstract interface for Monte Carlo transport engines. Interfaces are implemented for
* [GEANT3](https://github.com/vmc-project/geant3)
* [GEANT4](https://github.com/vmc-project/geant4_vmc)

each deriving from `TVirtualMC` implementing all necessary methods.

Before a user can instantiate an engine, an object deriving from `TVirtualMCApplication` needs to be present which has to be implemented by the user. It contains necessary hooks called from the `TVirtualMC`s depending on their internal state. At the same time it provides the bridge between the user code and VMC. For instance, the user code can contain the geometry construction routines somewhere which should be called from the implemented `UserApplication::ConstructGeometry()`.

Further general information on the VMC project can be found [here](https://root.cern/vmc)

## Running multiple different engines

The simulation of an event can be shared among multiple different engines deriving from `TVirtualMC` which are handled by a singleton `TMCManager` object. In such a scenario the user has to call `TVirtualMCApplication::RequestMCManager()` in the constructor of the user application. A pointer to the manager object is then available via the protected `TVirtualMCApplication::fMCManager` but can also be obtained using the static method `TMCManager *TMCManager::Instance()`.

`TMCManager` provides the following interfaces:

* `void SetUserStack(TVirtualMCStack* userStack)` notifies the manager on the user stack such that it will be kept up-to-date during the simulation. Running without having set the user stack is not possible and the `TMCManager` will abort in that case.
* `void ForwardTrack(Int_t toBeDone, Int_t trackId, Int_t parentId, TParticle* userParticle)`: The user is still the owner of all track objects (aka `TParticle`) being created. Hence, all engine calls to `TVirtualMCStack::PushTrack(...)` are forwarded to the user stack. This can then invoke the `ForwardTrack(..)` method of the manager to pass the pointer to the constructed `TParticle` object. If a particle should be pushed to an engine other than the one currently running the engine's id has to be provided as the last argument.
* `void TransferTrack(Int_t targetEngineId)`: E.g. during `TVirtualMCApplication::Stepping()` the user might decide that the current track should be transferred to another engine, for instance, if a certain volume is entered. Specifying the ID of the target engine the manager will take care of interrupting the track in the current engine, extracting the kinematics and geometry state and it will push this to the stack of the target engine.
* `template <typename F> void Apply(F f)` assumes `f` to implement the `()` operator and taking a `TVirtualMC` pointer as an arument. `f` will be then called for all engines.
* `template <typename F> void Init(F f)` works as `TMCManager::Apply` during the initialization of the engines. It can also be called without an argument such that no additional user routine is included.
* `void Run(Int_t nEvents)` steers a run for the specified number of events.
* `void ConnectEnginePointers(TVirtualMC *&mc)` gives the possibility for a user to pass a pointer which will always be set to point to the currently running engine.
* `TVirtualMC *GetCurrentEngine()` provides the user with the currently running engine.

An example of how the `TMCManager` is utilized in a multi-run can be found in `examples/EME` of the [GEANT4_VMC repository](https://github.com/vmc-project/geant4_vmc).

### Workflow

**Implementation**
1. Implement your application as you have done before. Request the `TMCManager` in your constructor if needed via `TVirtualMCApplication::RequestMCManager()`
2. Implement your user stack as you have done before. At an appropriate stage (e.g. in `UserStack::PushTrack(...)`) you should call `TMCManager::ForwardTrack(...)` to forward the pointers to your newly constructed `TParticle` objects.
3. Set your stack using `TMCManager::SetUserStack(...)`.

**Usage**
1. Instantiate your application
2. Instantiate the engines you want to use.
3. Call `TMCManager::Init(...)`.
4. Call `TMCManager::Run(...)`

**Further comments**

The geometry is built once centrally via the `TMCManager` calling

1. `TVirtualMCApplication::ConstructGeometry()`
2. `TVirtualMCApplication::MisalignGeometry()`
2. `TVirtualMCApplication::ConstructOpGeometry()`

so it is expected that these methods do not depend on any engine.

If multiple engines have been instantiated, never call `TVirtualMC::ProcessRun(...)` or other steering methods on the engine since that would bypass the `TMCManager`
