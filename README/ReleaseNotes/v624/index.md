% ROOT Version 6.24 Release Notes
% 2020-05-19
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.24/00 is scheduled for release in November 2020.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Guilherme Amadio, CERN/SFT,\
 Bertrand Bellenot, CERN/SFT,\
 Jakob Blomer, CERN/SFT,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Massimiliano Galli, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Hadrien Grasland, IJCLab/LAL,\
 Enrico Guiraud, CERN/SFT,\
 Claire Guyot, CERN/SFT,\
 Stephan Hageboeck, CERN/SFT,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak-Tadel, UCSD/CMS,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, CERN/SFT and UPV,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Andrea Sciandra, SCIPP-UCSC/Atlas, \
 Oksana Shadura, UNL/CMS,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Matevz Tadel, UCSD/CMS,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,\
 Stefan Wunsch, CERN/SFT

## Deprecation and Removal


## Core Libraries

Due to internal changes required to comply with the deprecation of Intel TBB's `task_scheduler_init` and related
interfaces in recent TBB versions, as of v6.24 ROOT will not honor a maximum concurrency level set with
`tbb::task_scheduler_init` but will require instead the usage of `tbb::global_control`:

```cpp
  //tbb::task_scheduler_init init(2); // does not affect the number of threads ROOT will use anymore

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 2);
  ROOT::TThreadExecutor p1;  // will use 2 threads
  ROOT::TThreadExecutor p2(/*nThreads=*/8); // will still use 2 threads
```

Note that the preferred way to steer ROOT's concurrency level is still through `[ROOT::EnableImplicitMT](https://root.cern/doc/master/namespaceROOT.html#a06f2b8b216b615e5abbc872c9feff40f)` or by passing the appropriate parameter to executors' constructors, as in `[TThreadExecutor::TThreadExecutor](https://root.cern/doc/master/classROOT_1_1TThreadExecutor.html#ac7783d52c56cc7875d3954cf212247bb)`.

See the discussion at [ROOT-11014](https://sft.its.cern.ch/jira/browse/ROOT-11014) for more context.

## I/O Libraries


## TTree Libraries

- TTree now supports the inclusion of leaves of types `long` and `unsigned long` (and therefore also `std::size_t` on most systems) also for branches in "leaflist mode". The corresponding leaflist letters are 'G' and 'g'.

## RDataFrame

- With [ROOT-10023](https://sft.its.cern.ch/jira/browse/ROOT-10023) fixed, RDataFrame can now read and write certain branches containing unsplit objects, i.e. TBranchObjects. More information is available at [ROOT-10022](https://sft.its.cern.ch/jira/browse/ROOT-10022).
- Snapshot now respects the basket size and split level of the original branch when copying branches to a new TTree.
- For some `TTrees`, RDataFrame's `GetColumnNames` method returns multiple valid spellings for a given column. For example, leaf `"l"` under branch `"b"` might now be mentioned as `"l"` as well as `"b.l"`, while only one of the two spellings might have been recognized before.
- Introduce `ROOT::RDF::RunGraphs`, which allows to compute the results of multiple RDataFrames concurrently while sharing the same thread pool. The computation may be more efficient than running the RDataFrames sequentially if an analysis consists of many RDataFrames, which don't have enough data to fully utilize the available resources.
- CSV files can now be opened and processed from HTTP(S) locations
- Certain RDF-related types in the ROOT::Detail and ROOT::Internal namespaces have been renamed, most notably `RCustomColumn` is now `RDefine`. This does not impact code that only makes use of entities in the public ROOT namespace, and should not impact downstream code unless it was patching or reusing internal RDataFrame types.
- Just-in-time compilation of string expressions passed to `Filter` and `Define` now generates functions that take fundamental types by const value (rather than by non-const reference as before). This will break code that was assigning to column values in string expressions: this is an intended side effect as we want to prevent non-expert users from performing assignments (=) rather than comparisons (==). Expert users can resort to compiled callables if they absolutely have to assign to column values (not recommended). See [ROOT-11009](https://sft.its.cern.ch/jira/browse/ROOT-11009) for further discussion.
- RDataFrame action results are now automatically mergeable thanks to the new interface provided by `ROOT::Detail::RDF::RMergeableValue` and derived, introduced in [#5552](https://github.com/root-project/root/pull/5552). A feature originally requested with [ROOT-9869](https://sft.its.cern.ch/jira/browse/ROOT-9869), it helps streamline RDataFrame workflow in a distributed environment. Currently only a subset of RDataFrame actions have their corresponding mergeable class, but in the future it will be possible to extend it to any action through the creation of a new `RMergeableValue` derived class.

## Histogram Libraries


## Math Libraries

- Update the definitions of the physical constants using the recommended 2018 values from NIST.
 - Use also the new SI definition of base units from 2019, where the Planck constant, the Boltzman constant , the elementary electric charge and the Avogadro constant are exact numerical values. See
 <https://en.wikipedia.org/wiki/2019_redefinition_of_the_SI_base_units>. Note that with this new definition the functions `TMath::HUncertainty()`, `TMath::KUncertainty()`,
 `TMath::QeUncertainty()` and `TMath::NaUncertainty()` all return a  `0.0` value.



## RooFit Libraries

### Unbiased binned fits
When RooFit performs binned fits, it takes the probability density at the bin centre as a proxy for the probability in the bin. This can lead to a bias.
To alleviate this, the new class [RooBinSamplingPdf](https://root.cern/doc/v624/classRooBinSamplingPdf.html) has been added to RooFit.

### Improved recovery from invalid parameters
When a function in RooFit is undefined (Poisson with negative mean, PDF with negative values, etc), RooFit can now pass information about the
"badness" of the violation to the minimiser. The minimiser can use this to compute a gradient to find its way out of the undefined region.
This can drastically improve its ability to recover when unstable fit models are used, for example RooPolynomial.

For details, see the RooFit tutorial [rf612_recoverFromInvalidParameters.C](https://root.cern/doc/v624/rf612__recoverFromInvalidParameters_8C.html).

## 2D Graphics Libraries

- Add the method `AddPoint`to `TGraph(x,y)` and `TGraph2D(x,y,z)`. equivalent to `SetPoint(g->GetN(),x,y)`and `SetPoint(g->GetN(),x,y,z)`
- Option `E0` draws error bars and markers are drawn for bins with 0 contents. Now, combined
  with options E1 and E2, it avoids error bars clipping.

## 3D Graphics Libraries


## Geometry Libraries


## Database Libraries


## Networking Libraries

### Multithreaded support for FastCGI
Now when THttpServer creates FastCGI engine, 10 worker threds used to process requests
received via FastCGI channel. This significantly increase a performance, especially when
several clients are connected.

### Better security for THttpServer with webgui
If THttpServer created for use with webgui widgets (RBrowser, RCanvas, REve), it only will
provide access to the widgets via websocket connection - any other kind of requests like root.json
or exe.json will be refused completely. Cobined with connection tokens and https protocol,
this makes usage of webgui components in public networks more secure.

### Enabled WLCG Bearer Tokens support in RDavix
Bearer tokens are part of WLCG capability-based infrastructure with capability-based scheme which uses an infrastructure that describes what the bearer is allowed to do as opposed to who that bearer is. Token discovery procedure are developed according WLCG Bearer Token Discovery specification document (https://github.com/WLCG-AuthZ-WG/bearer-token-discovery/blob/master/specification.md). Short overview:
   1. If the `BEARER_TOKEN` environment variable is set, then the value is taken to be the token contents.
   2. If the `BEARER_TOKEN_FILE` environment variable is set, then its value is interpreted as a filename. The contents of the specified file are taken to be the token contents.
   3. If the `XDG_RUNTIME_DIR` environment variable is set, then take the token from the contents of `$XDG_RUNTIME_DIR/bt_u$ID`(this additional location is intended to provide improved security for shared login environments as `$XDG_RUNTIME_DIR` is defined to be user-specific as opposed to a system-wide directory.).
   4. Otherwise, take the token from `/tmp/bt_u$ID`.

## GUI Libraries

### RBrowser improvments
- central factory methods to handle browsing, editing and drawing of different classes
- simple possibility to extend RBrowser on user-defined classes
- support of web-based geometry viewer
- better support of TTree drawing
- server-side handling of code editor and image viewer widgets
- rbrowser content is fully recovered when web-browser is reloaded
- load of widgets code only when really required (shorter startup time for RBrowser)


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings


## JavaScript ROOT

### Major JSROOT update to version 6
- update all used libraries `d3.js`, `three.js`, `MathJax.js`, openui5
- change to Promise based interface for all async methods, remove call-back arguments
- change scripts names, core scripts name now `JSRoot.core.js`
- unify function/methods naming conventions, many changes in method names
- provide central code loader via `JSROOT.require`, supporting 4 different loading engines
- many nice features and many bug fixes; see JSROOT v6 release notes


## Tutorials


## Class Reference Guide


## Build, Configuration and Testing Infrastructure

- a new cmake variable, `CMAKE_INSTALL_PYTHONDIR`, has been added: it allows customization of the installation directory of ROOT's python modules
- The developer build option `asserts` is introduced to enable/disable asserts via the `NDEBUG` C/CXX flag. Asserts are always enabled for `CMAKE_BUILD_TYPE=Debug` and `dev=ON`. The previous behavior of the builds set via the `CMAKE_BUILD_TYPE` variable has not changed.

The following builtins have been updated:

- VecCore 0.7.0

## PyROOT

- Deprecate `TTree.AsMatrix` in this release and mark for removal in v6.26. Please use instead `RDataFrame.AsNumpy`.
