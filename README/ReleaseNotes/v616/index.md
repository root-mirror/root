% ROOT Version 6.16 Release Notes
% 2018-06-25
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.16/00 is scheduled for release end of 2018.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Kim Albertsson, CERN/ATLAS,\
 Guilherme Amadio, CERN/SFT,\
 Bertrand Bellenot, CERN/SFT,\
 Iliana Betsou, CERN/SFT,\
 Brian Bockelman, UNL,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Gerri Ganis, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Enrico Guiraud, CERN/SFT,\
 Stephan Hageboeck, CERN/SFT,\
 Siddhartha Rao Kamalakara, GSOC, \
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak Tadel, UCSD/CMS,\
 Axel Naumann, CERN/SFT,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Oksana Shadura, UNL,\
 Ravi Kiran Selvam, GSOC, \
 Manos, Stergiadis, GSOC, \
 Matevz Tadel, UCSD/CMS,\
 Yuka Takahashi, Princeton,\
 Massimo Tumolo, Politecnico di Torino,\
 Mohammad Uzair, CERN/SFT, \
 Xavier Valls, CERN/SFT,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,\
 Stefan Wunsch, CERN/SFT

## Deprecation and Removal

### Ruby bindings

The ruby binding has been unmaintained for several years; it does not build with current ruby versions.
Given that this effectively meant that Ruby was dysfunctional and given that nobody (but package maintainers) has complained, we decided to remove it.

### Removal of previously deprecated or disabled packages

The packages `afs`, `chirp`, `glite`, `sapdb`, `srp` and `ios` have been removed from ROOT.
They were deprecated before, or never ported from configure, make to CMake.

### Remove GLUtesselator forward declaration from TVirtualX.h

It was never used in TVirtualX interfaces. If GLUtesselator forward declaration is required, use TGLUtil.h include instead.

## C++ Modules Technology Preview

ROOT has several features which interact with libraries and require implicit
header inclusion. This can be triggered by reading or writing data on disk,
or user actions at the prompt. Often, the headers are immutable and reparsing is
redundant. C++ Modules are designed to minimize the reparsing of the same
header content by providing an efficient on-disk representation of C++ Code.

This is an experimental feature which can be enabled by compiling ROOT with
`-Druntime_cxxmodules=On`. You can read more about the current state of the
feature [here](../../README.CXXMODULES.md).

## Core Libraries

### Automatic activation of implicit multithreading

Provided that the build configuration featured IMT (`-DIMT=ON`), implicit multithreading (`ROOT::EnableImplicitMT()`) is now automatically enabled at runtime for the `root` binary, under the following conditions:
it is  running on a MacOS or Windows machine;
or on Linux machines, if an `X11` (or `Xorg` etc) process is running.
The goal is to only enable IMT on desktop-style computers, that are not shared resources (batch systems, grid nodes, login clusters).

This can be disabled by explicitly passing `-t` or `-T` as command line arguments to turn IMT on or off, and by disabling the automatic IMT setting in `rootrc` (`~/.rootrc` or `system.rootrc`).

### New command line flag "--version" for root

`root --version` now displays ROOT version and build info and quits:

```
ROOT Version: 6.15/01
Built for linuxx8664gcc on Sep 20 2018, 11:04:35
From heads/master@v6-13-04-1273-gea3f4333a2
```

### Fish support for thisroot script

`. bin/thisroot.fish` sets up the needed ROOT environment variables for one of the ROOT team's favorite shells, the [fish shell](https://fishshell.com/).

### Change of setting the compression algorithm in `rootrc`

The previous setting called `ROOT.ZipMode` is now unused and ignored.
Instead, use `Root.CompressionAlgorithm` which sets the compression algorithm according to the values of [ECompression](https://root.cern/doc/master/Compression_8h.html#a0a7df9754a3b7be2b437f357254a771c):

* 0: use the default value of `R__ZipMode` (currently selecting ZLIB)
* 1: use ZLIB (the default until 6.12 and from 6.16)
* 2: use LZMA
* 3: legacy, please don't use
* 4: LZ4

### TRef

* Improve thread scalability of `TRef`. Creating and looking up a lot of `TRef` from the same `processID` now has practically perfect weak scaling.

### Parallelism
* Upgrade the built-in TBB version to 2019_U1.

### Type System
* Upgrade the `TClass::GetMissingDictionaries` method to support `std::unique_ptr`, `std::array` and `std::tuple` without getting trapped in the internal STL implementation details.

## I/O Libraries

* To allow for increase run-time performance and increase thread scalability the override ability of `TFile::GetStreamerInfoList` is replaced by an override of `TFile::GetStreamerInfoListImp` with updated return type and arguments.   If a class override `TFile::GetStreamerInfoList` you will now see a compilation error like:

```
/opt/build/root_builds/rootcling.cmake/include/TSQLFile.h:225:19: error: declaration of 'GetStreamerInfoList' overrides a 'final' function
virtual TList *GetStreamerInfoList();
^
/opt/build/root_builds/rootcling.cmake/include/TFile.h:231:24: note: overridden virtual function is here
virtual TList      *GetStreamerInfoList() final; // Note: to override behavior, please override GetStreamerInfoListImpl
^
```

Instead you need to override the protected method:

```
InfoListRet GetStreamerInfoListImpl(bool lookupSICache);
```

which can be implemented as

```
InfoListRet DerivedClass::GetStreamerInfoListImpl(bool /*lookupSICache*/) {
ROOT::Internal::RConcurrentHashColl::HashValue hash;
TList *infolist = nullptr;
//
// Body of the former Derived::GetStreamerInfoList with the
// return statement replaced with something like:

// The second element indicates success or failure of the load.
// (i.e. {nullptr, 0, hash} indicates the list has already been processed
//  {nullptr, 1, hash} indicates the list failed to be loaded
return {infolist, 0, hash};
}
```

See `TFile::GetStreamerInfoListImpl` implementation for an example on how to implement the caching.

* ZLIB (with compression level 1) is now the default compression algorithm for new ROOT files (LZ4 was default compression algorithm in 6.14). Because of reported "corner cases" for LZ4, we are working on the fix to be landed in a next release and return back LZ4 as a default compression algorithm.

* Introducing a possibility for ROOT to use generic compression algorithm/level/settings, by introducing new generic class RCompressionSetting together with new structs ELevel (compression level), EDefaults (default compression settings) and EAlgorithm (compression algorithm). These changes are the first step in generalization of setup of ROOT compression algorithm. It also provides correctness of resolution of compression level and compression algorithm from defined ROOT compression settings:

```
Attaching file hsimple.root as _file0...
(TFile *) 0x55f79b0e6360
root [1] _file0->GetCompressionAlgorithm()
(int) 1
root [2] _file0->GetCompressionLevel()
(int) 1
root [3] _file0->GetCompressionSettings()
(int) 101
root [4]
```

## TTree Libraries
### RDataFrame
  - Optimise the creation of the set of branch names of an input dataset, doing the work once and caching it in the RInterface.
  - Add [StdDev](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a482c4e4f81fe1e421c016f89cd281572) action.
  - Add [Display](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#aee68f4411f16f00a1d46eccb6d296f01) action and [tutorial](https://github.com/root-project/root/blob/master/tutorials/dataframe/df024_Display.C).
  - Add [Graph](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a804b466ebdbddef5c7e3400cc6b89301) action and [tutorial](https://github.com/root-project/root/blob/master/tutorials/dataframe/df021_createTGraph.C).
  - Improve [GetColumnNames](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a951fe60b74d3a9fda37df59fd1dac186) to have no redundancy in the returned names.
  - Add [Kahan Summation tutorial](https://github.com/root-project/root/blob/master/tutorials/dataframe/df022_useKahan.C) to subscribe a Kahan summation action to the RDataFrame.
  - Add [Aggregate tutorial](https://github.com/root-project/root/blob/master/tutorials/dataframe/df019_Cache.C) [Python version](https://github.com/root-project/root/blob/master/tutorials/dataframe/df019_Cache.py).
  - Add [Cache tutorial](https://github.com/root-project/root/commit/cd3e2fdc4baa99111f57240bf8012dcc5f1b5dc6).
  - Fix ambiguous call on Cache() with one or two columns as parameters.
  - Add [GetFilterNames](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a25026681111897058299161a70ad9bb2).
  - Improve RDF node ownership model. The net effect is that users do not have to worry about keeping the first node of a computation graph in scope anymore.
  - Make RResultPtr copy/move-assignable and copy/move-constructible.
  - Add [GetColumnType](https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#ad3ccd813d9fed014ae6a080411c5b5a8a) utility method to query the type of a RDF column (returned as a string).
  - Add [PassAsVec](https://root.cern/doc/master/namespaceROOT_1_1RDF.html#a1ecc8a41e8f12e65e1bf0d2e65aec36d) helper function.
  - Add [SaveGraph](https://root.cern/doc/master/namespaceROOT_1_1RDF.html#adc17882b283c3d3ba85b1a236197c533) helper function to write out the RDF computation graph as a graphviz file.
  - Add a [tutorial for RDataFrame helper functions](https://root.cern/doc/master/df020__helpers_8C.html).
  - Throw if name of a custom column is not a valid C++ name.
  - Allow every RDataFrame variable be cast to a common type `ROOT::RDF::RNode`.
  - Speed up just-in-time compilation (and therefore runtime) of Snapshots with a large number of branches.
  - Create names for histograms and graphs based on the input columns if no model is provided.
  - RCutFlowReport can print cumulative efficiency of cuts.
  - Reading and writing of columns holding `vector<bool>` instances and `bool` C arrays.
  - Support `rdfentry_` and `rdfslot_` implicit columns (`tdfentry_` and `tdfslot_` kept for backwards compatibility).
  - Remove `RDataFrame` from the 32-bit builds.
  - Speed up interpreted usage of RDataFrame (i.e. in macros or from ROOT prompt) by removing certain cling runtime safety checks.
  - Streamline and better document usage of multi-thread RDataFrame: edge cases in which processing of an event could start
    before processing of another event finished have been removed, making it easier for user to write safe parallel RDF operations.
    See the [relevant documentation](https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html#parallel-execution) for more information.

### TTreeProcessorMT
  - Parallelise search of cluster boundaries for input datasets with no friends or TEntryLists. The net effect is a faster initialization time in this common case.
  - Handle gracefully the presence of chains the files associated to which are corrupted.
  - Reduce number of expensive `TChain::LoadTree` calls by spawning nested TBB tasks to ensure clusters of a given file will be most likely processed by the same thread.

### TTree
  - TTrees can be forced to only create new baskets at event cluster boundaries.
    This simplifies file layout and I/O at the cost of memory.  Recommended for
    simple file formats such as ntuples but not more complex data types.  To
    enable, invoke `tree->SetBit(TTree::kOnlyFlushAtCluster)`.

## Histogram Libraries


## Math Libraries

### [Clad](https://github.com/vgvassilev/clad)
  - Enable experimental automatic differentiation techniques to compute
    derivatives and gradients of functions. Automatic differentiation is
    superior to the slow symbolic or often inaccurate numerical differentiation.
    It uses the fact that every computer program can be divided into a set of
    elementary operations (-,+,*,/) and functions (sin, cos, log, etc). By
    applying the chain rule repeatedly to these operations, derivatives of
    arbitrary order can be computed.
  - Implement experimental `TFormula::GradientPar` derivative engine which
    employs clad.

### VecOps
  - Add `All` helper: return true if all of the elements equate to true, return false otherwise.
  - Add `Any` helper: return true if any of the elements equates to true, return false otherwise.
  - Add `ArgSort` helper: return an RVec of indices that sort the vector.
  - Add `Combinations` helper which can:
    - return the indices which represent all combinations of the elements of two vectors.
    - return the indices which represent all unique n-tuple combinations of the elements of a given vector.
  - Add `Intersect` helper: return the intersection of elements of two RVecs.
  - Add `Nonzero` helper: return the indices of the elements which are not zero
  - Add `Reverse` helepr: return copy of reversed RVec.
  - Add `Sort` helper: return copy of vector with elements sorted in ascending order (also according to a user defined predicate)
  - Add `Take` helper which can:
    - return elements of a RVec at given indices.
    - return first elements or last elements of an RVec.
  - Add `Where` helper which can:
    - return the elements of v1 if the condition c is true and v2 if the condition c is false.
    - return the elements of v1 if the condition c is true and sets the value v2 if the condition c is false.
    - return the elements of v2 if the condition c is false and sets the value v1 if the condition c is true.
    - return a vector with the value v2 if the condition c is false and sets the value v1 if the condition c is true.

## RooFit Libraries
  - Add value printer for RooAbsArg and daughters.
  - Add a Python version for the majority of the Tutorials.

## TMVA Library

### Deep Learning

This release contains several fixes and improvement for the `MethodDL`. The `MethodDL` is also now the recommended class to use for Deep Learning in TMVA and is replacing the previous existing
`MethodDNN`, which is still available, but it has a limited functionality and it supports only dense layer.

The new features of `MethodDL` are:

 - Support training and evaluation of Convolutional layer on GPU
 - Several ML optimizers are now included and they can be used in addition to SGD. These are ADAM (the new default), ADAGRAD,
RMSPROP, ADADELTA. A new option, *Optimizer* has been added in the option string used to define the training strategy options.
 - Add support for regression in MethodDL
 - Use single precision (float) types as the fundamental type for the neural network architecture. Double precision could be enabled, but it will require recompiling TMVA. 
 - Support inference (network evaluation) in batch mode in addition to single event. Batch mode evaluation is now the default when used within the `TMVA::Factory` class (i.e. when calling
`Factory::TestAllMethod()` or `Factory::EvaluateAllMethods()`
 - Support splitting the overall training data in Train and Validation data. The train data is used for finding the optimal network weight and the validation data is used to monitor the validation
error. The weights which are giving a minimal validation error will be stored. For the splitting a new option, *ValidationSize* has been added to the global options for `MethodDL`.
The same option is also available in the `PyKeras` method of `PyMVA`
 - The fast tanh implementation from VDT is now used as activation function when training the network on CPU.
 - Using `Cblas` from the GSL library is supported for CPU training when no other Blas libraries are found. However, it is strongly recommended, to use an optimized Blas implementation such as `libopenblas`, that is
available in cvmfs.
 - Add several performance optimizations for both CPU and GPU versions of `MethodDL`. 


### Other New TMVA Features

 - Add a new option to the `DataLoader` to switch off computation of correlation matrix. The new option is called *CalcCorrelations* and it should be used when a large number of input variables are
  provided, otherwise TMVA will spend a long time in setting up the data set before training. 
 
 - Build configuration:
   - Add new cmake flags, `tmva-cpu` and `tmva-gpu`, which can be used to swicth on/off the CPU and GPU (based on CUDA) implementations of the TMVA Deep Learning module. `tmva-cpu` is enabled by
   default if a Blas or CBlas library is found in the system. `tmva-gpu` is enabled when the cmake flag `cuda` is enabled and a compatible Cuda library is found. 
   enabled if the corre
   - Add possibility to independently configure building of optional pymva part of tmva with flag `-Dpymva=ON|OFF`.

 - New Cross Validation features:
   - Add stratified splitting for cross validation.
   - New plotting option in cross validation, average ROC curve.

 - Bugfixes:
   - Fix bug in BDT training with imt=on
   - Improved handling of large event numbers in cross validation using deterministic splitting

 - Documentation:
   - Update TMVA Users' guide

## 2D Graphics Libraries

  - Highlight mode is implemented for `TH1` and for `TGraph` classes. When
    highlight mode is on, mouse movement over the bin will be represented
    graphically. Histograms bins or graph points will be highlighted. Moreover,
    any highlight emits signal `TCanvas::Highlighted()` which allows the user to
    react and call their own function. For a better understanding see also
    the tutorials `$ROOTSYS/tutorials/hist/hlHisto*.C` and
    `$ROOTSYS/tutorials/graphs/hlGraph*.C` .
  - Implement fonts embedding for PDF output. The "EmbedFonts" option allows to
    embed the fonts used in a PDF file inside that file. This option relies on
    the "gs" command (https://ghostscript.com).

    Example:

~~~ {.cpp}
   canvas->Print("example.pdf","EmbedFonts");
~~~
  - In TAttAxis::SaveAttributes` take into account the new default value for `TitleOffset`.
  - When the histograms' title's font was set in pixel the position of the
    `TPaveText` containing the title was not correct. This problem was reported
    [here](https://root-forum.cern.ch/t/titles-disappear-for-font-precision-3/).
  - In `TGraph2D` when the points were all in the same plane (along X or Y) at a
    negative coordinate, the computed axis limits were not correct. This was reported
    [here](https://root-forum.cern.ch/t/potential-bug-in-tgraph2d/29700/5).
  - Implemented the drawing of filled polygons in NDC space as requested
    [here](https://sft.its.cern.ch/jira/browse/ROOT-9523)
  - Implement the drawing of filled polygons in NDC space.
  - When drawing a histogram with the automatic coloring options (PMC, PLC etc ...)
    it was easy to forget to add a drawing option. This is now fixed. If no drawing
    option is specified the default drawing option for histogram is added.
  - When drawing a `TGraph` if `TH1::SetDefaultSumw2()` was on, then the underlying
    histogram used to draw the `TGraph` axis was created with errors and therefore
    the histogram painter painted it with errors which was a non sense in that
    particular case. This is now fixed. It was discussed
    [here](https://root-forum.cern.ch/t/horizontal-line-at-0-on-y-axis/30244/26)
  - Add `TGraph2D::GetPoint`, with similar interface and behaviour as `TGraph::GetPoint`

## 3D Graphics Libraries


## Web Graphics Libraries

  - Introduce web-based output for TCanvas. When ROOT compiled with cxx14 support `cmake -Dcxx14=ON ..` and
    started with web option like `root --web hsimple.C`, TCanvas will be displayed in web browser.
  - Most of histograms and graphs classes are supported. See 
    [JavaScript ROOT website](https://root.cern/js/latest/examples.htm) for list of suported classes
  - Also some classes with customize painters may be working - like TTree::Draw() into TParallelCoord object. These kind of
    objects handled on server side, web browser used just for display of produced primitives like polylines or text          


## Geometry Libraries


## Database Libraries


## Networking Libraries


## GUI Libraries


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings

### PyROOT
  - Fixed support for templated functions when in need of:
    - typedef resolution (`Foo<Float_t>` -> `Foo<float>`)
    - namespace addition (`Foo<vector<float>>` -> `Foo<std::vector<float>>`)
    - full name completion (`Foo<std::vector<float>>` -> `Foo<std::vector<float, std::allocator<float>>>`)

### Experimental PyROOT
  - Added pythonisations for `TTree` and its subclasses (e.g. `TChain`, `TNtuple`)
    - Pythonic iterator (`for event in tree:`)
    - Access tree branches as attributes (`mytree.mybranch`)
    - `TTree::Branch` pythonisation
    - `TTree::SetBranchAddress` pythonisation
  - Added pythonisations for `TDirectory` and its subclasses (e.g `TFile`, `TDirectoryFile`)
    - Access directories/objects in `TDirectory`/`TDirectoryFile`/`TFile` as attributes
    (`mydir1.mydir2.myhist`, `myfile.myhist`, `myfile.mydir.myhist`)
    - `TDirectory::Get` pythonisation
    - `TDirectory::WriteObject` pythonisation
    - `TFile::Open` pythonisation
  - Added pretty printing generic pythonisation for all classes
  - Added interoperability with NumPy arrays for STL vectors and `RVec`s (zero-copy wrapping of
  vectors and `RVec`s into NumPy arrays)

### Jupyter Notebook Integration
  - Make sure the ROOT C++ Jupyter kernel runs with the same Python version (major and minor) that ROOT
  was built with.
  - Make the Jupyter server started with `root --notebook` listen on all interfaces. This can be useful
  if the user wants to connect to the server remotely. It also fixes an issue observed when starting
  the Jupyter server inside a Docker container.

## JavaScript ROOT
  - Support of TWebCanvas functionality. Code for ROOT 6.16 will
    be maintained in v6-16-00-patches branch of JSROOT repository. 
  - Singificant speed up (factor 10) when drawing canvas with many subpads
  - Many small improvements and bug fixes, see JSROOT release notes for v5.4.2 - v5.6.2 

## Tutorials
  - Refurbish text in the `RDataFrame` tutorials category.

## Command line tools
  - Fixed `rooteventselector` when both applying a cut (based on branch values) and selecting only
  a subset of the branches. Previously, the size of the output file was bigger than expected.

## Class Reference Guide


## Build, Configuration and Testing Infrastructure
