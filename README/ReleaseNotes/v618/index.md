% ROOT Version 6.18 Release Notes
% 2019-05-28
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.18/00 is scheduled for release in June, 2019.

For more information, see [http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Kim Albertsson, CERN/ATLAS,\
 Guilherme Amadio, CERN/SFT,\
 Bertrand Bellenot, CERN/SFT,\
 Iliana Betsou, CERN/SFT,\
 Jakob Blomer, CERN/SFT,\
 Brian Bockelman, Nebraska,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Javier Cervantes Villanueva, CERN/SFT,\
 Olivier Couet, CERN/SFT,\
 Alexandra Dobrescu, CERN/SFT,\
 Giulio Eulisse, CERN/ALICE,\
 Gerri Ganis, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Enrico Guiraud, CERN/SFT,\
 Stephan Hageboeck, CERN/SFT,\
 Jan Knedlik, GSI,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak-Tadel, UCSD/CMS,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, Bicocca/SFT,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Henry Schreiner, Princeton,\
 Oksana Shadura, Nebraska,\
 Simon Spies, GSI,\
 Yuka Takahashi, Princeton and CERN/SFT,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Matevz Tadel, UCSD/CMS,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,\
 Zhe Zhang, Nebraska,\
 Stefan Wunsch, CERN/SFT


## Deprecation and Removal

### Deprecated packages

The Virtual Monte Carlo (VMC) interfaces have been deprecated for this release
and will be removed in a future release. It is no longer built by default, but
can still be enabled with the option `-Dvmc=ON` in the CMake configuration phase.
A standalone version of VMC is being developed at [https://github.com/vmc-project/vmc](https://github.com/vmc-project/vmc)
to replace the deprecated version in ROOT.

### Removed packages

Support for the following optional components of ROOT has been removed:

 * afdsmgrd (Dataset manager for PROOF-based analysis facilities)
 * bonjour (Avahi/Bonjour/Zeroconf)
 * castor (CERN Advanced STORage manager)
 * geocad (OpenCascade)
 * globus (Globus authentication)
 * hdfs (Hadoop Distributed File System)
 * krb5 (Kerberos 5 authentication)
 * ldap (OpenLDAP authentication)
 * memstat (legacy memory statistics utility)
 * qt, qtgsi, qtroot (Qt4-based GUI components)
 * rfio (Remote File IO for CASTOR)
 * table (libTable contrib library)

In addition, the following deprecated parts of ROOT components have been
removed:

 * PROOF's PQ2 module
 * `THttpServer::ExecuteHttp()` and `THttpServer::SubmitHttp` from `THttpServer`

### Other changes

The deprecation of the GraphViz integration has been reverted since the code is
still in use.

The ODBC interface, deprecated in ROOT 6.16, is no longer deprecated in ROOT 6.18.
It is the main option to support databases on Windows, so the decision to deprecate
it was reverted.

The `xft` option has been merged into `x11` and is no longer used (its value is
now ignored by ROOT).

## Preprocessor deprecation macros
### Deprecated Classes
  * `R__SUGGEST_ALTERNATIVE("Suggestion text")` macro allows to suggest alternatives to classes. It must be used after the class definition and before the final semicolon:
```{.cpp}
class DoNotUseClass {
} R__SUGGEST_ALTERNATIVE("Use ... instead.");
```
It is activated by the preprocessor defines `R__SUGGEST_NEW_INTERFACE`. The former is useful when deprecation warnings should be activated/deactivated at global level, for example for an entire project. This could be done by defining `R__SUGGEST_NEW_INTERFACE` in the build system. 
If the warning needs to be confined within single translation units, irrespective of the definition of `R__SUGGEST_NEW_INTERFACE`, the `R__ALWAYS_SUGGEST_ALTERNATIVE` macro can be used:
```{.cpp}
#ifndef DONOTUSECLASS_H
#define DONOTUSECLASS_H

class DoNotUseClass {
} R__ALWAYS_SUGGEST_ALTERNATIVE("Use ... instead.");

#endif
```
### Deprecated Functions
The same macro as for classes can be used for functions:
```{.cpp}
TIterator* createIterator() const
R__SUGGEST_ALTERNATIVE("begin(), end() and range-based for loops.") {
 return makeLegacyIterator();
}
```


### I/O Libraries

* The deprecrated `I/O` plugins for  `HDFS`, `Castor` and `RFIO` have been removed.

### THttpServer classes

The following methods were deprecated and removed:

   * `Bool_t THttpServer::SubmitHttp(THttpCallArg *arg, Bool_t can_run_immediately = kFALSE, Bool_t ownership = kFALSE)`
   * `Bool_t THttpServer::ExecuteHttp(THttpCallArg *arg)`
   * `Bool_t TRootSniffer::Produce(const char *path, const char *file, const char *options, void *&ptr, Long_t &length, TString &str)`
   * `TString THttpCallArg::GetPostDataAsString()`
   * `void THttpCallArg::FillHttpHeader(TString &buf, const char *header = nullptr)`
   * `void THttpCallArg::SetBinData(void *data, Long_t length)`

The methods could be replaced by equivalent methods with other signature:

   * `Bool_t THttpServer::SubmitHttp(std::shared_ptr<THttpCallArg> arg, Bool_t can_run_immediately = kFALSE)`
   * `Bool_t THttpServer::ExecuteHttp(std::shared_ptr<THttpCallArg> arg)`
   * `Bool_t TRootSniffer::Produce(const std::string &path, const std::string &file, const std::string &options, std::string &res)`
   * `const void *THttpCallArg::GetPostData() const`
   * `Long_t THttpCallArg::GetPostDataLength() const`
   * `std::string THttpCallArg::FillHttpHeader(const char *header = nullptr)`
   * `void THttpCallArg::SetContent(std::string &&cont)`

### Core Libraries
  * The `TStringLong` class is deprecated. Please use `std::string` (or, if needeed, `TString`) instead.


## Core Libraries
  - Prevent usage of non integer class id in `ClassDef(Inline)` macros with an error prompted at dictionary generation or compile time.

## I/O Libraries

* Added simpler way to retrieve object from `TDirectory` and `TFile`:
~~~ {.cpp}
auto obj = directory->Get<MyClass>("some object");
~~~

* Added support for read-only `TMemFile`s.

### TNetXNGFile
Added necessary changes to allow [XRootD local redirection](https://github.com/xrootd/xrootd/blob/8c9d0a9cc7f00cbb2db35be275c35126f3e091c0/docs/ReleaseNotes.txt#L14)
  - Uses standard VectorReadLimits and does not query a XRootD data server (which is unknown in local redirection), when it is redirected to a local file
  - Adds a new constructor with a `const char *lurl` to `TNetXNGFile` and passes it to `TFile`, if set. This allows redirection to files that have a different name in the local file system and is important to allow derivation (for example to `TAlien` and `TJAlienFile`) while still keeping functionality via `TArchiveFile` when the file name in the local file system does not match `*.zip`

### TBufferJSON
Add possibility to convert STL `std::map`, `std::multimap`, `std::unordered_map`,
`std::unordered_multimap` classes into JSON object. This only possible when key typename
is `std::string` (or compatible) and contains only valid JSON identifiers. By default these classes converted
into JSON array of `std::pair` objects. To enable new feature, compact parameter should be 5:

~~~ {.cpp}
std::map<std::string,int> obj;
obj["name1"] = 11;
obj["name1"] = 22;
auto json = TBufferJSON::ToJSON(&obj, 5);
// {"_typename": "map<string,int>", "name1": 11, "name2": 22}
auto dflt_json = TBufferJSON::ToJSON(&obj);
// [{"$pair" : "pair<string,int>", "first" : "name1", "second" : 11}, {"$pair" : "pair<string,int>", "first" : "name2", "second" : 22}]
~~~

Also one could put "JSON_object" string in class-member comments to enable this feature:

~~~ {.cpp}
class Container {
   int field{5};
   std::unordered_map<std::string,double> data;  ///< JSON_object indicates conversion
};
~~~

Now one could disable storage of type information - `_typename` field. For that compact parameter
has to include value 100. Be aware that such JSON representation may not be recognized by JSROOT.
Maximal compression of JSON can be achieved now with compact parameter 128 = 100 + 20 + 5 + 3:
   3 - remove all spaces and new lines
   5 - convert map->object (when applicable)
   20 - special compression of large arrays (auto-detected in JSROOT)
   100 - suppressing `_typename` for all classes


## TTree Libraries

### RDataFrame
  - Use TPRegexp instead of TRegexp to interpret the regex used to select columns
    in the invocation of `Cache` and `Snapshot`. This allows usage of a larger set
    of regular expressions to specify which columns should be written out.
  - Speed up jitting of Filter and Define expressions passed as string
  - Speed up event loop, improve scaling in the presence of a large amount of Defines
  - Allow Filter expressions to return types convertible to bool, rather than only bool
  - Add `GetNSlots` method to easily retrieve the number of slots that will be used by
    `DefineSlot`, `ForeachSlot`, `OnPartialResultSlot`, ...
  - Add support for TTrees/TChains with TEntryLists (currently only for single-thread event loops)
  - Add `HasColumn` method to check whether a column is available to a given RDF node
  - PyROOT: add `AsRNode` helper function to convert RDF nodes to the common RNode type
  - PyROOT: add `AsNumpy` method to export contents of a RDataFrame as a dictionary of numpy arrays
    * Such dictionary of numpy arrays can also be used to create a pandas DataFrame
  - Experimental PyROOT: add `MakeNumpyDataFrame` factory to process data owned by numpy arrays with RDataFrame
  - The `Stats` method has been added, allowing to retrieve a `TStatistic` object filled with the values of a column and, optionally, the values of a second column to be used as weights.

### TLeafF16 and TLeafD32
  - New leaf classes allowing to store `Float16_t` and `Double32_t` values using the truncation methods from `TBuffer`
    (See for example `TBuffer::WriteDouble32`)
  - The new types can be specified using the type characters `'f'` (`Float16_t`) and `'d'` (`Double32_t`)
  - It is also possible to specify a range and a number of bits to be stored using the syntax from `TStreamerElement::GetRange`.
    Therefore the range and bits specifier has to be attached to the type character.
  - All functionalities of the other datatypes have been reimplemented.
  - The documentation of `TTree` and `TBuffer` has been updated accordingly.
  - The following example shows how to use the new features:
~~~ {.cpp}
Float16_t  floatVal;
Float16_t  floatArray[7];
Double32_t doubleVal;
Double32_t doubleArray[5];
TTree *tree = new TTree("tree", "An example tree using the new data types");
tree->Branch("floatVal",   &floatVal,    "floatVal/f");               // Float16_t value with default settings
tree->Branch("floatArray",  floatArray,  "floatArray[7]/f[0,100]");   // Float16_t array with range from 0 to 100
tree->Branch("doubleVal",  &doubleVal,   "doubleVal/d[0,1000,20]");   // Double32_t value with range from 0 to 1000 and 20 bits
tree->Branch("doubleArray", doubleArray, "doubleArray[5]/d[0,0,18]"); // Double32_t array without range and 18 bits
~~~

### Bulk I/O
  - The new `TBulkBranchRead` class (inside the `ROOT::Experimental::Internal` namespace) provides
    a mechanism for reading, in a single library call, many events' worth of simple data (primitive types,
    arrays of primitives, split structures) stored in a `TTree`.  This allows for extremely fast delivery
    of event data to the process.  This is meant as an internal interface that allows the ROOT team to
    implement faster high-level interface.
  - The `TTreeReaderFast ` class (inside the `ROOT::Experimental::Internal` namespace) provides a simple
    mechanism for reading ntuples with the bulk IO interface.

## Histogram Libraries

### TH1
  - Add a search range to the `TH1::FindFirstBinAbove(..)` and `TH1::FindLastBinAvove(..)` functions

### TH2Poly
  - Add implementation of `SetBinError` and fix a bug in `GetBinError` in case of weighted events.

### TF1
  - The implementation of `TF1::GetX` has been improved. In case of the presence of multiple roots, the function will return the root with the lower x value. In case of no-roots a NaN will be returned instead of returning a random incorrect value.

### TKDE
  - Add support for I/O


## Math Libraries
  - Add `TComplex` value printer for printing the value of object at the root prompt and in python
  - Add to the documentation of `TLorentzVector` a link to `ROOT::Math::LorentzVector`, which is a superior tool.
  - Add new implementation of `TStatistic::Merge` able to deal silently with empty TStatistic objects. This implementation is useful when filling TStatistics with one of ROOT's implicitly parallelised utilities such as `RDataFrame` or `TThreadExecutor`.
  - Add `T RVec<T>::at>(size_t, T)` method to allow users to specify a default value to be returned in case the vector is shorter than the position specified. No exception is thrown.
  - Add the `Concatenate` helper to merge the content of two `RVec<T>` instances.
  - Generalise the `VecOps::Map` utility allowing to apply a callable on a set of RVecs and not only to one.
  - Add the `DeltaR2`, `DeltaR` and `DeltaPhi` helpers for RVec.
  - Add the `InvariantMass(es)` helpers computing the invariant mass from particle kinematics stored in RVecs.
  - Add the `Max`, `Min`, `ArgMax`, and `ArgMin` helpers for RVec.
  - Add the `Construct` helper to build an `RVec<T>` starting from N `RVec<P_i>`, where a constructor `T::T(P_0, P_1, ..., P_Nm1)` exists.
  - Experimental PyROOT: Add `AsRVec` helper to adopt memory owned by numpy arrays with RVecs

### [Clad](https://github.com/vgvassilev/clad)
  - Upgrade Clad to 0.5 The new release includes some improvements in both
    Forward and Reverse mode:
    * Extend the way to specify a dependent variables. Consider function,
      `double f(double x, double y, double z) {...}`, `clad::differentiate(f, "z")`
      is equivalent to `clad::differentiate(f, 2)`. `clad::gradient(f, "x, y")`
      differentiates with respect to `x` and `y` but not `z`. The gradient results
      are stored in a `_result` parameter in the same order as `x` and `y` were
      specified. Namely, the result of `x` is stored in `_result[0]` and the result
      of `y` in `_result[1]`. If we invert the arguments specified in the string to
      `clad::gradient(f, "y, x")` the results will be stored inversely.
    * Enable recursive differentiation.
    * Support single- and multi-dimensional arrays -- works for arrays with constant
      size like `double A[] = {1, 2, 3};`, `double A[3];` or `double A[1][2][3][4];`

## RooFit Libraries
### RooJohnson PDF
The Johnson SU PDF has been added to RooFit. It comes with an analytical integral and a generator function,
which make it superior (faster and more accurate) than implementing it manually with an interpreted/compiled formula.

### HistFactory
hist2workspace performance optimisations. For a large, ATLAS-style Higgs-->bb workspace with > 100 systematic uncertainties and more than 10 channels, the run time for converting histograms to a fit model decreases by a factor 11 to 12.

### Faster, STL-like Collections in RooFit
RooFit's collections `RooArgSet` and `RooArgList` have been made more STL-like. The underlying implementation used to be the `RooLinkedList`, but now both collections work with `std::vector`. The collections have an STL-like interface concerning iterators such that iterations over the two collections that looked like
```
TIterator* depIter = intDepList.createIterator() ;
RooAbsArg* arg;
while((arg=(RooAbsArg*)depIter->Next())) {
  ...
}
delete depIter;
```
now look like:
```
for (auto arg : intDepList) {
  ...
}
```
Depending on how many elements are iterated, RooFit will be between 10 and 20% faster if the new iterators are used. Heavily using old iterators might slow it down by 5 to 10%. Iterators in key classes have been updated, such that many workflows in RooFit are 10 - 20% faster.

#### Legacy iterators
The (three kinds) of legacy iterators in RooFit are still supported, such that old code will not break, but they are slower than `begin(), end()` and range-based for loops.

**Important caveat**:
The old RooFit collections could be modified while iterating. The STL-like iterators do not support this (as for a *e.g.* std::vector)! Using the legacy iterators with the new collections (*i.e.* in existing code), mutating the collection is still possible in the following cases:
- Inserting/deleting elements **after** the current iterator.
- Changing an element at a position **other than** the current iterator
- **But not** inserting/deleting before/at the current iterator position. With a debug build (with assertions), the legacy iterators will check that the collection is not mutated. In a release build, elements might be skipped or be iterated twice.

#### Moving away from the slower iterators
The legacy iterators have been flagged with a special deprecation macro that can be used help the user use the recommended ROOT interface. Defining one of the [deprecation macros](#preprocessor-deprecation-macros) (either in a single translation unit or in the build system), and creating a legacy iterator will trigger a compiler warning such as:
```
<path>/RooChebychev.cxx:66:34: warning: 'createIterator' is deprecated: There is a superior alternative: begin(), end() and range-based for loops. [-Wdeprecated-declarations]
  TIterator* coefIter = coefList.createIterator() ;
                                 ^
1 warning generated.
```


## TMVA

This release provides a consolidation and several fixes of the new machine learning tools provided in TMVA such as the Deep Learning module.
The method `TMVA::Types::kDL` should be used now for building Deep Learning architecture in TMVA, while `TMVA::Types::kDNN` is now deprecated. `TMVA::Types::kDL` provides all the functionality of `TMVA::Types::kDNN`, i.e building fully connected dense layer, but in addition supports building convolutional and recurrent neural network architectures.
These release contains improvements in the `MethodDL` such as:
  - fix droput support for dense layer
  - add protection to avoid returning NaN in the cross-entropy loss function

In addition we have :

  - New `TMVA::Executor` class to control the multi-thread running of TMVA. By default now MT running will be enabled only when `ROOT::EnabledImplicitMT()` is called. But we can take the control of the threads by using `TMVA::gConfig().EnableMT(...)` and `TMVA::gConfig().DisableMT()`


### PyMVA
  - add support when using the Tensorflow backend in Keras to control the number of threads
  - add possibility to control options for configuring GPU running. FOr example we can now set the mode to allocate memory only as needed. This is required when using the new RTX gaming cards from NVIDIA




## 2D Graphics Libraries

  - In the statistics painting for 2D histograms, the central cell of
    the underflow/overflow grid was not properly rendered for very large contents.
    This problem was reported [here](https://root-forum.cern.ch/t/stat-box-for-th2/).
  - The automatic placement of legend now "sees" TMultiGraph and THStack.
  - Improve and simplify the drawing the 2D histogram's option "ARR".
  - The option ARR can be combined with the option COL or COLZ.
  - `TBox::DistancetoPrimitive` and `TBox::ExecuteEvent` now work in log scales (by Jérémie Dudouet).
  - Take the line attributes into account when drawing a histogram with option bar or hbar.
    They were ignored until now.
  - The new draw option MIN0 makes same effect as gStyle->SetHistMinimumZero(1), but can be specified
    individually for each histogram.
  - Improve the line clipping when a histogram is drawn with option "L". The following
    example shows the improvement.
~~~ {.cpp}
  auto h = new TH1F("h","h",5,0.5,5.5);
  h->SetBinContent(1,100000);
  h->SetBinContent(2,10000);
  h->SetBinContent(3,1000);
  h->SetBinContent(4,100);
  h->SetBinContent(5,10);
  h->SetMinimum(50.);
  h->SetMaximum(40000);
  h->Draw("L*");
  gPad->SetLogy();
~~~
  - `ChangeLabel` is now available for alphanumeric axis.
  - Implement transparency for lines, texts and markers in the TeX output.

## 3D Graphics Libraries

  - Make sure a TF3 is painted the same way in GL and non-GL mode.
    The mismatch was reported in [this post](https://root-forum.cern.ch/t/how-to-specify-the-level-value-in-isosurface-drawing-with-tf3-and-gl/32179)

## Geometry Libraries


## Database Libraries

The CMake module `FindOracle.cmake` was updated to support version 18.x
of the Oracle client libraries.

## Networking Libraries


## GUI Libraries


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings


## JavaScript ROOT

### New functionality from 5.7.0 release

   - Add support of `TProfile2Poly` class
   - Add support of `TGeoOverlap` class
   - Add support of `TGeoHalfSpace` for composites
   - Implement update of `TF2` drawings, see `tutorials/graphics/anim.C`
   - Improve windows handling in flex(ible) layout
   - Provide special widget for object inspector
   - Use `requestAnimationFrame` when do monitoring, improves performance
   - Better position for text in `TH2Poly` drawings
   - Support eve7 geometry viewer - render data generated in ROOT itself
   - Provide initial WebVR support, thanks to Diego Marcos
   - Use `gStyle` attributes to draw histogram title
   - Enable projections drawing also with `TH2` lego plots
   - Many adjustment with new `TWebCanvas` - interactivity, attributes/position updates, context menus
   - Upgrade three.js 86 -> 102, use `SoftwareRenderer` instead of `CanvasRenderer`
   - Upgrade d3.js 4.4.4 -> 5.7.0
   - Fix - support clipping for tracks and points in geo painter
   - Fix - drawing of TGeoNode with finder
   - Fix - key press events processed only in active pad (ROOT-9128)
   - Fix - use X0/Y0 in xtru shape, thanks to @altavir


### New files location

JSROOT sources were moved from `etc/http/` into `js/` subfolder in ROOT sources tree.
OpenUI5 files were moved to `ui5/` subfolder. After ROOT compilation they can be found in
`$ROOTSYS/js/` and `$ROOTSYS/ui5/` subfolders respectively.


## Tutorials
  - Add `RSqliteDS` examples.
  - Make RCsvDS and RLazyDS tutorials standalone, i.e. downloading input csv directly using `TFile::Cp` rather than relying on CMake.

## Class Reference Guide


## Build, Configuration and Testing Infrastructure

### CMake build system requirements and updates

The minimum required version of CMake has been updated to 3.9 or newer to be
able to take advantage of new features such as native support for the CUDA
language, among other things. Please refer to CMake's release notes for further
information.

The method to select the C++ standard has changed. Now the recommended way
to select the C++ standard is via the option `-DCMAKE_CXX_STANDARD=XX`, which
is the idiomatic way to do it in CMake. The old options still work, but have
been deprecated and will be removed in a future release.

Build option descriptions have been updated to indicate which builtins require
an active network connection during the build. You can inspect the list of
options and their descriptions by running `cmake -LH $PWD` in the build
directory.

The build system has been updated to remove most file globbing to improve
the reliability of incremental builds when source files are added or removed.

A new check has been added to make ROOT fail during the configuration step
if incompatible versions of the Python interpreter and its libraries are
selected.

The `all=ON` option now tries to enable more options. Some options had their
default value toggled to disabled, which affected `all=ON`. Now all options
are listed explicitly so that they are enabled regardless of their default
value.

### Builtins

The following builtins had their versions updated for this release:

* VecCore 0.5.2
* Vc 1.4.1
* XRootD 4.8.5
* OpenSSL 1.0.2q
* PCRE 8.42

### Header location and `ROOT_GENERATE_DICTIONARY` / `ROOT_STANDARD_LIBRARY_PACKAGE`

A change in the argument handling of `ROOT_GENERATE_DICTIONARY` and `ROOT_STANDARD_LIBRARY_PACKAGE` might need your attention:
these macros now respect whether a header file was passed with its full relative path (the common case), or with a full path.
The latter allows to find headers at runtime - at the cost of a loss of relocatability: you cannot move the library containing
that dictionary to a different directory, because the header location is stored in the dictionary. This case is used by roottest
but should likely not be used by anything but test suites.

Instead pass relative paths, together with a `-I` option to find the headers, plus setting `ROOT_INCLUDE_PATH` for finding the
headers back at runtime. The ROOT stress suite is now updated to follow this behavior; see for instance the injection of
`$ROOTSYS/test` in `test/stressMathCore.cxx`, allowing ROOT to find the header at runtime, whether interpreted
(`R__ADD_INCLUDE_PATH`) or compiled (`TROOT::AddExtraInterpreterArgs({"-I..."})` before interpreter construction).

If you called `ROOT_GENERATE_DICTIONARY(Dict ${CMAKE_CURRENT_SOURCE_DIR}/subdir/Header1.h LINKDEF LinkDef.h)` then update that
call to `ROOT_GENERATE_DICTIONARY(Dict Header1.h OPTIONS -I subdir LINKDEF LinkDef.h)` *if* the header is usually included as
`#include "Header1.h"`, or to `ROOT_GENERATE_DICTIONARY(Dict subdir/Header1.h LINKDEF LinkDef.h)` *if* the header is usually
included as `#include "subdir/Header1.h"`. I.e. the general rule is: pass to `ROOT_GENERATE_DICTIONARY` (or
`ROOT_STANDARD_LIBRARY_PACKAGE`) the spelling as `#include`ed.

As an important side-effect, `ROOT_GENERATE_DICTIONARY` and thus `ROOT_STANDARD_LIBRARY_PACKAGE` now *require* the header to
be found at configuration time. We have seen too many cases where the header location was mis-stated, and as a consequence,
CMake did not generate the proper dependencies. If the header should not be taken into account for dependencies and / or if
the header will not be found (e.g. the standard library's `vector`) please pass the header through the `NODEPHEADERS` option
to `ROOT_GENERATE_DICTIONARY` or `ROOT_STANDARD_LIBRARY_PACKAGE`.

We believe that this simplification / regularization of behavior, and the additional checks are worth the possible changes
on the user side.


## PyROOT

If the fix or new feature is a pythonization related to a C++ class, the change is added to the respective section above.

### Current PyROOT

- Fix compatibility with Python3.7 (ROOT-9922, ROOT-9871, ROOT-9809)
- Fix lookup for templated methods (ROOT-9789)
- Fix lookup for templated free functions (ROOT-9836)

### Experimental PyROOT

- All pythonisations from current PyROOT already migrated (`TTree` and subclasses, `TDirectory` and subclasses,
    `TCollection` and subclasses, `TObject`, `TClass`, `TString`, `TObjString`, `TIter`, `TStyle`, `TH1`, `TFX`, `TMinuit`, `TVector3`,
    `TVectorT`, `TArray`, `TCollection`, `TSeqCollection`, `TClonesArray`, `TComplex`, `TGraph`, `RooDataHist`) - ROOT-9510
- Cppyy updated to cppyy 1.4.7, cppyy-backend 1.8.1 (clingwrapper), CPyCppyy 1.7.1
  * Includes fixed template support, fixed overload resolution, Windows fixes and other
- Merged Cppyy's patch to support using namespace declarations (PR-3579)
- Add `DeclareCppCallable` decorator, which allows to call Python callables from C++, e.g., in an RDataFrame workflow:
~~~ {.python}
@ROOT.DeclareCppCallable(["float"], "float")

def f(x):
   return 2.0 * x

ROOT.CppCallable.f(21.0)
# Returns 42.0

df = ROOT.ROOT.RDataFrame(4).Define("x", "CppCallable::f(rdfentry_)")

df.AsNumpy()
# Returns {'x': numpy.array([0., 2., 4., 6.], dtype=float32)}
~~~

