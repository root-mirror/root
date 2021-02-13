% ROOT Version 6.20 Release Notes
% 2019-05-29
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.20/00 is scheduled for release in November 2019.

For more information, see:

[http://root.cern](http://root.cern)

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
 Desislava Kalaydjieva, CERN/SFT,\
 Jan Knedlik, GSI,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak-Tadel, UCSD/CMS,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, CERN/SFT and UPV,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Otto Schaile, Uni-Muenchen,\
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

## ROOT

### Splash screen

The venerable splash screen is now disabled by default to make ROOT's startup
faster. Many users already use `root -l` to start ROOT, but this also hides the
useful text banner with version information along with the splash screen. With
this new default, starting up ROOT as just `root` will show only the text banner
instead of the splash screen. The splash screen can still be seen with `root -a`
or in `TBrowser` by opening `Browser Help → About ROOT`.

## Deprecation and Removal
 * rootcling flags `-cint`, `-gccxml`, `-p`, `-r` and `-c` have no effect
   and will be removed. Please remove them from the rootcling invocations.
 * rootcling legacy cint flags `+P`, `+V` and `+STUB` have no effect and will be
   removed. Please remove them from the rootcling invocations.
 * genreflex flag `--deep` has no effect and will be removed. Please remove it
   from the genreflex invocation.
 * rootcling warns if it sees and unrecognized flag (usually coming from the
   CXXFLAGS of the build system). Please remove them from the invocation because
   the warning will become a hard error in the next releases.
 * The empty headers `Gtypes.h` and `Htypes.h` are deprecated. Please include
   `Rtypes.h`
 * TInterpreter::EnableAutoLoading currently does nothing and is deprecated.

### Deprecated packages

### Removed packages

## Core Libraries

* Speed-up startup, in particular in case of no or poor network accesibility, by avoiding
a network access that was used as input to generate a globally unique ID for the current
process.
 * This network access is replaced by a passive scan of the network interface. This
reduces somewhat the uniqueness of the unique ID as the IP address is no longer
guaranteed by the DNS server to be unique.   Note that this was already the case when
the network access (used to look up the hostname and its IP address) failed.


## I/O Libraries

* TFile: A new bit `TFile::kReproducible` was introduced. It can be enabled by
  specifying the `"reproducible"` url option when creating the file:
~~~ {.cpp}
      TFile *f = TFile::Open("name.root?reproducible","RECREATE","File title");
~~~
   Unlike regular `TFile`s, the content of such file has reproducible binary
   content when writing exactly same data. This achieved by writing pre-defined
   values for creation and modification date of TKey/TDirectory objects and null
   value for TUUID objects inside TFile. As drawback, TRef objects stored in such
   file cannot be read correctly.
* Significantly improved the scaling of hadd tear-down/cleanup-phase in the presence
of large number histograms and in the presence of large number of directories.
* TMemFile: Apply customization of minimal block size also to the first block.
* Add renaming rule for instances of the math classes from `genvector` and `smatrix` to
instance for one floating point type (`float`, `double`, `Double32_t`, `Float16_t`) to
instances for any other floating point type.
* Corrected the application of  `I/O customization rules` when the target classes contained
typedefs (in particular `Double32_t`)
* Prevent splitting of objects when a `Streamer Function` was explicitly attached to their
`TClass`.
* In hadd fix verbose level arg parsing
* Allow user to change the type of the content of a TClonesArray.
* Avoid deleted memory access in `MakeProject` and in handling of
`I/O customization rules`.

## TTree Libraries

* Prevent a situation in `TTreeFormula` when stale cached information was re-used.
* Prevent a noticeable memory leak when reading uncompressed TTree.

## Histogram Libraries

* Allow reading v5 TF1 that were stored memberwise in a TClonesArray.

## Math Libraries


## RooFit Libraries


## 2D Graphics Libraries

 * Provide support of NDC coordinates for TArrow.
 * Fix interactive movement of TLine/TArrow objects when NDC coordinates are used
 * Provide TGraph::MovePoints() method
 * New options `RX`and `RY` for TMultiGraph in order to draw reverse axis along X and Y.
 * Combined with the option "Z" the option "CJUST" allows to draw the color palette
   with axis labels justified on the color boundaries (implemented by Otto Schaile).
 * The `TCanvas` Event Status Bar now displays the date and time when the mouse cursor
   is moved over a time axis (implemented by Otto Schaile).
 * Negative values were not painted with option "TEXT" for TH2Poly.


## 3D Graphics Libraries


## Geometry Libraries


## Database Libraries


## Networking Libraries


## GUI Libraries


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings

### Jupyter Notebook Integration
- When starting Jupyter server with `root --notebook arg1 arg2 ...`, extra arguments can be provided.
  All these arguments delivered as is to jupyter executable and can be used for configuration.
  Like server binding to specific host `root --notebook --ip=hostname`
- Remove `c.NotebookApp.ip = '*'` from default jupyter config. One has to provide ip address for server
  binding using `root --notebook --ip=<hostaddr>` arguments
- Now Jupyter Notebooks will use JSROOT provided with ROOT installation. This allows to use notebooks
  without internet connection (offline).


## JavaScript ROOT
- Provide monitoring capabilities for TGeoManager object. Now geomtry with some tracks can be displayed and
  updated in web browser, using THttpServer monitoring capability like histogram objects.


## Tutorials
- Add the "Legacy" category collecting the old tutorials which do not represent any more best practices


## Class Reference Guide
- Images in tutorials can now be displayed à JavaScript thanks to the (js) option
  added next to the directive `\macro_image`
- As the tutorial `palettes.C` is often hit when searching the keyword `palette`
  in the reference guide, a direct link from this example to the full list of
  predefined palettes given in `TColor` has been added.
- Revisited the TSpectrum2 documentation. All the static images have been replaced
  by macros generating images at reference guide build time. These macros have
  been added in the tutorial section of the reference guide.
- The Reference Guide can now be accessed directly from the ROOT prompt thanks to
  a great extension (implemented by Desislava Kalaydjieva) of the `.help` command.
  For example to access the Reference Guide for `TTree` it is enough to type:
~~~ {.cpp}
   root[0] .help TTree
~~~
  To open the reference guide for a function/member:
~~~ {.cpp}
   root[0] .help TTree::Draw
~~~

## Build, Configuration and Testing Infrastructure

- Make MLP optional via the `-Dmlp={OFF,ON}` switch for CMake
- Make Spectrum optional via the `-Dspectrum={OFF,ON}` switch for CMake
- ROOT now fails to configure when any package is missing
  when `-Dfail-on-missing=ON` is passed to CMake
- The `-Dall=ON` now switches the default value of all optional packages to `ON`
- The options `astiff`, `cling`, `pch`, `thread`, and `explicitlink` have been
  removed and are now ignored. They either had no effect (their value was not
  being used in the build system), or could not be disabled (like `cling` and
  `explicitlink`).
- ROOT library targets now export which C++ standard they were built with via
  the target compile features `cxx_std_11`, `cxx_std_14`, and `cxx_std_17`.
- The file `RootNewMacros.cmake` has been renamed to `RootMacros.cmake`.
  Including the old file by name is deprecated and will generate a warning.
  Including `RootMacros.cmake` is not necessary, as now it is already included
  when calling `find_package(ROOT)`. If you still need to inherit ROOT's compile
  options, however, you may use `include(${ROOT_USE_FILE})` as before.
- ROOT's internal CMake modules (e.g. CheckCompiler.cmake, SetUpLinux.cmake, etc)
  are no longer installed with `make install`. Only the necessary files by
  dependent projects are installed by default now, and they are installed
  directly into the cmake/ directory, not cmake/modules/ as before.
- The macro `ROOT_GENERATE_DICTIONARY()` can now attach the generated source
  file directly to a library target by using the option `MODULE <library>`, where
  `<library>` is an existing library target. This allows the dictionary to inherit
  target properties such as compile options and include directories from the library
  target, even when they are added after the call to `ROOT_GENERATE_DICTIONARY()`.
- The macros `REFLEX_GENERATE_DICTIONARY()` and `ROOT_GENERATE_DICTIONARY()` can
  now have custom extra dependencies added with the options `DEPENDS` and
  `EXTRA_DEPENDENCIES`, respectively.

The following builtins have been updated:

- FFTW3 3.3.8
- GSL 2.5
- Intel TBB 2019 U8
- PCRE 8.43
- OpenSSL 1.0.2s
- Vdt 0.4.3
- VecCore 0.6.0
- XRootD 4.10.0
