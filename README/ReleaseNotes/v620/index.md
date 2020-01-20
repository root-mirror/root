% ROOT Version 6.20 Release Notes
% 2020-01-10
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.20/00 is scheduled for release in January 2020.

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
 Massimiliano Galli, CERN/SFT and Unibo,\
 Gerri Ganis, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Hadrien Grasland, CNRS,\
 Enrico Guiraud, CERN/SFT,\
 Stephan Hageboeck, CERN/SFT,\
 Desislava Kalaydjieva, CERN/SFT,\
 Jan Knedlik, GSI,\
 Philip Leindecker, CERN/SFT,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Emmanouil Michalainas, AUTh,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak-Tadel, UCSD/CMS,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, Bicocca/SFT,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Otto Schaile, Uni-Muenchen,\
 Henry Schreiner, Princeton,\
 Oksana Shadura, Nebraska,\
 Simon Spies, GSI,\
 Matevz Tadel, UCSD/CMS,\
 Yuka Takahashi, Princeton and CERN/SFT,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/ATLAS,\
 Stefan Wunsch, CERN/SFT,\
 Zhe Zhang, Nebraska

## General

### Splash screen

The venerable splash screen is now disabled by default to make ROOT's startup
faster. Many users already use `root -l` to start ROOT, but this also hides the
useful text banner with version information along with the splash screen. With
this new default, starting up ROOT as just `root` will show only the text banner
instead of the splash screen. The splash screen can still be seen with `root -a`
or in `TBrowser` by opening `Browser Help → About ROOT`.

## Deprecation and Removal

 * rootcling flags `-cint`, `-gccxml`, `-p`, `-r` and `-c` have no effect
   and will be removed in a future release. Please remove them from the rootcling invocations.
 * rootcling legacy cint flags `+P`, `+V` and `+STUB` have no effect and will be
   removed in a future release. Please remove them from the rootcling invocations.
 * genreflex flag `--deep` has no effect and will be removed in a future release. Please remove it
   from the genreflex invocation.
 * rootcling warns if it sees and unrecognized flag (usually coming from the
   CXXFLAGS of the build system). Please remove them from the invocation because
   the warning will become a hard error in the next releases.
 * The empty headers `Gtypes.h` and `Htypes.h` are deprecated. Please include
   `Rtypes.h`

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

       TFile *f = TFile::Open("name.root?reproducible","RECREATE","File title");
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

### RDataFrame

- Improved CSV data source
  * Support for 1.0e+10 syntax for type inference of double columns
  * Empty lines are now skipped


## Histogram Libraries

* Allow reading v5 TF1 that were stored memberwise in a TClonesArray.
* Make RHist bin iteration order consistent with that of THx.
* Add RCoordArray constructor from std::array.
* Remove thread-unsafe accessors of RHistConcurrentFill.
* Fix Flush logic of buffered RHist wrappers.
* New class TGraphMultiErrors: A TGraph with asymmetric error bars and multiple y error
  dimensions (author: Simon Spies).

## Math Libraries
* [ROOT::Math::KahanSum](https://root.cern/doc/master/classROOT_1_1Math_1_1KahanSum.html) can use SIMD instructions for faster accumulation

## RooFit Libraries
* **Documentation** Many improvements to the doxygen documentation of RooFit classes.

* **Automatic legends for RooPlot** RooPlot now supports [BuildLegend](https://root.cern.ch/doc/master/classRooPlot.html#acd16cf22aca843f08ef405c93753c512) as a good starting point for a legend.

* **Short prefits** In unbinned fits that take long, good starting values for parameters can be found
by running a short prefit before the final fit. Passing `PrefitDataFraction(0.1)` to [`fitTo()`](https://root.cern.ch/doc/master/classRooAbsPdf.html#a8f802a3a93467d5b7b089e3ccaec0fa8)
will *e.g.* run a prefit on 1% of the data.

* **Iterating over categories** Category classes deriving from [RooAbsCategory](https://root.cern.ch/doc/master/classRooAbsCategory.html), *e.g.* [RooCategory](https://root.cern.ch/doc/master/classRooCategory.html)
now support "natural" iterating using range-based for loops in C++ or Python loops:

      import ROOT
      cat = ROOT.RooCategory("cat", "cat")
      cat.defineType("1Lep", 1)
      cat.defineType("2Lep", 2)
      for state in cat:
        print(state.getVal(), state.GetName())

      (1, '1Lep')
      (2, '2Lep')
* **Asymptotically correct parameter uncertainties**
Added computation of asymptotically correct parameter uncertainties in likelihood fits with event weights. See [arXiv 1911.01303](https://arxiv.org/abs/1911.01303) and [rf611_weightedfits.C](https://root.cern/doc/master/rf611__weightedfits_8C.html).

* **Barlow-Beeston tutorial**
The tutorial [rf709_BarlowBeeston.C](https://root.cern/doc/master/rf709__BarlowBeeston_8C.html) has
been added to RooFit to demonstrate how to incorporate Monte Carlo statistics as systematic
uncertainties into a fit. It demonstrates both the "full" and the "light" method.

* **RooFitMore for GSL** All parts of RooFit that use the GSL (PDFs and integrators) have been moved into the library
`RooFitMore`. It gets enabled automatically with the `MathMore` library (`-Dmathmore=ON`, default).
Note that `-lRooFitMore` might now be required when linking against RooFit.

### Fast function evaluation and vectorisation
A `BatchMode` for faster unbinned fits has been added. By loading data more efficiently,
**unbinned likelihood computations can run about 3x faster** if the PDFs support it.
To enable it, use
```
pdf.fitTo(data, RooFit::BatchMode());
```
Most unbinned PDFs that are shipped with RooFit have been updated to support this mode.

In addition, if ROOT is compiled for a specific architecture, SIMD instructions can be used in PDF
computations. This requires ROOT to be compiled with `-march=native` or *e.g.* `-mavx2` if the hardware
supports it. For maximal performance, ROOT should also be configured with `-Dvdt=ON`.
[VDT](https://github.com/dpiparo/vdt) is a library of fast math functions, which will automatically
be used in RooFit when available.
Depending on the compiler, on the instruction set supported by the CPU and on what kind of PDFs are used,
**PDF evaluations will speed up 5x to 16x**.
For details see [CHEP 2019](https://indico.cern.ch/event/773049/contributions/3476060/).

### New RooFormulaVar / RooGenericPdf
RooFormula has been updated to use ROOT's [TFormula](https://root.cern.ch/doc/master/classTFormula.html).
This means that expressions passed to RooFormulaVar / RooGenericPdf are compiled with `-O2` optimisation
level before they are evaluated. For complicated functions, this might improve the speed of the
computation. Further, functions that are known to the interpreter can be used in the expression passed
to a PDF. The following macro *e.g.* prints the expected result 5.4:
```
  double func(double x, double a) {
    return a*x*x + 1.;
  }

  void testRooFormulaWithClingFunction() {
    RooRealVar x("x", "x", 2, -10, 10);
    RooRealVar a("a", "a", 1.1, -10, 10);

    RooGenericPdf pdf("pdfWithExternalFunc", "func(x, a)", {a, x});
    std::cout << pdf.getVal() << std::endl;
  }
```

#### New PDFs added to RooFit
* [RooHypatia2](https://root.cern.ch/doc/master/classRooHypatia2.html), featuring a hyperbolic, Crystal-ball-like core and two adjustable tails.
* [RooWrapperPdf](https://root.cern.ch/doc/master/classRooWrapperPdf.html), a class that can convert any function (= not normalised) into a PDF by integrating and normalising it.


### RooStats / HistFactory
* To facilitate plotting of pre-fit model uncertainties, gamma parameters now have reasonable pre-fit uncertainties.
* RooStats global config switches are now all accessible in one place,
[`RooStats::GetGlobalRooStatsConfig()`](https://root.cern.ch/doc/master/namespaceRooStats.html#a827f04b74fab219d613f178fa24d0bc9).
* Tools like ToyMCSampler and HypoTestInverter have been stabilised and better validate their inputs to
prevent infinite loops or crashes.


## 2D Graphics Libraries

 * Provide support of NDC coordinates for `TArrow`.
 * Fix interactive movement of `TLine/TArrow` objects when NDC coordinates are used.
 * Provide `TGraph::MovePoints()` method.
 * New options `RX`and `RY` for TMultiGraph in order to draw reverse axis along X and Y.
 * Combined with the option "Z" the option "CJUST" allows to draw the color palette
   with axis labels justified on the color boundaries (implemented by Otto Schaile).
 * The `TCanvas` Event Status Bar now displays the date and time when the mouse cursor
   is moved over a time axis (implemented by Otto Schaile).
 * Negative values were not painted with option "TEXT" for `TH2Poly`.
 * The Z title was not properly set in `TEfficiency`.
 * Implement `TImage:ls()`
 * `TGraph2D`: X Y and Z titles were not saved by SavePrimitive.
 * In some cases vertical hatches drawn by `TPad::PaintHatches` disappeared.
 * Save the X and Y axis limits in `TMultiGraph::SavePrimitive`.
 * Some markers did not match between PDF and PNG.
 * The MaxDigits attribute was not imported from TAxis to TGaxis.


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
- JSROOT graphics are now supported in the JupyterLab interface. They are activated in the same way as in
  the classic Jupyter, i.e. by typing at the beginning of a notebook cell:
~~~ {.python}
%jsroot on
~~~

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

      root[0] .help TTree

  To open the reference guide for a function/member:

      root[0] .help TTree::Draw

- Change the layout of the ROOT reference.

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


## PyROOT

### Current PyROOT

- Several changes for forward compatibility with experimental PyROOT and new Cppyy have been added:
  * Template instantiation can be done with square brackets. The parenthesis syntax is deprecated.
~~~ {.python}
my_templated_function['int','double']()  # new syntax
my_templated_function('int','double')()  # old sytax, throws a deprecation warning
~~~
  * When converting `None` to a null pointer in C++, a deprecation warning is issued.
  * When using `buffer.SetSize` a deprecation warning is issued, the forward compatible alternative is `buffer.reshape`.
  * When using `ROOT.Long` or `ROOT.Double`, a deprecation warning is issued in favour of their equivalent `ctypes` types
(`c_long`, `c_int`, `c_double`)
  * Added the forward compatible names `as_cobject` and `bind_object` for `AsCObject` and `BindObject`, respectively.
  * nullptr is also accessible as `cppyy.nullptr`, not only as `cppyy.gbl.nullptr`.
  * Pythonization functions (e.g. `add_pythonization`) are accessible via `cppyy.py`.
  * Some attributes of Python proxies have been added with the name they have in the new Cppyy
(`__creates__`, `__mempolicy__`, `__release_gil__` for function proxies, `__smartptr__` for object proxies).
- The support for enums (both scoped and non-scoped) was improved. Now, when creating an enum from Python,
its underlying type is obtained by PyROOT.
- Added support for non-ASCII Python strings (e.g. UTF-8) to `std::string` and `C string`.
- Added converters from `bytes` to `std::string` and to C string.
- Added compatibility of STL iterators with GCC9.
- Added support for templated methods with reference parameters.
- Introduced two teardown modes: soft and hard. The soft one only clears the proxied objects, while the hard one also
shuts down the interpreter.

### Experimental PyROOT

- MultiPython: build and install PyROOT with multiple Python versions.
  - Build:
  ~~~ {.bash}
  cmake -DPYTHON_EXECUTABLE=/path/to/first/Python/installation /path/to/ROOT/source
  cmake --build .
  cmake -DPYTHON_EXECUTABLE=/path/to/second/Python/installation /path/to/ROOT/source
  cmake --build .
  ( ... )
  ~~~
  - Source a specific built version (last one picked as default):
  ~~~ {.bash}
  ROOT_PYTHON_VERSION=X.Y source /path/to/bin/thisroot.sh
  ~~~
  - PyROOT installation directory can be customized:
  ~~~ {.bash}
  cmake -DCMAKE_INSTALL_PYROOTDIR=/path/to/PyROOT/install/dir /path/to/ROOT/source
  cmake --build .
  make install
  ~~~
- Updated cppyy packages to the following versions:
  * cppyy: cppyy-1.5.3
  * cppyy_backend: clingwrapper-1.10.3
  * CPyCppyy: CPyCppyy-1.9.3
- Introduced two teardown modes: soft and hard. The soft one only clears the proxied objects, while the hard one also
shuts down the interpreter.
- A few changes for backward compatibility with current PyROOT have been added:
  * Added `MakeNullPointer(klass)` as `bind_object(0,klass)`
  * Provided `BindObject` and `AsCObject`
- `ROOT.Long` and `ROOT.Double` are no longer in the API of PyROOT, since new Cppyy requires to use their equivalent
`ctypes` types (`c_long`, `c_int`, `c_double`).
- Added TPython.
- Added support for `from ROOT import *` in Python2 (only).
