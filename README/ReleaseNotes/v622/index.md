% ROOT Version 6.22 Release Notes
% 2020-05-19
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.22/00 is scheduled for release in May, 2020.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Guilherme Amadio, CERN/SFT,\
 Bertrand Bellenot, CERN/SFT,\
 Jakob Blomer, CERN/SFT,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Enrico Guiraud, CERN/SFT,\
 Stephan Hageboeck, CERN/SFT,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Alja Mrak-Tadel, UCSD/CMS,\
 Jan Musinsky, SAS Kosice,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, CERN/SFT and UPV,\
 Danilo Piparo, CERN/SFT,\
 Timur Pocheptsoff, Qt Company,\
 Renato Quagliani, LPNHE, CNRS/IN2P3, Sorbonne Université,\
 Fons Rademakers, CERN/SFT,\
 Oksana Shadura, Nebraska,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Matevz Tadel, UCSD/CMS,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,\
 Stefan Wunsch, CERN/SFT

## Deprecation and Removal

- `ROOT::GetImplicitMTPoolSize` has been deprecated in favor of the newly added `ROOT::GetThreadPoolSize` and
  will be removed in v6.24.
- Manually setting `TThreadedObject::fgMaxSlots` is deprecated: TThreadedObject now increases the number of slots
  on-demand rather than running out and throwing an exception


## Core Libraries

- ROOT comes with C++ Modules enabled. More details about the technology found [here](../../README.CXXMODULES.md).
- The `ACLiC` can be configured to pass options to the `rootcling` invocation by enabling in the `.rootrc` the `ACLiC.ExtraRootclingFlags [-opts]` line.
- A call to `ROOT::EnableThreadSafety` is not required before using `TThreadExecutor` or `TTreeProcessorMT` anymore
- `TTreeProcessorMT` does not silently activate implicit multi-threading features anymore. An explicit call to
  `ROOT::EnableImplicitMT` is required instead
- `TTreeProcessorMT` now has a constructor argument to set the number of threads for its thread-pool


## I/O Libraries


## TTree Libraries

- A new status bit was added to `TTree`: `kEntriesReshuffled`, which indicates a `TTree` that is the output of the
  processing of another tree during which its entry order has been changed (this can happen, for instance, when
  processing a tree in a multi-thread application). To avoid silent entry number mismatches, trees with this bit set
  cannot add friend trees nor can be added as friends, unless the friend `TTree` has an appropriate `TTreeIndex`.


## Histogram Libraries


## Math Libraries


## RooFit Libraries

### RooWorkspace::Import() for Python
`RooWorkspace.import()` cannot be used in Python, since it is a reserved keyword. Users therefore had to resort
to
    getattr(workspace, 'import')(...)
Now,
    workspace.Import(...)
has been defined for the new PyROOT, which makes calling the function easier.

### Modernised category classes
RooFit's categories were modernised. Previously, the class RooCatType was used to store category states. It stores
two members, an integer for the category index, and up to 256 characters for a category name. Now, such states are
stored only using an integer, and category names can have arbitrary length. This will use 4 instead of 288 bytes
per category entry in a dataset, and make computations that rely on category states faster.

The interface to define or manipulate category states was also updated. Since categories are mappings from state names
to state index, this is now reflected in the interface. Among others, this is now possible:
| ROOT 6.22                                      | Before (still supported)                                       |
|------------------------------------------------|----------------------------------------------------------------|
| `RooCategory cat("cat", "Lepton flavour");`    | `RooCategory cat("cat", "Lepton flavour");`                    |
| `cat["electron"] = 1;`                         | `cat.defineType("electron", 1);`                               |
| `cat["muon"] = 2;`                             | `cat.defineType("muon", 2);`                                   |

See also [Category reference guide](https://root.cern.ch/doc/master/classRooCategory.html).

### Type-safe proxies for RooFit objects
RooFit's proxy classes have been modernised. The class `RooTemplateProxy` allows for access to other RooFit objects
similarly to a smart pointer. In older versions of RooFit, the objects held by *e.g.* `RooRealProxy` had to be
accessed like this:
    RooAbsArg* absArg = realProxy.absArg();
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(absArg);
    assert(pdf); // This *should* work, but the proxy doesn't have a way to check
    pdf->fitTo(...);
That is, a `RooRealProxy` stores a pointer to a RooAbsArg, and this pointer has to be cast. There was no type
safety, *i.e.*, any object deriving from RooAbsArg could be stored in that proxy, and the user had to take care
of ensuring that types are correct.
Now, if one uses
    RooTemplateProxy<RooAbsPdf> pdfProxy;
instead of
    RooRealProxy realProxy;
in RooFit classes, the above code can be simplified to
    pdfProxy->fitTo(...);

Check the [doxygen reference guide](https://root.cern.ch/doc/master/classRooTemplateProxy.html) for `RooTemplateProxy` for
more information on how to modernise old code.

### HistFactory

#### Switch default statistical MC errors to Poisson
When defining HistFactory samples with statistical errors from C++, e.g.
    Sample background1( "background1", "background1", InputFile );
    background1.ActivateStatError();
statistical MC errors now have Poisson instead of Gaussian constraints. This better reflects the uncertainty of the MC simulations.
This can be reverted as follows:
    // C++:
    Channel chan("channel1");
    chan.SetStatErrorConfig( 0.05, "Gauss" );
    // Within <Channel ... > XML:
    <StatErrorConfig RelErrorThreshold="0.05" ConstraintType="Gauss" />

#### Less verbose HistFactory
HistFactory was very verbose, writing to the terminal with lots of `cout`. Now, many HistFactory messages are going
into RooFit's message stream number 2. The verbosity can therefore be adjusted using
    RooMsgService::instance().getStream(2).minLevel = RooFit::PROGRESS;

`hist2workspace` is also much less verbose. The verbosity can be restored with `hist2workspace -v` or `-vv`.

## 2D Graphics Libraries

  - Universal time (correct time zone and daylight saving time) in PDF file. Implemented by
    Jan Musinsky.
  - The crosshair type cursor type did not work on MacOS Catalina. This has been fixed by
    Timur Pocheptsoff.
  - Take into account the Z errors when defining the frame to paint a TGraph2DErrors.
  - Implement the of "F" in `TPad::RedrawAxis` to allow the plot's frame redrawing when
    erased.
  - Implement `TCanvas::SetRealAspectRatio` to resize a canvas so that the plot inside is
    shown in real aspect.
  - New graphics style "BELLE2" from Martin Ritter.

## 3D Graphics Libraries


## Geometry Libraries

### Geometry drawing in web browser

When ROOT compiled with -Droot7=ON flag, one can enable geometry drawing in web browser.
Just apply --web option when starting root like: `root --web tutorials/geom/rootgeom.C`
Not all features of TGeoPainter are supported - only plain drawing of selected TGeoVolume


## Database Libraries


## Networking Libraries


## GUI Libraries


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings


## JavaScript ROOT


## Tutorials


## Class Reference Guide


## Build, Configuration and Testing Infrastructure

By default, ROOT now falls back to the built-in version of xrootd if it can't find it in the system.
This means that passing `-Dbuiltin_xrootd=ON` is not necessary anymore to build ROOT with xrootd support.
Note that built-in xrootd requires a working network connection.

### Experimental address sanitizer build configuration
Added a build flag `asan` that switches on address sanitizer. It's experimental, so expect problems. For example, when building with gcc,
manipulations in global variables in llvm will abort the build. Such checks can be disabled using environment variables. Check the address
sanitizer documentation or the link below for details. In clang, which allows to blacklist functions, the build will continue.

See [core/sanitizer](https://github.com/root-project/root/tree/master/core/sanitizer) for information.


### Optimization of ROOT header files

Many (but intentionally not all) unused includes were removed from ROOT header files. For instance, `#include "TObjString.h"` and
`#include "ThreadLocalStorage.h"` were removed from `TClass.h`. Or `#include "TDatime.h"` was removed from
`TDirectory.h` header file . Or `#include "TDatime.h"` was removed from `TFile.h`.
This change may cause errors during compilation of ROOT-based code. To fix it, provide missing the includes
where they are really required.
This improves compile times and reduces code inter-dependency; see https://github.com/include-what-you-use/include-what-you-use/blob/master/docs/WhyIWYU.md for a good overview of the motivation.

Even more includes will be "hidden" when ROOT configured with `-Ddev=ON` build option.
In that case ROOT uses `#ifdef R__LESS_INCLUDES` to replace unused includes by class forward declarations.
Such `dev` builds can be used to verify that ROOT-based code really includes all necessary ROOT headers.


## RDataFrame

- Starting from this version, when `RSnapshotOptions.fMode` is `"UPDATE"` (i.e. the output file is opened in "UPDATE"
  mode), Snapshot will refuse to write out a TTree if one with the same name is already present in the output file.
  Users can set the new flag `RSnapshotOption::fOverwriteIfExists` to `true` to force the deletion of the TTree that is
  already present and the writing of a new TTree with the same name. See
  [ROOT-10573](https://sft.its.cern.ch/jira/browse/ROOT-10573) for more details.
- RDataFrame changed its error handling strategy in case of unreadable input files. Instead of simply logging an error
  and skipping the file, it now throws an exception if any of the input files is unreadable (this could also happen in
  the middle of an event loop). See [ROOT-10549](https://sft.its.cern.ch/jira/browse/ROOT-10549) for more details.
- New analysis examples based on the recent ATLAS Open Data release ([`Higgs to two photons`](https://root.cern/doc/master/df104__HiggsToTwoPhotons_8py.html), [`W boson analysis`](https://root.cern/doc/master/df105__WBosonAnalysis_8py.html), [`Higgs to four leptons`](https://root.cern/doc/master/df106__HiggsToFourLeptons_8py.html))
- An exception is now thrown in case the size of ROOT's thread-pool changes between RDataFrame construction time and the time the event loop begins.
- Just-in-time compilation of large portions of the computation graph has been optimized, and it is now much faster. Please report any regressions you might encounter on [our issue tracker](https://sft.its.cern.ch/jira/projects/ROOT).
- `MakeRootDataFrame` is now a safe way to construct RDFs. It used to return RDFs with more limited functionality.


## PyROOT

- Introduce the `ROOT.Numba.Declare` decorator which provides a simple way to call Python callables from C++. The Python callables are
  just-in-time compiled with [numba](http://numba.pydata.org/), which ensures a runtime performance similar to a C++ implementation.
  The feature is targeted to improve the performance of Python based analyses, e.g., allows seamless integration into `RDataFrame` workflows.
  See the tutorial [`pyroot004_NumbaDeclare.py`](https://root.cern/doc/master/pyroot004__NumbaDeclare_8py.html) for further information.
