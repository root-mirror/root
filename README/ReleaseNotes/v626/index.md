% ROOT Version 6.26 Release Notes
% 2021-03-03
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.26/00 is scheduled for release in May, 2021.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Bertrand Bellenot, CERN/SFT,\
 Josh Bendavid, CERN/CMS,\
 Jakob Blomer, CERN/SFT,\
 Rene Brun, CERN/SFT,\
 Will Buttinger, STFC/ATLAS,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Gerri Ganis, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Enrico Guiraud, CERN/SFT,\
 Jonas Hahnfeld, CERN/SFT,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Axel Naumann, CERN/SFT,\
 Max Orok, U Ottawa,\
 Alexander Penev, University of Plovdiv,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Jonas Rembser, CERN/SFT,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Oksana Shadura, UNL/CMS,\
 Matevz Tadel, UCSD/CMS,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,\
 Stefan Wunsch, CERN/SFT

## Deprecation, Removal, Backward Incompatibilities

- The "Virtual MonteCarlo" facility VMC (`montecarlo/vmc`) has been removed from ROOT. The development of this package has moved to a [separate project](https://github.com/vmc-project/). ROOT's copy of VMC was deprecated since v6.18.
- `TTreeProcessorMT::SetMaxTasksPerFilePerWorker` has been removed. `TTreeProcessorMT::SetTasksPerWorkerHint` is a superior alternative.
- `TTree::GetEntry()` and `TTree::GetEvent()` no longer have 0 as the default value for the first parameter `entry`. We are not aware of correct uses of this function without providing an entry number. If you have one, please simply pass `0` from now on.


## Core Libraries


## I/O Libraries


## TTree Libraries

## RDataFrame

### New features

- Add `Redefine` to the `RDataFrame` interface, which allows to overwrite the value of an existing column.
- Add `Describe` to the `RDataFrame` interface, which allows to get useful information, e.g. the columns and their types.
- Add `DescribeDataset` to the `RDataFrame` interface, which allows to get information about the dataset (subset of the output of Describe()).
- `Book` now suports just-in-time compilation, i.e. it can be called without passing the column types as template parameters (with some performance penalty, as usual).

## Histogram Libraries

- Implement the `SetStats` method for `TGraph` to turn ON or OFF the statistics box display
  for an individual `TGraph`.

## Math Libraries


## RooFit Libraries
### Creating RooFit datasets from RDataFrame
RooFit now contains two RDataFrame action helpers, `RooDataSetHelper` and `RooDataHistHelper`, which allow for creating RooFit datasets by booking an action:
```c++
  RooRealVar x("x", "x", -5.,   5.);
  RooRealVar y("y", "y", -50., 50.);
  auto myDataSet = rdataframe.Book<double, double>(
    RooDataSetHelper{"dataset",          // Name   (directly forwarded to RooDataSet::RooDataSet())
                    "Title of dataset",  // Title  (                   ~ " ~                      )
                    RooArgSet(x, y) },   // Variables to create in dataset
    {"x", "y"}                           // Column names from RDataFrame
  );
```
For more details, consult the tutorial [rf408_RDataFrameToRooFit](https://root.cern/doc/v626/rf408__RDataFrameToRooFit_8C.html).

### Changes in `RooAbsPdf::fitTo` behaviour for multi-range fits

The `RooAbsPdf::fitTo` and `RooAbsPdf::createNLL` functions accept a command argument to specify the fit range.
One can also fit in multiple ranges simultaneously.
The definition of such multi-range likelihoods for non-extended fits changes in this release.
Previously, the individual likelihoods were normalized separately in each range, which meant that the relative number of events in each sub-range was not used to estimate the PDF parameters.
From now on, the likelihoods are normalized by the sum of integrals in each range. This implies that the likelihood takes into account all inter-range and intra-range information.

### Deprecation of the `RooMinuit` class

The `RooMinuit` class was the old interface between RooFit and minuit. With ROOT version 5.24, a the more general `RooMinimizer` adapter was introduced, which became the default with ROOT 6.08. 

Before 6.26, it was possible to still use the `RooMinuit` by passing the `Minimizer("OldMinuit", "minimizer")` command argument to `RooAbsPdf::fitTo()`. This option is now removed.

## 2D Graphics Libraries

- Implement the option `X+` and `Y+` for reverse axis on TGraph.
- Offsets for axis titles with absolute-sized fonts (size%10 == 3) are now relative only to the font size (i.e. no longer relative to pad dimensions).

## 3D Graphics Libraries


## Geometry Libraries


## Database Libraries


## Networking Libraries


## GUI Libraries

## WebGUI Libraries

- provide `--web=server` mode, which only printout window URLs instead of starting real web browser.
  Dedicated for the case when ROOT should be running as server application, providing different RWebWindow instances for connection.


## Montecarlo Libraries


## PROOF Libraries


## Language Bindings


## JavaScript ROOT


## Jupyter lab

- Let use created notebooks with viewers like https://nbviewer.jupyter.org/
- Fix problem with using of local JSROOT version


## Tutorials


## Class Reference Guide

- Images for ROOT7 tutorials can be generated, in json format, using the directive using
  `\macro_image (json)` in the macro header.


## Build, Configuration and Testing Infrastructure
