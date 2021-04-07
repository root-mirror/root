% ROOT Version 6.26 Release Notes
% 2021-03-03
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.26/00 is scheduled for release in May, 2021.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this new version:

 Bertrand Bellenot, CERN/SFT,\
 Jakob Blomer, CERN/SFT,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Gerri Ganis, CERN/SFT,\
 Andrei Gheata, CERN/SFT,\
 Jonas Hahnfeld, CERN/SFT,\
 Sergey Linev, GSI,\
 Pere Mato, CERN/SFT,\
 Lorenzo Moneta, CERN/SFT,\
 Axel Naumann, CERN/SFT,\
 Danilo Piparo, CERN/SFT,\
 Fons Rademakers, CERN/SFT,\
 Jonas Rembser, CERN/SFT,\
 Enric Tejedor Saavedra, CERN/SFT,\
 Oksana Shadura, UNL/CMS,\
 Matevz Tadel, UCSD/CMS,\
 Vassil Vassilev, Princeton/CMS,\
 Wouter Verkerke, NIKHEF/Atlas,

## Deprecation and Removal

- `TTreeProcessorMT::SetMaxTasksPerFilePerWorker` has been removed. `TTreeProcessorMT::SetTasksPerWorkerHint` is a superior alternative.


## Core Libraries


## I/O Libraries


## TTree Libraries


## Histogram Libraries


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

## 2D Graphics Libraries


## 3D Graphics Libraries


## Geometry Libraries


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


