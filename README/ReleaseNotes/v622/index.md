% ROOT Version 6.22 Release Notes
% 2020-06-14
<a name="TopOfPage"></a>

## Introduction

ROOT version 6.22/00 was released on June 14, 2020.

For more information, see:

[http://root.cern](http://root.cern)

The following people have contributed to this version:

 Guilherme Amadio, CERN/SFT,\
 Bertrand Bellenot, CERN/SFT,\
 Jakob Blomer, CERN/SFT,\
 Rene Brun, CERN/SFT,\
 Philippe Canal, FNAL,\
 Olivier Couet, CERN/SFT,\
 Surya Dwivedi, GSOC/SFT, \
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
 Jan Musinsky, SAS Kosice,\
 Axel Naumann, CERN/SFT,\
 Vincenzo Eduardo Padulano, Bicocca/SFT,\
 Danilo Piparo, CERN/SFT,\
 Timur Pocheptsoff, Qt Company,\
 Renato Quagliani, LPNHE, CNRS/IN2P3, Sorbonne Université,\
 Fons Rademakers, CERN/SFT,\
 Oksana Shadura, UNL/CMS,\
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


## TTree Libraries

- A new status bit was added to `TTree`: `kEntriesReshuffled`, which indicates a `TTree` that is the output of the
  processing of another tree during which its entry order has been changed (this can happen, for instance, when
  processing a tree in a multi-thread application). To avoid silent entry number mismatches, trees with this bit set
  cannot add friend trees nor can be added as friends, unless the friend `TTree` has an appropriate `TTreeIndex`.
- `TTree::GetBranch` has been updated to execute its search trough the set of branches breadth first instead of depth first.  In some cases this changes the branch that is returned.  For example with:

```
    struct FloatInt {
       float f;
       int x;
    };
```
and

```q
   int x = 1;
   FloatInt s{2.1, 3};
   TTree t("t", "t");
   t.Branch("i", &i); // Create a top level branch named "i" and a sub-branch name "x"
   t.Branch("x", &x);
```
the old version of `t.GetBranch("x")` was return the `i.x` sub-branch while the new version return the `x` top level branch.

- `TTree::GetBranch` was also upgraded to properly handle being give a full branch name.  In prior version, with the example above, `GetBranch("i.x")` would return nullptr.  With the new version it returns the address of the branch `i.x`.

- `TChain` now properly carries from `TTree` to `TTree`, each TBranch `MakeClass` (also known as `DecomposedObj[ect]` mode, if it was set (implicitly or explicitly).

- `TTree::SetBranchAddress` and `TChain::SetBranchAddress` now automatically detect when a branch needs to be switched to MakeClass` (also known as `DecomposedObj[ect]` mode).

-  `TChain` now always checks the branch address; even in the case the address is set before any TTree is loaded. Specifically in addition to checking the address type in the case:

```
chain->LoadTree(0);
chain->SetBranchAdress(branch_name, &address);
```

now `TChain` will also check the addresses in the case:

```
chain->SetBranchAdress(branch_name, &address);
chain->LoadTree(0);
```

- `TTree` and `TChain` no longet override user provided sub-branch addresses.
Prior to this release doing:

```
// Not setting the top level branch address.
chain->SetBranchAdress(sub_branch_name, &address);
chain->GetEntry(0);
```

Resulted in the address set to be forgotten.  Note that a work-around was:

```
// Not setting the top level branch address.
chain->GetEntry(0);
chain->SetBranchAdress(sub_branch_name, &address);
```

But also the address needed to (in most cases) also be set again after each new tree was loaded.

Note that, the following:

```
chain->SetBranchAdress(sub_branch_name, &address);
chain->SetBranchAdress(top_level_branch_name, &other_address);
chain->GetEntry(0);
```

will result (as one would expect) with the first SetBranchAddress being ignored/over-ridden.


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


## Histogram Libraries

- When fitting an histogram using the fit option `W` (set all bin errors to be equal to 1),  renormalize the obtained fit  parameter errors using the obtained chi2 fit value. This is the same procedure performed already when fitting a `TGraph`.
- Applying some fixes when performing histogram fits using bin integral in multi-thread mode. 

## Math Libraries

### Minuit2

- Fix the case when a new better minimum is found while running MINOS
- add some fixed in handling the message printing
- use std::numeric_limits<double> for the default value of the internal precision

### Matrix

- Improve `TDecompQRH` , QR decomposition in an orthogonal matrix Q. Add a function `TDecompQR::GetOrthogonalMatrix()` returning the
Q matrix found by the decomposition. 


## RooFit Libraries

### RooWorkspace::Import() for Python
`RooWorkspace.import()` cannot be used in Python, since it is a reserved keyword. Users therefore had to resort
to
```
getattr(workspace, 'import')(...)
```
Now,
```
workspace.Import(...)
```
has been defined for the new PyROOT, which makes calling the function easier.


### Modernised category classes
RooFit's categories were modernised. Previously, the class RooCatType was used to store category states. It stores
two members, an integer for the category index, and up to 256 characters for a category name. Now, such states are
stored only using an integer, and category names can have arbitrary length. This will use 4 instead of 288 bytes
per category entry in a dataset, and paves the way for faster computations that rely on category states.

The interface to define or manipulate category states was also updated. Since categories are mappings from state names
to state index, this is now reflected in the interface. Among others, this is now possible:
| ROOT 6.22                                      | Before (still supported)                                       |
|------------------------------------------------|----------------------------------------------------------------|
| `RooCategory cat("cat", "Lepton flavour");`    | `RooCategory cat("cat", "Lepton flavour");`                    |
| `cat["electron"] = 1;`                         | `cat.defineType("electron", 1);`                               |
| `cat["muon"] = 2;`                             | `cat.defineType("muon", 2);`                                   |

See also the [Category reference guide](https://root.cern.ch/doc/v622/classRooCategory.html) or different
[RooFit tutorials](https://root.cern.ch/doc/v622/group__tutorial__roofit.html),
specifically [rf404_categories](https://root.cern.ch/doc/v622/rf404__categories_8C.html).


### Type-safe proxies for RooFit objects
RooFit's proxy classes have been modernised. The class [RooTemplateProxy](https://root.cern.ch/doc/v622/classRooTemplateProxy.html)
allows for access to other RooFit objects
similarly to a smart pointer. In older versions of RooFit, the objects held by *e.g.* `RooRealProxy` had to be
accessed like this:
```
  RooAbsArg* absArg = realProxy.absArg();
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(absArg);
  assert(pdf); // This *should* work, but the proxy doesn't have a way to check
  pdf->fitTo(...);
```
That is, because the old `RooRealProxy` stored a pointer to a `RooAbsArg`, and this pointer had to be cast. There was no type
safety, *i.e.*, any object deriving from RooAbsArg could be stored in that proxy, and the user had to take care
of ensuring that types are correct.
Now, if one uses
```
  RooTemplateProxy<RooAbsPdf> pdfProxy; // instead of "RooRealProxy realProxy;"
```
the above code can be simplified to
```
  pdfProxy->fitTo(...);
```

Check the [doxygen reference guide](https://root.cern.ch/doc/v622/classRooTemplateProxy.html) for `RooTemplateProxy` for
more information on how to modernise old code.

### Documentation updates
- The documentation is now structured better. The functional parts like roofitcore, roofit, roostats, histfactory
can now be found in [doxygen](https://root.cern.ch/doc/v622/group__Roofitmain.html) with a similar structure as
the [directory structure in git](https://github.com/root-project/root/tree/master/roofit).
- Especially the hard-to-find command arguments that can be passed to various functions in RooFit now have
[their own doxygen page](https://root.cern.ch/doc/v622/group__CmdArgs.html).
- Clarified and extended documentation for plotting PDFs:
  - Since plotting PDFs can be tricky when special normalisation is required (e.g. because of blinding or fitting in side bands),
    the documentation of [RooAbsPdf::plotOn()](https://root.cern.ch/doc/v622/classRooAbsPdf.html#a7f01ccfb4f3bc3ad3b2bc1e3f30af535)
    and the arguments it accepts was updated and clarified.
  - The new tutorial [rf212_plottingInRanges_blinding](https://root.cern.ch/doc/v622/rf212__plottingInRanges__blinding_8C.html)
    shows how to normalise a PDF if it has been fit only in a part of the plotting range.
- Started grouping of interfaces that are public but not necessarily useful for "normal" in doxygen.
Examples are [RooCategory](https://root.cern.ch/doc/v622/classRooCategory.html) or the very broad public interface of
[RooAbsArg](https://root.cern.ch/doc/v622/classRooAbsArg.html),
where the client-server interface that's only relevant for RooFit or the legacy interfaces are now split off.
- Completely overhaul the documentation of [RooWorkspace's factory syntax](https://root.cern.ch/doc/v622/classRooWorkspace.html#a0ddded1d65f5c6c4732a7a3daa8d16b0).
- Make [SPlot](https://root.cern.ch/doc/v622/classRooStats_1_1SPlot.html) documentation understandable.
- Clarify different constructors of [RooRealVars](https://root.cern.ch/doc/v622/classRooRealVar.html#aa1483bd1c3e7b5792765043f6ac0a50f),
  and how to make them constant or initialise them.
- Add a tutorial about plotting with limited ranges [rf212_plottingInRanges_blinding](https://root.cern.ch/doc/v622/rf212__plottingInRanges__blinding_8C.html)
- Improve the tutorial for weighted fits, [rf611_weightedfits](https://root.cern.ch/doc/v622/rf611__weightedfits_8C.html).

### Etc.
- Warn and fail better. Convert things that silently give the wrong results into things that fail noisily.
  - When RooDataSet is read from a file, and there's a read error, the function now returns `nullptr` instead of an incomplete dataset.
  - Better checking of invalid formulae in RooGenericPdf / RooFormulaVar. Using e.g. `RooFormulaVar(... "x+y", x);` silently
    assumed that `y` is zero and constant.
  - Ensure that the right number of coefficients is given for [RooAddPdf](https://root.cern.ch/doc/v622/classRooAddPdf.html) (depends on
    recursive/non-recursive fractions), and complain if that's not the case.
- Add automatic conversions when reading values for variables from Trees. New are e.g. `ULong64_t, Long64_t, Short_t, UShort_t`.
- [RooResolutionModel](https://root.cern.ch/doc/v622/classRooResolutionModel.html) now also takes [RooLinearVar](https://root.cern.ch/doc/v622/classRooLinearVar.html) as observables. With this, it is possible to easily scale and shift the observable used in the convolutions.
- Add a [helper](https://root.cern.ch/doc/v622/classRooHelpers_1_1LocalChangeMsgLevel.html) to locally/temporarily change RooFit's message
level. This also includes an [object that can hijack](https://root.cern.ch/doc/v622/classRooHelpers_1_1HijackMessageStream.html)
all messages sent by a specific object or with specific topics.

### HistFactory

#### Switch default statistical MC errors to Poisson
When defining HistFactory samples with statistical errors from C++, e.g.
```
  Sample background1( "background1", "background1", InputFile );
  background1.ActivateStatError();
```
statistical MC errors now default to Poisson instead of Gaussian constraints. This better reflects the uncertainty of the MC simulations,
and now actually implements what is promised in the [HistFactory paper](https://cds.cern.ch/record/1456844). For large number of entries
in a bin, this doesn't make a difference, but for bins with low statistics, it better reflects the "counting" nature of bins.
This can be reverted as follows:
```
  // C++:
  Channel chan("channel1");
  chan.SetStatErrorConfig( 0.05, "Gauss" );
```
```
  // Within a <Channel ... > XML:
  <StatErrorConfig RelErrorThreshold="0.05" ConstraintType="Gauss" />
```

#### Less verbose HistFactory
HistFactory was very verbose, writing to the terminal with lots of `cout`. Now, many HistFactory messages are going
into [RooFit's message stream](https://root.cern.ch/doc/v622/classRooMsgService.html) number 2.
The verbosity can therefore be adjusted using
```
  RooMsgService::instance().getStream(2).minLevel = RooFit::PROGRESS;
```

`hist2workspace` is also much less verbose. The verbosity can be restored with `hist2workspace -v` for intermediate verbosity
or `-vv` for what it printed previously.

## TMVA

### Deep Learning module

- Extend Deep Learning module by adding support for *LSTM* and *GRU* Recurrent layers.
- Add implementation for all recurrent layers (*LSTM*, *GRU*, and simple RNN) for GPU using the CUDA *cuDNN* library
- Add support for Batch Normalization implemented in both CPU and GPU
- Deprecate `MethodDNN` , defined by the `TMVA::kDNN` enumeration for building neural networks with fully connected dense layers.  It is replaced by `MethodDL` (`TMVA::kDL`)
- Deprecate also the option `Architecture=Standard` for `MethodDL`. Have only `CPU` or `GPU` architecture options. When a BLAS implementation is not found use the CPU architecture with matrix
operations provided by the ROOT TMatrix classes.
- Add new tutorials showing usage new Deep Learning functionalities
- All these new updates are documented in the new version of the **TMVA Users Guide**, available in github
[here](https://github.com/root-project/root/blob/master/documentation/tmva/UsersGuide/TMVAUsersGuide.pdf).

### PyMVA

- Add support for Tensorflow version 2 as backend for Keras in `MethodPyKeras`. Note that it requires still an independent Keras installation with a Keras version >= 2.3
- Add also some fixes needed when using Python  version 3. 

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

## Geometry Libraries

### Geometry drawing in web browser

When ROOT compiled with -Droot7=ON flag, one can enable geometry drawing in web browser.
Just apply --web option when starting root like: `root --web tutorials/geom/rootgeom.C`
Not all features of TGeoPainter are supported - only plain drawing of selected TGeoVolume


## Language Bindings

### PyROOT

ROOT 6.22 makes the new (experimental) PyROOT its default. This new PyROOT is designed on top of the new cppyy, which
provides more and better support for modern C++. The documentation for cppyy and the new features it provides can be
found [here](https://cppyy.readthedocs.io).

For what concerns new additions to PyROOT, this is the list:
- The `ROOT.Numba.Declare` decorator provides a simple way to call Python callables from C++. The Python callables are
  just-in-time compiled with [numba](http://numba.pydata.org/), which ensures a runtime performance similar to a C++ implementation.
  The feature is targeted to improve the performance of Python based analyses, e.g., allows seamless integration into `RDataFrame` workflows.
  See the tutorial [`pyroot004_NumbaDeclare.py`](https://root.cern/doc/master/pyroot004__NumbaDeclare_8py.html) for further information.
- Multi-Python builds are now supported (see more info in the `Build, Configuration and Testing Infrastructure` section).
- The two PyROOT teardown modes (soft and hard) are now fully functional. In the soft mode, at teardown of the Python interpreter,
  all the proxied C++ objects will be cleaned; in the hard mode, in addition to the proxy cleaning, gROOT.EndOfProcessCleanups() will be invoked
  and thus the C++ interpreter will be shut down. The hard mode is the default, but if the application just runs a PyROOT script as part of
  a longer process and it would like ROOT to still be usable after the script finishes, the soft mode can be activated by adding this to the
  PyROOT script:
~~~ {.python}
import ROOT
ROOT.PyConfig.ShutDown = False
~~~

On the other hand, there are some backwards-incompatible changes of the new PyROOT with respect to the new one, listed next:
- Instantiation of function templates must be done using square brackets instead of parentheses. For example, if we consider the following
  code snippet:
~~~ {.python}
> import ROOT

> ROOT.gInterpreter.Declare("""
template<typename T> T foo(T arg) { return arg; }
""")

> ROOT.foo['int']  # instantiation
 cppyy template proxy (internal)

> ROOT.foo('int')  # call (auto-instantiation from arg type)
'int'
~~~

Note that the above does not affect class templates, which can be instantiated either with parenthesis or square brackets:
~~~ {.python}
> ROOT.std.vector['int'] # instantiation
<class cppyy.gbl.std.vector<int> at 0x5528378>

> ROOT.std.vector('int') # instantiation
<class cppyy.gbl.std.vector<int> at 0x5528378>
~~~

- Overload resolution in new cppyy has been significantly rewritten, which sometimes can lead to a different overload choice
(but still a compatible one!). For example, for the following overloads of `std::string`:
~~~ {.cpp}
string (const char* s, size_t n)                           (1)
string (const string& str, size_t pos, size_t len = npos)  (2)
~~~
when invoking `ROOT.std.string(s, len(s))`, where `s` is a Python string, the new PyROOT will pick (2) whereas the old
would pick (1).

- The conversion between `None` and C++ pointer types is not allowed anymore. Instead, `ROOT.nullptr` should be used:
~~~ {.python}
> ROOT.gInterpreter.Declare("""
class A {};
void foo(A* a) {}
""")

> ROOT.foo(ROOT.nullptr)  # ok

> ROOT.foo(None)          # fails
TypeError: void ::foo(A* a) =>
TypeError: could not convert argument 1
~~~

- Old PyROOT has `ROOT.Long` and `ROOT.Double` to pass integer and floating point numbers by reference. In the new
PyROOT, `ctypes` must be used instead.
~~~ {.python}
> ROOT.gInterpreter.Declare("""
void foo(int& i) { ++i; }
void foo(double& d) { ++d; }
""")

> import ctypes

> i = ctypes.c_int(1)
> d = ctypes.c_double(1.)

> ROOT.foo(i); i
c_int(2)

> ROOT.foo(d); d
c_double(2.0)
~~~

- When a character array is converted to a Python string, the new PyROOT only considers the characters before the
end-of-string character:
~~~ {.python}
> ROOT.gInterpreter.Declare('char MyWord[] = "Hello";')

> mw = ROOT.MyWord

> type(mw)
<class 'str'>

> mw  # '\x00' is not part of the string
'Hello'
~~~

- Any Python class derived from a base C++ class now requires the base class to define a virtual destructor:
~~~ {.python}
> ROOT.gInterpreter.Declare("class CppBase {};")
 True
> class PyDerived(ROOT.CppBase): pass
Traceback (most recent call last):
File "<stdin>", line 1, in <module>
TypeError: CppBase not an acceptable base: no virtual destructor
~~~

- There are several name changes for what concerns cppyy APIs and proxy object attributes, including:

| Old PyROOT/cppyy                          | New PyROOT/cppyy                |
|-------------------------------------------|---------------------------------|
| cppyy.gbl.MakeNullPointer(klass)          | cppyy.bind\_object(0, klass)    |
| cppyy.gbl.BindObject / cppyy.bind\_object | cppyy.bind\_object              |
| cppyy.AsCObject                           | libcppyy.as\_cobject            |
| cppyy.add\_pythonization                  | cppyy.py.add\_pythonization     |
| cppyy.compose\_method                     | cppyy.py.compose\_method        |
| cppyy.make\_interface                     | cppyy.py.pin\_type              |
| cppyy.gbl.nullptr                         | cppyy.nullptr                   |
| cppyy.gbl.PyROOT.TPyException             | cppyy.gbl.CPyCppyy.TPyException |
| buffer.SetSize(N)                         | buffer.reshape((N,))            |
| obj.\_\_cppname\_\_                       | obj.\_\_cpp\_name\_\_           |
| obj.\_get\_smart\_ptr                     | obj.\_\_smartptr\_\_            |
| callable.\_creates                        | callable.\_\_creates\_\_        |
| callable.\_mempolicy                      | callable.\_\_mempolicy\_\_      |
| callable.\_threaded                       | callable.\_\_release\_gil\_\_   |

- New PyROOT does not parse command line arguments by default anymore. For example, when running
`python my_script.py -b`, the `-b` argument will not be parsed by PyROOT, and therefore the
batch mode will not be activated. If the user wants to enable the PyROOT argument parsing again,
they can do so by starting their Python script with:
~~~ {.python}
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
~~~

- In new PyROOT, `addressof` should be used to retrieve the address of fields in a struct,
for example:
~~~ {.python}
> ROOT.gInterpreter.Declare("""
struct MyStruct {
  int a;
  float b;
};
""")

> s = ROOT.MyStruct()

> ROOT.addressof(s, 'a')
94015521402096L

> ROOT.addressof(s, 'b')
94015521402100L

> ROOT.addressof(s)
94015521402096L
~~~
In old PyROOT, `AddressOf` could be used for that purpose too, but its behaviour was inconsistent.
`AddressOf(o)` returned a buffer whose first position contained the address of object `o`, but
`Address(o, 'field')` returned a buffer whose address was the address of the field, instead
of such address being contained in the first position of the buffer.
Note that, in the new PyROOT, `AddressOf(o)` can still be invoked, and it still returns a buffer
whose first position contains the address of object `o`.

- In old PyROOT, there were two C++ classes called `TPyMultiGenFunction` and `TPyMultiGradFunction`
which inherited from `ROOT::Math::IMultiGenFunction` and `ROOT::Math::IMultiGradFunction`, respectively.
The purpose of these classes was to serve as a base class for Python classes that wanted to inherit
from the ROOT::Math classes. This allowed to define Python functions that could be used for fitting
in Fit::Fitter.
In the new PyROOT, `TPyMultiGenFunction` and `TPyMultiGradFunction` do not exist anymore, since their
functionality is automatically provided by new cppyy: when a Python class inherits from a C++ class,
a wrapper C++ class is automatically generated. That wrapper class will redirect any call from C++
to the methods implemented by the Python class. Therefore, the user can make their Python function
classes inherit directly from the ROOT::Math C++ classes, for example:
~~~ {.python}
import ROOT

from array import array

class MyMultiGenFCN( ROOT.Math.IMultiGenFunction ):
    def NDim( self ):
        return 1

    def DoEval( self, x ):
        return (x[0] - 42) * (x[0] - 42)

    def Clone( self ):
        x = MyMultiGenFCN()
        ROOT.SetOwnership(x, False)
        return x

def main():
    fitter = ROOT.Fit.Fitter()
    myMultiGenFCN = MyMultiGenFCN()
    params = array('d', [1.])
    fitter.FitFCN(myMultiGenFCN, params)
    fitter.Result().Print(ROOT.std.cout, True)

if __name__ == '__main__':
    main()
~~~

- Inheritance of Python classes from C++ classes is not working in some cases. This is described
in [ROOT-10789](https://sft.its.cern.ch/jira/browse/ROOT-10789) and
[ROOT-10582](https://sft.its.cern.ch/jira/browse/ROOT-10582).
This affects the creation of GUIs from Python, e.g. in the
[Python GUI tutorial](https://root.cern.ch/doc/master/gui__ex_8py.html), where the inheritance
from `TGMainFrame` is not working at the moment. Future releases of ROOT will fix these
issues and provide a way to program GUIs from Python, including a replacement for TPyDispatcher,
which is no longer provided.

- When iterating over an `std::vector<std::string>` from Python, the elements returned by
the iterator are no longer of type Python `str`, but `cppyy.gbl.std.string`. This is an
optimization to make the iteration faster (copies are avoided) and it allows to call
modifier methods on the `std::string` objects.
~~~ {.python}
> import cppyy

> cppyy.cppdef('std::vector<std::string> foo() { return std::vector<std::string>{"foo","bar"};}')

> v = cppyy.gbl.foo()

> type(v)
<class cppyy.gbl.std.vector<string> at 0x4ad8220>

> for s in v:
...   print(type(s))  # s is no longer a Python string, but an std::string
...
<class cppyy.gbl.std.string at 0x4cd41b0>
<class cppyy.gbl.std.string at 0x4cd41b0>
~~~


## Tutorials

- Add 3 new TMVA tutorials:
   - ``TMVA_Higgs_Classification.C`` showing the usage of the TMVA deep neural network in a typical classification problem,
   - ``TMVA_CNN_Classification.C`` showing the usage of CNN in TMVA with `MethodDL` and in PyMVA using `MethodPyKeras`,
   - ``TMVA_RNN_Classification.C`` showing the usage of RNN in both TMVA and Keras.

## Class Reference Guide

- A new version switching mechanism has been introduced; see https://github.com/root-project/doxyvers


## Build, Configuration and Testing Infrastructure

- By default, ROOT now falls back to the built-in version of xrootd if it can't find it in the system.
This means that passing `-Dbuiltin_xrootd=ON` is not necessary anymore to build ROOT with xrootd support.
Note that built-in xrootd requires a working network connection.
- Glew library was moved from ROOT graf3d sources to the proper ROOT builtins mechanism (GLEW library and headers are still installed together with ROOT for backward compatibility)
- ROOT uses header files from source directories to build own ROOT libraries; header files from $ROOTSYS/include directory only required to compile users code.
- Updated version of Glew used as buitin mechanism in ROOT: 2.1.0 (http://glew.sourceforge.net/)
- Updated version of googletest library: 1.10.0

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

### Multi-Python PyROOT

In 6.22, the new (experimental) PyROOT is built by default. In order to build with the old PyROOT instead, the option
`-Dpyroot_legacy=ON` can be used. This is a summary of the PyROOT options:
* `pyroot`: by default `ON`, it enables the build of PyROOT.
* `pyroot_legacy`: by default `OFF`, it allows the user to select the old PyROOT (legacy) to be built instead
of the new one.
* `pyroot_experimental`: this option is **deprecated** in 6.22 and should no longer be used. If used, it triggers
a warning.

This new PyROOT also introduces the possibility of building its libraries for both Python2 and Python3 in a single
ROOT build (CMake >= 3.14 is required). If no option is specified, PyROOT will be built for the most recent
Python3 and Python2 versions that CMake can find. If only one version can be found, PyROOT will be built for only
that version. Moreover, for a given Python installation to be considered, it must provide both the Python interpreter
(binary) and the development package.

The user can also choose to build ROOT only for a given Python version, even if multiple installations exist in
the system. For that purpose, the option `-DPYTHON_EXECUTABLE=/path/to/python_exec` can be used to point to the
desired Python installation.

If the user wants to build PyROOT for both Python2 and Python3, but not necessarily with the highest versions that
are available on the system, they can provide hints to CMake by using `-DPython2_ROOT_DIR=python2_dir` and/or
`-DPython3_ROOT_DIR=python3_dir` to point to the root directory of some desired Python installation. Similarly,
`Python2_EXECUTABLE` and/or `Python3_EXECUTABLE` can be used to point to particular Python executables.

When executing a Python script, the Python version used will determine which version of the PyROOT libraries will
be loaded. Therefore, once the ROOT environment has been set (e.g. via `source $ROOTSYS/bin/thisroot.sh`), the user
will be able to use PyROOT from any of the Python versions it has been built for.

Regarding `TPython`, its library is built only for the highest Python version that PyROOT is built with.
Therefore, in a Python3-Python2 ROOT build, the Python code executed with `TPython` must be Python3-compliant.


### Bugs and Issues fixed in this release

* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-8799'>ROOT-8799</a>] - Automatically convert Python iterables into C++ collections
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10346'>ROOT-10346</a>] - [DF] Warn on mismatch between slot pool size and effective number of slots
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-2566'>ROOT-2566</a>] - `RooNDKeysPdf` has no default constructor which causes a crash when reading it from a workspace
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-2936'>ROOT-2936</a>] - Problem plotting slices of `RooSimultaneous`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-3579'>ROOT-3579</a>] - `RooTreeDataStore` not Cloning the tree properly (and const correctness)
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-4060'>ROOT-4060</a>] - histfactory - print error and exit if channel name begins with numeric char
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-4312'>ROOT-4312</a>] - NLL problem using setData
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-4580'>ROOT-4580</a>] - `RooDataSet` `reduce(Cut(mycut))` and `Draw(myvar,mycut)` give different results
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7485'>ROOT-7485</a>] - Incorrect normalization when using named ranges in fitTo() and fitting several times
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7520'>ROOT-7520</a>] - Segmentation fault when loading a `RooLinearVar` from a `RooWorkspace`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7578'>ROOT-7578</a>] - HistFactory: invalid Sample names cause segfault
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7921'>ROOT-7921</a>] - component(s) ... already in the workspace when using PROD after EDIT
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7983'>ROOT-7983</a>] - Missing `RooDataSet::createHistogram` methods
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-7986'>ROOT-7986</a>] - Severe memory leak reading `RooFitResult` and `HypoTestInverterResult` from disk
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-9535'>ROOT-9535</a>] - Unexpected crash when cloning `TH2` histogram
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-9556'>ROOT-9556</a>] - [DF] Users might silently get wrong results when they process two trees as friends, one of which has been produced by a multi-thread Snapshot of the other
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-9768'>ROOT-9768</a>] - Char array branches not correctly retrieved as Python strings
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-9966'>ROOT-9966</a>] - Thread-safety problem in `TFile` from use of `TSystem`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10109'>ROOT-10109</a>] - Wrong overload resolution for `GetBinErrorUp` in `TH2D`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10146'>ROOT-10146</a>] - `TH2` `FitSlicesX()` and `TH3` `FitSlicesZ()` don't use the axis range set by user
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10159'>ROOT-10159</a>] - from ROOT import <something> on python cmd line fails if IPython is already imported
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10175'>ROOT-10175</a>] - Notebook `jsmva on` package util not found
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10178'>ROOT-10178</a>] - `TTreeProcessorMT` can't deal with trees with different names in the same `TChain`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10209'>ROOT-10209</a>] - ROOT 6.18/00 does not compile with CUDA available in C++17 mode
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10277'>ROOT-10277</a>] - Crash when generating a large number of datasets
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10397'>ROOT-10397</a>] - [TTreeReaderArray] Segfault when accessing data from leaflist
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10417'>ROOT-10417</a>] - TMVA test bug: rtensor.cxx:303
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10458'>ROOT-10458</a>] - `RDataFrame`: sub-par diagnostics for interpreted types
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10488'>ROOT-10488</a>] - `TClass::GetListOfMethods` does not include constructor for `std::hash`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10494'>ROOT-10494</a>] - `TChain::Add` and `TChain::AddFile` broken with URL that contains double slashes
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10503'>ROOT-10503</a>] - `RooAbsPdf::chi2FitTo` overloads get lost in PyROOT
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10508'>ROOT-10508</a>] - `RDataFrame` created from `TChain` after calling `TChain::AddFriend()` crashes when trying to plot a variable if MultiThread is enabled
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10517'>ROOT-10517</a>] - `RooFit::CutRange` is buggy with multiple range
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10518'>ROOT-10518</a>] - `chi2FitTo` and `fitTo` are not working properly with multiple ranges
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10524'>ROOT-10524</a>] - Performance anomaly when running `TChain::AddFriend()`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10527'>ROOT-10527</a>] - `RDataFrame` `Display()` crashes with column of `std::string`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10529'>ROOT-10529</a>] - `Error in <TClass::ReadRules()>: Cannot find rules etc/class.rules`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10533'>ROOT-10533</a>] - [TTreeProcessorMT] Support friend chains with differently named trees
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10545'>ROOT-10545</a>] - `.qqqqqq` does not exit root 
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10549'>ROOT-10549</a>] - `RDataFrame` warns but does not throw when skipping files after network glitch
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10550'>ROOT-10550</a>] - `plotOn` normalization crazy and error message when doing loops of toys
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10561'>ROOT-10561</a>] - [TTreeProcessorMT] Breaks without `EnableImplicitMT` on
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10563'>ROOT-10563</a>] - `RDataFrame::Cache` fails when used with `RDataFrame::Alias`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10568'>ROOT-10568</a>] - Wrong version in master branch
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10572'>ROOT-10572</a>] - roofit `plotOn` reproducibly crashes when having done many toys before, in a specific configuration
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10573'>ROOT-10573</a>] - `Snapshot` fails to update `TTree` if Multithreading enabled
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10575'>ROOT-10575</a>] - Issues with enabling `pyroot_experimental`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10576'>ROOT-10576</a>] - Error getting list of methods of std in cxx17
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10578'>ROOT-10578</a>] - `TMVA::Experimental` Intel compiler support
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10583'>ROOT-10583</a>] - Incorrect OpenGL rendering
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10587'>ROOT-10587</a>] - [Doc] Python RDF tutorials not shown correctly in doxygen
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10595'>ROOT-10595</a>] - Bias obtained in `TH1::KolmogorovTest` when using option X
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10596'>ROOT-10596</a>] - [DF] Deprecate/hide `ROOT::MakeRootDataFrame`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10597'>ROOT-10597</a>] - ROOT 6.20.0 fails to compile with `-Dcuda=yes` when cudnn not installed
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10599'>ROOT-10599</a>] - `rootls -t` doesn't work with multiple cycles for a tree
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10601'>ROOT-10601</a>] - Problem Drawing Histograms after `TH2::FitSlicesX`/`Y`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10605'>ROOT-10605</a>] - wrong slice in plotting simultaneous pdf built with `RooSimWSTool`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10606'>ROOT-10606</a>] - current master fails basic `RooFit` functionality due to cppyy changes
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10619'>ROOT-10619</a>] - [DF] Crash with unused jitted filters
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10624'>ROOT-10624</a>] - [Graphics][JSROOT] Broken histogram fill color with `%jsroot on`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10629'>ROOT-10629</a>] - clang CodeGen assertion "popping exception stack when not empty"
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10630'>ROOT-10630</a>] - Failed to run python code: `copy_string=history.history.keys()[0]`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10647'>ROOT-10647</a>] - Missing implementation of "section" in `RooArgSet::writeToStream`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10648'>ROOT-10648</a>] - Occasional streamer type mismatch in multi-thread `TTree` writing
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10654'>ROOT-10654</a>] - `rootmv` doesn't work with multiple cycles `TTrees` in `TFiles`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10656'>ROOT-10656</a>] - `TTreeProcessorMT` leaves implicit multi-threading enabled after a call to `Process`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10658'>ROOT-10658</a>] - wrong `-I` flags provided for dictionary compilation
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10661'>ROOT-10661</a>] - `rootcling` produces uncompilable dictionary if class depends on boost headers
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10663'>ROOT-10663</a>] - ROOT dictionary for ATLAS persistent class can't load correctly
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10668'>ROOT-10668</a>] - `RooFit` with Asymptotically Correct approach will segfault if `RooRealVar` `Name != Title`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10670'>ROOT-10670</a>] - `hesse` does not know which minimizer to use
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10672'>ROOT-10672</a>] - segfault when retrieving `TTreeCache` from `TChain`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10676'>ROOT-10676</a>] - Constructing a `RooDataSet` from `TTree` with cuts causes warnings
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10680'>ROOT-10680</a>] - [PyROOT] numpy vs `RVec<bool>` comparison results in infinite loop
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10681'>ROOT-10681</a>] - New `cppyy` only provides `addressof`, not `AddressOf`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10689'>ROOT-10689</a>] - cling misidentifies called lambda
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10693'>ROOT-10693</a>] - [PyROOT] Infinite iteration over `RVec`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10696'>ROOT-10696</a>] - `TDirectory::GetObject()` with wrong type leaks memory
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10702'>ROOT-10702</a>] - [TTree] Wrong data could be silently written if data-member of object has same name as another branch
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10706'>ROOT-10706</a>] - Fit to `TMultiGraph` with single `TGraphAsymErrors` different from single graph fit
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10709'>ROOT-10709</a>] - RooFit depends from MathMore when ROOT compiled with `-Dmathmore=ON`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10711'>ROOT-10711</a>] - X, Y axis scale changes during allinging X,Y axis title as centered in `TMultiGraph`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10712'>ROOT-10712</a>] - Failing template instantiation during tear down.
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10713'>ROOT-10713</a>] - Several memory leaks when reading `RooStats::HypoTestInverterResult`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10718'>ROOT-10718</a>] - `TMath::Poisson(n,mu)` for `n = 0` and `mu < 0` returns wrong value !
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10721'>ROOT-10721</a>] - Python tutorials are not compatible with Python 3
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10723'>ROOT-10723</a>] - [TChain] Potential access to invalid `fTree` in `TChain` destructor
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10731'>ROOT-10731</a>] - [PyROOT] Broken conversion from boolean numpy array to `bool*`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10732'>ROOT-10732</a>] - `hadd` broken
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10733'>ROOT-10733</a>] - gtest-tree-treeplayer-test-treeprocessormt fails without xrootd
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10734'>ROOT-10734</a>] - pymva tutorials not disable when pymva not available
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10736'>ROOT-10736</a>] - roofit tutorials failing without PyRoot-experimental.
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10737'>ROOT-10737</a>] - `RooAbsPdf` is distorted outside its "Range" used previously with `fitTo`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10751'>ROOT-10751</a>] - ROOT Fails to find C++ standard library path on non-standard locale
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10753'>ROOT-10753</a>] - [TTreeReader] Wrong entries are loaded in case of `TChain`+`TEntryList`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10754'>ROOT-10754</a>] - [cling] Miscompilation with `make_shared`/`shared_ptr`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10764'>ROOT-10764</a>] - `RooAbsAnaConvPdf` crashes when being imported and retrieved from `RooWorkspace`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10765'>ROOT-10765</a>] - C++14 now requires root7
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10766'>ROOT-10766</a>] - Missing overload in LHCb's `FuncOps__::__le__`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10768'>ROOT-10768</a>] - Document backwards incompatible changes & new features in the new PyROOT
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10774'>ROOT-10774</a>] - `rootbrowse` freezes after startup
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10777'>ROOT-10777</a>] - Cannot py-import Gaudi: assert "partial specialization scope specifier in SFINAE context"
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10779'>ROOT-10779</a>] - HistFactory models that are written to a file, then retrieved with updated histograms find only old histograms
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10782'>ROOT-10782</a>] - With gcc10 STL headers don't include implicitly `stdexcept`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10784'>ROOT-10784</a>] - Mistake in what is reported in documentation
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10788'>ROOT-10788</a>] - PyROOT failed to load cmssw libraires
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10790'>ROOT-10790</a>] - [DF] Single-thread `Snapshot` into a directory also creates a spurious `TTree` outside of it 
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10792'>ROOT-10792</a>] - [DF] `Snapshot` of `TClonesArrays` read via `TTreeReaderArray` is broken
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10798'>ROOT-10798</a>] - ROOT 6.20/02 can't generate a dictionary for ATLAS's `EventLoop` package
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10804'>ROOT-10804</a>] - assertion in `clang::Sema::LookupSpecialMember`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10810'>ROOT-10810</a>] - Segmentation fault in pickling of weighted RooFit datasets 
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10820'>ROOT-10820</a>] - Circular includes `Rtypes.h` and `TGenericClassInfo.h`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10822'>ROOT-10822</a>] - [DF] `RVec`s of non-split branches can read from invalid addresses
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10833'>ROOT-10833</a>] - PyROOT FutureWarning in `ROOT.py`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10837'>ROOT-10837</a>] - `hadd` crashes when slow merging file with multiple array with same index
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10839'>ROOT-10839</a>] - Missing lock guard in `THashTable`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-6882'>ROOT-6882</a>] - RooFit doesn't support `ULong64_t`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10341'>ROOT-10341</a>] - Support I/O for `RooNDKeysPDF`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10534'>ROOT-10534</a>] - [RDataFrame] Provide report about number of event loops run
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-2558'>ROOT-2558</a>] - `SPlot` crashes when yields are not the top level fitparameters
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-5311'>ROOT-5311</a>] - `HistFactory` output is too verbose and cannot be controlled
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-8033'>ROOT-8033</a>] - Poisson constraint for statistical error in `HistFactory`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-9365'>ROOT-9365</a>] - Member functions' definition locations missing in Class `HistoToWorkspaceFactoryFast`
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10513'>ROOT-10513</a>] - [PyROOT exp] Use consistent names for interoperability functions
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10521'>ROOT-10521</a>] - `RooAbsData::getRange()` argument var should be const
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10657'>ROOT-10657</a>] - [DF] Speed up jitting
* [<a href='https://sft.its.cern.ch/jira/browse/ROOT-10665'>ROOT-10665</a>] - Make `TFile::ReadProcessID` thread safe


## HEAD of the v6-22-00-patches branch

These changes will be part of a future 6.22/02.

* None so far.
