# C++ Modules in ROOT. Technology Preview

*Vassil Vassilev, Oksana Shadura, Yuka Takahashi and Raphael Isemann*

## Overview

ROOT has several features which interact with libraries and require implicit
header inclusion. This can be triggered by reading or writing data on disk,
or user actions at the prompt. Often, the headers are immutable and reparsing is
redundant. C++ Modules are designed to minimize the reparsing of the same
header content by providing an efficient on-disk representation of C++ Code.

The ROOT v6.16 release comes with a preview of the module technology;
dedicated binaries have been built and can be reproduced by passing
`-Druntime_cxxmodules=On` as configure flag. The goals of this technology
preview are:
  * Gain feedback from early adoption -- the technology is being long anticipated
  by some of the users of ROOT. It improves correctness of ROOT and improves
  performance when carefully adopted.
  * Study performance bottlenecks -- the feature is designed with performance
  considerations in mind. In this document we describe the current performance
  bottlenecks and trade-offs.
  * Understand if the gradual migration policy is sufficient -- C++ Modules in
  ROOT support gradual migration. In particular, ROOT can enable C++ Modules for
  itself and still run in legacy mode for the third-party code (generating
  rootmap files and other scaffolding).


C++ Modules are here and we would like to give a brief introduction of how the
feature works, what are its pros and cons, what's the current state of the
implementation and how third-party code can use it.

Read more [[1]].

## Design Goals

  * Coherence with standard C++ -- C++ Modules TS is advancing and will be
  likely part the upcoming C++20 standard;
  * Performance -- provide performance that is competitive to ROOT with PCH and
  advance further the implementation of the C++ Modules in clang to optimize
  memory footprint and execution time;
  * Incremental adoption -- provide third-party code with an incremental
  migration process for their codebases.

## Motivation

An implementation of the modules concepts exists in the LLVM frontend Clang used
as a library by ROOT [[2]]. Clang supports the Modules TS and hosts modules 
research and development work. The implementation encourages incremental, 
bottom-up [[3]] adoption of the modules feature. Modules in Clang are designed
to work for C, C++, ObjectiveC, ObjectiveC++ and Swift. Users can enable the
modules feature without modifications in header files. The LLVM compiler allows
users to specify module interfaces in dedicated file, called *module maps files*.
A module map file expresses the mapping between a module file and a collection
of header files. If the compiler finds such file in the include paths it
automatically generates, imports and uses module files. The module map files can
be mounted using the compiler's virtual file system overlay mechanism to
non-writable production library installations.

In practice, a non-invasive *modularization*  can be done easily by introducing
a module map file.

```cpp
// A.h   
int pow2(int x) {
  return x * x;
}
```
```cpp
// B.cpp
#include "A.h" // clang rewires this to import A.
int main() {
  return pow2(42);
}
```

```cpp
// A.h module interface, aka module map file
module A {
  header "A.h"
  export * // clang exports the contents of A.h as part of module A.
}
```

A.h defines *pow2*, the module map file instructs clang to create *A.pcm* and
import it in B.cpp.

In a number of cases the module map files can be automatically generated if the
build system knows about the list of header files in every package. 


### Header parsing in ROOT

```cpp
// A.h
#include <string>
#include <vector>
template <class T, class U = int> struct AStruct {
  void doIt() { /*...*/ }
  std::string Name; 
  std::vector<U> Collection;
  // ...
};

template<class T, class U = AStruct<T>>
inline void freeFunction() { /* ... */ }
inline void do(unsigned N = 1) { /* ... */ }

```      
The associated with libA header files form libA's full descriptor. A.h,
potentially only part of the descriptor of libA, expands to more than 26000
lines of code.

```cpp
// Main.cpp
#include "A.h"
int main() {
  do();
  return 0;
}

```
Main.cpp, reuses code from libA by including libA's descriptor and links against
libA. The full descriptor can contain thousands of files expanding to millions
of lines of code -- a common case for framework libraries, for instance.

ROOT goes further and enhances C++ by allowing the following code to work without
explicitly requiring to `#include <A.h>`. Currently, ROOT's lack of support of
line `#5` is a long-standing, known limitation that is lifted with modules.


```cpp
  // ROOT prompt
 root [] AStruct<float> S0;     // #1: implicit loading of libA. Full descriptor required.
 root [] AStruct<float>* S1;    // #2: implicit loading of libA. No full descriptor required.
 root [] if (gFile) S1->doIt(); // #3: implicit loading of libA. Full descriptor required.
 root [] gSystem->Load("libA"); // #4: explicit loading of libA. No full descriptor required.
 root [] do();                  // #5: error: implicit loading of libA is currently unsupported.

```

This pattern is not only used in the ROOT prompt but in I/O hotspots such as
`ShowMembers` and `TClass::IsA`.

A naive implementation of this feature would require inclusion of all reachable
library descriptors (aka header files) at ROOT startup time. Of course this is
not feasible and ROOT inserts a set of optimizations to fence itself from the
costly full header inclusion. Unfortunately, several of them are home-grown and
in a few cases inaccurate (eg line #5) causing a noticeable technical debt.

Here we will briefly describe the three common layers of optimizations: ROOT PCH,
ROOTMAP and RDICT.

The ROOT precompiled header (PCH) reduces the CPU and memory cost for ROOT's
most used libraries. The precompiled header technology is well-understood since
decades [[4]]. It is an efficient on-disk representation of the state of the
compiler after parsing a set of headers. It can be loaded before starting the
next instance to avoid doing redundant work. At build time, rootcling (ROOT's
dictionary generator) creates such PCH file which is attached at ROOT startup
time. Its major drawback is the fact that if third-party users want to include
their libraries, they have to recompile it every time there is a change.

RDICT files store some useful information (in particular about class offsets) in
ROOT files to avoid the potentially expensive call to the interpreter if the
information is not the PCH. For example, ROOT's libGeom and other third-party
code. This is done to circumvent the costly call to `ShowMembers` which will
require parsing.

ROOTMAP files reduce parsing for code which is not in the PCH. Consider
`foo::bar` and `S` are defined in `libFoo`'s `Foo.h`:
```cpp
// Foo.h
namespace foo { struct bar{}; }
struct S{};
```

```bash
# libFoo.rootmap
{ decls }
namespace foo { }
struct S;
 
[ libFoo.so ]
# List of selected classes
class bar
struct S
```

```cpp
// G__Foo.cxx (aka libFoo dictionary)
namespace {
  void TriggerDictionaryInitialization_libFoo_Impl() {
    static const char* headers[] = {"Foo.h"}
    // More scaffolding
    extern int __Cling_Autoloading_Map;
    namespace foo{struct __attribute__((annotate("$clingAutoload$Foo.h"))) bar;}
    struct __attribute__((annotate("$clingAutoload$Foo.h"))) S;
    // More initialization scaffolding.
}
```

The code snippet bellow demonstrates the efforts which ROOT does to
avoid parsing redundant code.

```cpp
// ROOT prompt
root [] S *s;           // #1: does not require a definition.
root [] foo::bar *baz1; // #2: does not require a definition.
root [] foo::bar baz2;  // #3: requires a definition.
```

When starting up ROOT, it will locate all files with extensions \*.rootmap. It
parses the code in section {decls} and creates an internal map for the entities
defined in \[libFoo.so\] section. Upon seeing an unknown identifier, the
implementation searches in the database if this is a known entity.

Line #1 does not require a definition and the forward declaration consumed at
startup is sufficient. Parsing of `Foo.h` is not required. This comes at a cost
of having some non-trivial patches in clang to merge default function arguments
and default template arguments. The design of the the ROOTMAP infrastructure
requires the default arguments to be attached to more than one declaration which
is not allowed by standard C++. The behavior of line #1 is equivalent to:
```cpp
// ROOT prompt
root [] namespace foo { };struct S;
root [] S *s;
```

Line #2 does not require a definition, however, the second identifier lookup
fails. The implementation knows that `foo::bar` is in *libFoo*. It `dlopen`s
libFoo which in turn, during its static initialization, inserts annotated forward
declaration as shown in `G__Foo.cxx`. In turn, this resolves `foo::bar` and
parsing of `Foo.h` is again avoided at relatively small overhead. However, this
is very hard to measure because the dictionary of each library can have different
amount of content. In the case where the library is big and the annotated
forward declarations are many, and we want to include a relatively small header
file it may not pay off. Moreover, the loading of the annotated forward
declarations can happen at any time during parsing. This is nick-named
"recursive parsing" and is a code path that exists only in ROOT, never exercised
by clang itself and is thus not well tested. The behavior of line #2 is
equivalent to:
```cpp
// ROOT prompt
root [] namespace foo { };struct S;
root [] foo::bar/*store parsing state*/
        gSystem->Load("Foo");
        // More scaffolding.
        extern int __Cling_Autoloading_Map;
        namespace foo{struct __attribute__((annotate("$clingAutoload$Foo.h"))) bar;}
        struct __attribute__((annotate("$clingAutoload$Foo.h"))) S;
        // More initialization scaffolding.
        /*restore parsing state*/ *baz1;
```

Line #3 requires a definition and the implementation behaves exactly as in #2.
Then it is informed that a definition is required, it reads the information in
the annotation and parses `Foo.h`. The recursive parsing happens at two places
making this code path error prone.
```cpp
// ROOT prompt
root [] namespace foo { };struct S;
root [] foo::bar/*store parsing state*/
        gSystem->Load("Foo");
        // More scaffolding.
        extern int __Cling_Autoloading_Map;
        namespace foo{struct __attribute__((annotate("$clingAutoload$Foo.h"))) bar;}
        struct __attribute__((annotate("$clingAutoload$Foo.h"))) S;
        // More initialization scaffolding.
        /*restore parsing state*/ baz1 /*store parsing state*/
        #include <Foo.h>/*restore parsing state*/;
```

To recap, unfortunately, ROOT PCH is not extendable; ROOTMAP requires a lot of
maintenance and goes on a very untested codepath, while RDICT has a very limited
scope. The three features require a lot of mechanisms to work together and the
corner cases are very many. The interaction between some of the features often
break design and introduce layering violations.


## From C++ Modules to Dictionaries

C++ Modules have native capabilities to avoid reparsing. It combines all
home-grown solutions to avoid the costly operation at industry quality.

Currently, when ROOT is built with `-Druntime_cxxmodules=On` it gives priority to
C++ Module files (real *pcm* files). If such a file is present it reads all
necessary information from it. If no such file is present ROOT falls back to the
standard information flow.

### Observable differences from 'standard' ROOT

As always, ROOT is (mostly) API and ABI compatible. C++ Modules-aware ROOT is no
different. There are several differences which can be noticed:
  * \*modulemap files in $ROOTSYS/include -- those files are used by rootcling to
  put a set of header files in a single pcm file. For example, all related
  headers of *libGeom* are persisted in *Geom.pcm*. There are a few notable
  examples, which are specific to the way we build ROOT. In certain cases we
  want some header files to be compiled within C context or with RTTI on/off.
  That's mostly for bootstrapping ROOT (aka rootcling stage1).
  * modulemap.overlay.yaml -- automatically created virtual filesystem overlay
  file. This file introduces C++ Modules for external dependencies.
  For example, to 'modularize' glibc for ROOT we would need to place a modulemap
  file in (usually) `/usr/include`. This folder is not writable on many
  platforms. The vfs file tells the compiler to pretend there is a file at a
  specific location. This way we 'mount' `/usr/include/module.modulemap`
  non-invasively. The reasons why we need to extend the C++ modules support
  beyond ROOT is described bellow.
  * rootcling creates a new binary artifact *Name.pcm* after the library name --
  this is a temporary solution for the current technology preview. Once we
  advance further the implementation we will only create Name.pcm without the
  other 2 artifacts. At a final stage, ROOT might be able to integrate the
  Name.pcm with the shared library itself.
  * Preloads all \*pcm files at start up time -- this currently is the only
  remaining bottleneck which introduces a relatively small performance overhead
  at startup time and is described bellow. It will be negligible for third-
  party code (dominated by header parsing).
  * Improved correctness in number of cases -- in a few cases ROOT is more
  correct. In particular, when resolving global variables and function
  declarations which are not part of the ROOT PCH.
  * Enhanced symbol resolution mechanisms, bloom filters -- standard ROOT relies
  on information in ROOTMAP files to react when the llvm JIT issues an
  unresolved symbol callback. C++ Modules-aware ROOT relies on a behavior much
  closer to the standard linker behavior. In particular, we start searching on
  the LD_LIBRARY_PATH descending to the system libraries. The algorithm is very
  efficient because it uses bloom filters[[5]]. This in turn allows ROOT symbol
  to be extended to system libraries.
  
### Supported Platforms

  We support all platforms with glibc++ versions: 5.2, 6.2 and 7.2 and 8.

  We have experimental support for OSX XCode 10.

## Changes required by the users
  * Self-contained header files -- every header file should be able to compile
  on its own. For instance, `gcc -fsyntax-only -xc++ header.h`
  * Enable it in `rootcling` -- rootcling can produce a C++ Modules-aware
  dictionary when it is invoked with `-cxxmodule` flag.
  * Modularization of external dependencies -- if a header file is not explicitly
  nominated as part of a module and it is transitively included in two modules,
  both modules contain that header file content. In other words, the header is
  duplicated. In turn, this leads to performance regressions. If a dictionary
  depends on a header (directly or indirectly) from a external library (e.g.
  libxml) it needs to be modularized. As part of our ongoing efforts to move
  CMSSW to use C++ Modules [[6]] we have implemented a helper tool [[7]]. The
  tool detects (based on the include paths of the compiler) dependencies and
  tries to generate the relevant vfs file.

## State of the union

C++ Modules-aware ROOT preloads all modules at start up time. Our motivating
example:

```cpp
// ROOT prompt
root [] S *s;           // #1: does not require a definition.
root [] foo::bar *baz1; // #2: does not require a definition.
root [] foo::bar baz2;  // #3: requires a definition.
```

becomes equivalent to

```cpp
// ROOT prompt
root [] import ROOT.*;
root [] import Foo.*;
root [] S *s;           // #1: does not require a definition.
root [] foo::bar *baz1; // #2: does not require a definition.
root [] foo::bar baz2;  // #3: requires a definition.
```

The implementation avoids recursive actions and relies on a well-defined (by
the C++ standard) behavior. Currently, this comes with a constant performance
overhead which we go in details bellow.


### Current limitations
  * Relocatability issues -- we have fixed a few of the relocatability issues we
  found. We are aware of an obscure relocatability issue when ROOT is copied in
  another folder and we are rebuild. ROOT picks up both modulemap files in
  seemingly distinct locations.
  * Building pcms with rootcling -- in rare cases there might be issues when
  building pcm files with rootcling. The easiest will be to open a bug report
  to clang, however, reproducing a failure outside of rootcling is very difficult
  at the moment.
  * Generation of RooFit dictionary hangs -- on some platforms (depending on the
  version of libstdc++) the generation of the RooFit dictionary goes in an
  infinite loop. We have fixed a number of such situations. Please contact us if
  you see such behavior or disable roofit (`-Droofit=Off`).

### Performance
This section compares ROOT PCH technology with C++ Modules which is important but
unfair comparison. As we noted earlier, PCH is very efficient, it cannot be
extended to the experiments’ software stacks because of its design constraints.
On the contrary, the C++ Modules can be used in third-party code where the PCH
is not available.

The comparisons are to give a good metric when we are ready to switch ROOT to use
C++ Modules by default. However, since it is essentially the same technology,
optimizations of C++ Modules also affect the PCH. We have a few tricks up in
the slaves to but they come with given trade-offs. For example, we can avoid
preloading of all modules at the cost of introducing recursive behavior in
loading. This requires to build a global module index which is an on-disk
hash table. It will contain information about the mapping between an
identifier and a module name. Upon failed identifier lookup we will use the
map to decide which set of modules should be loaded. Another optimization
includes building some of the modules without `-fmodules-local-submodule-visibility`.
In turn, this would flatten the C++ modules structure and give us performance
comparable to the ROOT PCH. The trade-off is that we will decrease the
encapsulation and leak information about implementation-specific header files.

The main focus for this technology preview was not in performance due to
time considerations. We have invested some resources in optimizations and
we would like to show you (probably outdated) preliminary performance
results:

  * Memory footprint -- mostly due to importing all C++ Modules at startup
  we see overhead which depends on the number of preloaded modules. For
  ROOT it is between 40-60 MB depending on the concrete configuration.
  When the workload increases we notice that the overall memory performance
  decreases in number of cases.
  * Execution times -- likewise we have an execution overhead. For 
  workflows which take ms the slowdown can be 2x. Increasing of the work
  to seconds shows 50-60% slowdowns.

The performance of the technology preview is dependent on many factors such
as configuration of ROOT and workflow. You can read more at our Intel
IPCC-ROOT Showcase presentation here (pp 25-33)[[8]].

You can visit our continuous performance monitoring tool where we compare
the performance of the technology preview with respect to 'standard' ROOT[[9]].
*Note: if you get error 400, clean your cache or open a private browser session.*

## How to use
  Compile ROOT with `-Druntime_cxxmodules=On`. Enjoy.

## Citing ROOT's C++ Modules
```latex
% Peer-Reviewed Publication
%
% 22nd International Conference on Computing in High Energy and Nuclear Physics (CHEP)
% 8-14 October, 2016, San Francisco, USA
%
@inproceedings{Vassilev_ROOTModules,
  author = {Vassilev,V.},
  title = {{Optimizing ROOT's Performance Using C++ Modules}},
  journal = {Journal of Physics: Conference Series},
  year = 2017,
  month = {oct},
  volume = {898},
  number = {7},
  pages = {072023},
  doi = {10.1088/1742-6596/898/7/072023},
  url = {https://iopscience.iop.org/article/10.1088/1742-6596/898/7/072023/pdf},
  publisher = {{IOP} Publishing}
}
```
  
# Acknowledgement

We would like to thank the ROOT team.

We would like to thank Liz Sexton-Kennedy (FNAL) in particular for supporting
this project.

We would like to thank Axel Naumann for early feedback on this document.

This work has been supported by an Intel Parallel Computing Center grant, by U.S.
National Science Foundation grants PHY-1450377 and PHY-1624356, and by the U.S.
Department of Energy, Office of Science.


# References
(1): [Vassilev, V., 2017, October. Optimizing ROOT's Performance Using C++ Modules. In Journal of Physics: Conference Series (Vol. 898, No. 7, p. 072023). IOP Publishing.][1]

(2): [Clang Modules, Official Documentation][2]

(3): [Manuel Klimek, Deploying C++ Modules to 100s of Millions of Lines of Code, 2016, CppCon][3]

(4): [Precompiled Header and Modules Internals, Official Documentation][4]

(5): [Bloom Filter][5]

(6): [C++ Modules support (based on Clang), GitHub Repo][5]

(7): [Make your third party libraries modular][6]

(8): [Vassilev, V., Shadura, O., Takahashi, Y., IPCC-ROOT Showcase Presentation, Nov, 2018][7]

(9): [ROOT Continuous Performance Monitoring System][8]



[//]: # (Links)
[1]: https://www.researchgate.net/profile/Vassil_Vassilev3/publication/319717664_Optimizing_ROOT%27s_Performance_Using_C_Modules/links/59bad690aca272aff2d01c1c/Optimizing-ROOTs-Performance-Using-C-Modules.pdf "Vassilev, V., 2017, October. Optimizing ROOT’s Performance Using C++ Modules. In Journal of Physics: Conference Series (Vol. 898, No. 7, p. 072023). IOP Publishing."

[2]: https://clang.llvm.org/docs/Modules.html "Clang Modules"

[3]: https://cppcon2016.sched.com/event/7nM2/deploying-c-modules-to-100s-of-millions-of-lines-of-code "Deploying C++ Modules to 100s of Millions of Lines of Code"

[4]: https://clang.llvm.org/docs/PCHInternals.html "Precompiled Header and Modules Internals"

[5]: https://en.wikipedia.org/wiki/Bloom_filter "Bloom Filter"

[6]: https://github.com/cms-sw/cmssw/issues/15248 "C++ Modules support (based on Clang)"

[7]: https://github.com/Teemperor/ClangAutoModules "Make your third party libraries modular"

[8]: https://ipcc-root.github.io/downloads/20181108-ipcc-princeton-showcase-presentation.pdf "IPCC-ROOT Showcase Presentation"

[9]: http://root-bench.cern.ch/ "ROOT Continuous Performance Monitoring System"

