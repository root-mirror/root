<img src="https://root-forum.cern.ch/uploads/default/original/2X/3/3fb82b650635bc6d61461f3c47f41786afad4548.png" align="right"  height="50"/>

## About

The ROOT system provides a set of OO frameworks with all the functionality
needed to handle and analyze large amounts of data in a very efficient way.
Having the data defined as a set of objects, specialized storage methods are
used to get direct access to the separate attributes of the selected objects,
without having to touch the bulk of the data. Included are histograming
methods in an arbitrary number of dimensions, curve fitting, function
evaluation, minimization, graphics and visualization classes to allow
the easy setup of an analysis system that can query and process the data
interactively or in batch mode, as well as a general parallel processing
framework, PROOF, that can considerably speed up an analysis.

Thanks to the built-in C++ interpreter cling, the command, the
scripting and the programming language are all C++. The interpreter
allows for fast prototyping of the macros since it removes the time
consuming compile/link cycle. It also provides a good environment to
learn C++. If more performance is needed the interactively developed
macros can be compiled using a C++ compiler via a machine independent
transparent compiler interface called ACliC.

The system has been designed in such a way that it can query its databases
in parallel on clusters of workstations or many-core machines. ROOT is
an open system that can be dynamically extended by linking external
libraries. This makes ROOT a premier platform on which to build data
acquisition, simulation and data analysis systems.

[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL%20v2.1+-blue.svg)](https://www.gnu.org/licenses/lgpl.html) [![Test coverage](https://root.cern/files/img/coverage-badge.svg)](https://epsft-jenkins.cern.ch/job/root-nightly-master-coverage/cobertura)

## Build Status
| Branch | Nightly build status |
|--------|------------|
| master | [![Build Status](https://epsft-jenkins.cern.ch/buildStatus/icon?job=root-nightly-master)](https://epsft-jenkins.cern.ch/view/ROOT/job/root-nightly-master/) |
| v6-18-00-patches | [![Build Status](https://epsft-jenkins.cern.ch/buildStatus/icon?job=root-nightly-v6-18-00-patches)](https://epsft-jenkins.cern.ch/view/ROOT/job/root-nightly-v6-18-00-patches/) |
| v6-16-00-patches | [![Build Status](https://epsft-jenkins.cern.ch/buildStatus/icon?job=root-nightly-v6-16-00-patches)](https://epsft-jenkins.cern.ch/view/ROOT/job/root-nightly-v6-16-00-patches/) |

## Cite
We are [![DOI](https://zenodo.org/badge/10994345.svg)](https://zenodo.org/badge/latestdoi/10994345)

Please cite us as

    Rene Brun and Fons Rademakers, ROOT - An Object Oriented Data Analysis Framework,
    Proceedings AIHENP'96 Workshop, Lausanne, Sep. 1996,
    Nucl. Inst. & Meth. in Phys. Res. A 389 (1997) 81-86.
    See also "ROOT" [software], Release vX.YY/ZZ, dd/mm/yyyy,
    (Select the right link for your release here: https://zenodo.org/search?page=1&size=20&q=conceptrecid:848818&all_versions&sort=-version).

## Live Demo for CERN Users
[![](https://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png)](http://cern.ch/swanserver/cgi-bin/go?projurl=https://github.com/cernphsft/rootbinder.git)

## Screenshots
These screenshots shows some of the plots (produced using ROOT) presented when the Higgs boson discovery was [announced at CERN](http://home.cern/topics/higgs-boson):

![CMS Data MC Ratio Plot](https://d35c7d8c.web.cern.ch/sites/d35c7d8c.web.cern.ch/files/CMS04_1.png)

![Atlas P0 Trends](https://d35c7d8c.web.cern.ch/sites/d35c7d8c.web.cern.ch/files/Atlas06_0.png)

See more screenshots on our [gallery](https://root.cern/gallery).

## Download and Getting Started
See [root.cern download page](https://root.cern/downloading-root) for the latest binary releases.

[Getting started with ROOT.](https://root.cern/getting-started)

## Building
Clone the repo

    $ git clone https://github.com/root-project/root.git

Make a directory for building

    $ mkdir build
    $ cd build

Run cmake and make

    $ cmake ../root
    $ make -j8

Setup and run ROOT

    $ source bin/thisroot.sh
    $ root

[More information](https://root.cern/building-root) regarding building.

## Help and Support
- [Forum](https://root.cern/forum/)
- [Issue tracker](https://sft.its.cern.ch/jira/projects/ROOT/issues/ROOT-5820?filter=allopenissues)
- [Report a bug](https://root.cern/bugs) (Requires a [CERN lightweight account](https://account.cern.ch/account/Externals/RegisterAccount.aspx))
- [Mailing lists](https://groups.cern.ch/group/root-dev/default.aspx)
- [Documentation](https://root.cern/guides/reference-guide)
- [Tutorials](https://root.cern/doc/master/group__Tutorials.html)

## Contribution Guidelines
- [How to contribute](https://github.com/root-project/root/blob/master/CONTRIBUTING.md)
- [Bug reporting guidelines](https://root.cern/guidelines-submitting-bug)
- [Coding conventions](https://root.cern/coding-conventions)
- [Meetings](https://root.cern/meetings)

