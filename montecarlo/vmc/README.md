# Virtual Monte Carlo

The Virtual Monte Carlo (VMC) allows to run different simulation Monte Carlo without changing the user code and therefore the input and output format as well as the geometry and detector response definition.

The core of the VMC is the category of classes [**vmc**](https://root.cern.ch/doc/master/group__vmc.html). It provides a set of interfaces which completely decouple the dependencies between the user code and the concrete Monte Carlo:

- [TVirtualMC](https://root.cern.ch/doc/master/classTVirtualMC.html): Interface to the concrete Monte Carlo program
- [TVirtualMCApplication](https://root.cern.ch/doc/master/classTVirtualMCApplication.html): Interface to the user's Monte Carlo application
- [TVirtualMCStack](https://root.cern.ch/doc/master/classTVirtualMCStack.html): Interface to the particle stack
- [TVirtualMCDecayer](https://root.cern.ch/doc/master/classTVirtualMCDecayer.html): Interface to the external decayer
- [TVirtualMCSensitiveDetector](https://root.cern.ch/doc/master/classTVirtualMCSensitiveDetector.html): Interface to the user's sensitive detector

The implementation of the TVirtualMC interface is provided for two Monte Carlo transport codes, GEANT3 and [Geant4](http://geant4.web.cern.ch/geant4/), with the VMC packages listed below. The implementation for the third Monte Carlo transport code, FLUKA, has been discontinued by the FLUKA team in 2010.

The other three interfaces are implemented in the user application.The user has to implement two mandatory classes: the **MC application** (derived from TVirtualMCApplication) and the **MC stack** (derived from TVirtualMCStack), optionally an external decayer (derived from TVirtualMCDecayer) can be introduced. The user VMC application is independent from concrete transport codes (GEANT3, Geant4, FLUKA). The transport code which will be used for simulation is selected at run time - when processing a ROOT macro where the concrete Monte Carlo is instantiated.

The relationships between the interfaces and their implementations are illustrated in the class diagrams: User MC application , Virtual MC , demonstarting the decoupling between the user code and the concrete transport code.

## VMC and TGeo

The VMC is fully integrated with the Root geometry package, [TGeo](http://root.cern.ch/root/htmldoc/GEOM_GEOM_Index.html), and users can easily define their VMC application with TGeo geometry and this way of geometry definition is recommended for new users.

It is also possible to define geometry via Geant3-like functions defined in the VMC interface, however this way is kept only for backward compatibility and should not be used by new VMC users.

## Available VMCs

### For GEANT3 - geant3
Geant3 VMC (C++) is provided within a single package together with GEANT3 (Fortran) - geant3 .

### For Geant4 - geant4_vmc
[Geant4 VMC](http://root.cern.ch/drupal/content/geant4-vmc) is provided within a package geant4_vmc , that, in difference from geant3, naturally does not include Geant4 itself and you will need the Geant4 installation to run your VMC application with Geant4 VMC.

## Multiple VMCs

Since the development version the simulation can be shared among multiple different engines deriving from [TVirtualMC](https://root.cern.ch/doc/master/classTVirtualMC.html) which are handled by a singleton [TMCManager](https://root.cern.ch/doc/master/classTMCManager.html) object.

See more detailed description in [the dedicated README](README.multiple.md).

## Documentation

[https://root.cern.ch/vmc](https://root.cern.ch/vmc)

## Authors

The concept of Virtual MonteCarlo has been developed by the [ALICE Software Project](http://aliceinfo.cern.ch/Offline/).<br>
Authors: R. Brun<sup>1</sup>, F. Carminati<sup>1</sup>, A. Gheata<sup>1</sup>, I. Hrivnacova<sup>2</sup>, A. Morsch<sup>1</sup>, B. Volkel<sup>1</sup>;<br>
<sup>1</sup>European Organization for Nuclear Research (CERN), Geneva, Switzerland;<br>
<sup>2</sup>Institut de Physique Nucléaire dʼOrsay (IPNO), Université Paris-Sud, CNRS-IN2P3, Orsay, France 

Contact: root-vmc@root.cern.ch

VMC pages maintained by: Ivana Hrivnacova <br>
*Last update: 12/04/2019*
