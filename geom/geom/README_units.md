# Selecting the System of Units in ROOT #

Historically the system of units in ROOT was based on the three basic units
centimeters, seconds and GigaElectronVolts. 
For the LHC era in Geant4 collaboration decided that a besic unit system based
on millimeters, nanoseconds and MegaElectronVolts was better suited for the LHC
experiments. All LHC experiments use Geant4 and effectively adopted this
convention for all areas of data processing: simulation, reconstruction and
data analysis. Hence experiments using the ROOT geometry toolkit to describe
the geometry had two different system of units in the application code.

To allow users having the same system of units in the geometry description and the
application it is now possible to choose the system of units at startup of the
application:

``` {.cpp}
TGeoManager::SetDefaultUnits(xx);  xx = kG4Units, kRootUnits
```

To ensure backwards compatibility ROOT's default system of units is - as it was before -
based on centimeters, seconds and GigaElectronVolts, ie. the defaults are equivalent to:

``` {.cpp}
TGeoManager::SetDefaultUnits(kRootUnits);
```

To avoid confusion between materials described in ROOT units and materials described
in Geant4 units, this switch should by all means be set once, before any element or
material is constructed. If for whatever reason it is necessary to change the
system of units later, this is feasible disabling the otherwise fatal exception:

``` {.cpp}
TGeoManager::LockDefaultUnits(kFALSE);
```

followed later by a corresponding call to again lock the system of units:

``` {.cpp}
TGeoManager::LockDefaultUnits(kTRUE);
```
