/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// ----------------------------------------------------------------------
// HEP coherent system of Units
//
// This file has been provided to CLHEP by Geant4 (simulation toolkit for HEP).
// Adapted to TGeo units base by Marko Petric
//
// The basic units are :
// millimeter              (millimeter)
// nanosecond              (nanosecond)
// Mega electron Volt      (MeV)
// positron charge         (eplus)
// degree Kelvin           (kelvin)
// the amount of substance (mole)
// luminous intensity      (candela)
// radian                  (radian)
// steradian               (steradian)
//
// Below is a non exhaustive list of derived and pratical units
// (i.e. mostly the SI units).
// You can add your own units.
//
// The SI numerical value of the positron charge is defined here,
// as it is needed for conversion factor : positron charge = e_SI (coulomb)
//
// The others physical constants are defined in the header file :
// PhysicalConstants.h
//
// Authors: M.Maire, S.Giani
//
// History:
//
// 06.02.96   Created.
// 28.03.96   Added miscellaneous constants.
// 05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
// 20.05.98   names: meter, second, gram, radian, degree
//            (from Brian.Lasiuk@yale.edu (STAR)). Added luminous units.
// 05.08.98   angstrom, picobarn, microsecond, picosecond, petaelectronvolt
// 01.03.01   parsec
// 31.01.06   kilogray, milligray, microgray
// 29.04.08   use PDG 2006 value of e_SI
// 03.11.08   use PDG 2008 value of e_SI
// 19.08.15   added liter and its sub units (mma)
// 12.01.16   added symbols for microsecond (us) and picosecond (ps) (mma)
// 02.10.17   addopted units from CLHEP 2.3.4.3 and converted to TGeo unit base

#ifndef TGEANT4_SYSTEM_OF_UNITS_H
#define TGEANT4_SYSTEM_OF_UNITS_H

namespace TGeant4Unit {

  //
  // TGeo follows Geant3 convention as specified in manual
  // "Unless otherwise specified, the following units are used throughout the program:
  // centimeter, second, degree, GeV"

  //
  //
  //
  static constexpr double pi = 3.14159265358979323846;
  static constexpr double twopi = 2 * pi;
  static constexpr double halfpi = pi / 2;
  static constexpr double pi2 = pi * pi;

  //
  // Length [L]
  //
  static constexpr double millimeter = 1.0;
  static constexpr double millimeter2 = millimeter * millimeter;
  static constexpr double millimeter3 = millimeter * millimeter * millimeter;

  static constexpr double centimeter = 10. * millimeter; // Base unit
  static constexpr double centimeter2 = centimeter * centimeter;
  static constexpr double centimeter3 = centimeter * centimeter * centimeter;

  static constexpr double meter = 1000. * millimeter;
  static constexpr double meter2 = meter * meter;
  static constexpr double meter3 = meter * meter * meter;

  static constexpr double kilometer = 1000. * meter;
  static constexpr double kilometer2 = kilometer * kilometer;
  static constexpr double kilometer3 = kilometer * kilometer * kilometer;

  static constexpr double parsec = 3.0856775807e+16 * meter;

  static constexpr double micrometer = 1.e-6 * meter;
  static constexpr double nanometer = 1.e-9 * meter;
  static constexpr double angstrom = 1.e-10 * meter;
  static constexpr double fermi = 1.e-15 * meter;

  static constexpr double barn = 1.e-28 * meter2;
  static constexpr double millibarn = 1.e-3 * barn;
  static constexpr double microbarn = 1.e-6 * barn;
  static constexpr double nanobarn = 1.e-9 * barn;
  static constexpr double picobarn = 1.e-12 * barn;

  // symbols
  static constexpr double nm = nanometer;
  static constexpr double um = micrometer;

  static constexpr double mm = millimeter;
  static constexpr double mm2 = millimeter2;
  static constexpr double mm3 = millimeter3;

  static constexpr double cm = centimeter;
  static constexpr double cm2 = centimeter2;
  static constexpr double cm3 = centimeter3;

  static constexpr double liter = 1.e+3 * cm3;
  static constexpr double L = liter;
  static constexpr double dL = 1.e-1 * liter;
  static constexpr double cL = 1.e-2 * liter;
  static constexpr double mL = 1.e-3 * liter;

  static constexpr double m = meter;
  static constexpr double m2 = meter2;
  static constexpr double m3 = meter3;

  static constexpr double km = kilometer;
  static constexpr double km2 = kilometer2;
  static constexpr double km3 = kilometer3;

  static constexpr double pc = parsec;

  //
  // Angle
  //
  static constexpr double degree = 1.0; // Base unit
  static constexpr double radian = (180.0 / pi) * degree;
  static constexpr double milliradian = 1.e-3 * radian;

  static constexpr double steradian = 1.;

  // symbols
  static constexpr double rad = radian;
  static constexpr double mrad = milliradian;
  static constexpr double sr = steradian;
  static constexpr double deg = degree;

  //
  // Time [T]
  //
  static constexpr double nanosecond = 1.0;
  static constexpr double second = 1.e+9 * nanosecond; // Base unit
  static constexpr double millisecond = 1.e-3 * second;
  static constexpr double microsecond = 1.e-6 * second;
  static constexpr double picosecond = 1.e-12 * second;

  static constexpr double hertz = 1. / second;
  static constexpr double kilohertz = 1.e+3 * hertz;
  static constexpr double megahertz = 1.e+6 * hertz;

  // symbols
  static constexpr double ns = nanosecond;
  static constexpr double s = second;
  static constexpr double ms = millisecond;
  static constexpr double us = microsecond;
  static constexpr double ps = picosecond;

  //
  // Electric charge [Q]
  //
  static constexpr double eplus = 1.;             // positron charge
  static constexpr double e_SI = 1.602176487e-19; // positron charge in coulomb
  static constexpr double coulomb = eplus / e_SI; // coulomb = 6.24150 e+18 * eplus

  //
  // Energy [E]
  //
  static constexpr double megaelectronvolt = 1.0;
  static constexpr double electronvolt = 1.e-6 * megaelectronvolt;
  static constexpr double kiloelectronvolt = 1.e-3 * megaelectronvolt;
  static constexpr double gigaelectronvolt = 1.e+3 * megaelectronvolt; // Base unit
  static constexpr double teraelectronvolt = 1.e+6 * megaelectronvolt;
  static constexpr double petaelectronvolt = 1.e+9 * megaelectronvolt;

  static constexpr double joule = electronvolt / e_SI; // joule = 6.24150 e+12 * MeV

  // symbols
  static constexpr double MeV = megaelectronvolt;
  static constexpr double eV = electronvolt;
  static constexpr double keV = kiloelectronvolt;
  static constexpr double GeV = gigaelectronvolt;
  static constexpr double TeV = teraelectronvolt;
  static constexpr double PeV = petaelectronvolt;

  //
  // Mass [E][T^2][L^-2]
  //
  static constexpr double kilogram = joule * second * second / (meter * meter);
  static constexpr double gram = 1.e-3 * kilogram;
  static constexpr double milligram = 1.e-3 * gram;

  // symbols
  static constexpr double kg = kilogram;
  static constexpr double g = gram;
  static constexpr double mg = milligram;

  //
  // Power [E][T^-1]
  //
  static constexpr double watt = joule / second; // watt = 6.24150 e+3 * MeV/ns

  //
  // Force [E][L^-1]
  //
  static constexpr double newton = joule / meter; // newton = 6.24150 e+9 * MeV/mm

  //
  // Pressure [E][L^-3]
  //
#define pascal hep_pascal                             // a trick to avoid warnings
  static constexpr double hep_pascal = newton / m2;     // pascal = 6.24150 e+3 * MeV/mm3
  static constexpr double bar = 100000 * pascal;        // bar    = 6.24150 e+8 * MeV/mm3
  static constexpr double atmosphere = 101325 * pascal; // atm    = 6.32420 e+8 * MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  static constexpr double ampere = coulomb / second; // ampere = 6.24150 e+9 * eplus/ns
  static constexpr double milliampere = 1.e-3 * ampere;
  static constexpr double microampere = 1.e-6 * ampere;
  static constexpr double nanoampere = 1.e-9 * ampere;

  //
  // Electric potential [E][Q^-1]
  //
  static constexpr double megavolt = megaelectronvolt / eplus;
  static constexpr double kilovolt = 1.e-3 * megavolt;
  static constexpr double volt = 1.e-6 * megavolt;

  //
  // Electric resistance [E][T][Q^-2]
  //
  static constexpr double ohm = volt / ampere; // ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  static constexpr double farad = coulomb / volt; // farad = 6.24150e+24 * eplus/Megavolt
  static constexpr double millifarad = 1.e-3 * farad;
  static constexpr double microfarad = 1.e-6 * farad;
  static constexpr double nanofarad = 1.e-9 * farad;
  static constexpr double picofarad = 1.e-12 * farad;

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  static constexpr double weber = volt * second; // weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  static constexpr double tesla = volt * second / meter2; // tesla =0.001*megavolt*ns/mm2

  static constexpr double gauss = 1.e-4 * tesla;
  static constexpr double kilogauss = 1.e-1 * tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  static constexpr double henry = weber / ampere; // henry = 1.60217e-7*MeV*(ns/eplus)**2

  //
  // Temperature
  //
  static constexpr double kelvin = 1.;

  //
  // Amount of substance
  //
  static constexpr double mole = 1.;

  //
  // Activity [T^-1]
  //
  static constexpr double becquerel = 1. / second;
  static constexpr double curie = 3.7e+10 * becquerel;
  static constexpr double kilobecquerel = 1.e+3 * becquerel;
  static constexpr double megabecquerel = 1.e+6 * becquerel;
  static constexpr double gigabecquerel = 1.e+9 * becquerel;
  static constexpr double millicurie = 1.e-3 * curie;
  static constexpr double microcurie = 1.e-6 * curie;
  static constexpr double Bq = becquerel;
  static constexpr double kBq = kilobecquerel;
  static constexpr double MBq = megabecquerel;
  static constexpr double GBq = gigabecquerel;
  static constexpr double Ci = curie;
  static constexpr double mCi = millicurie;
  static constexpr double uCi = microcurie;

  //
  // Absorbed dose [L^2][T^-2]
  //
  static constexpr double gray = joule / kilogram;
  static constexpr double kilogray = 1.e+3 * gray;
  static constexpr double milligray = 1.e-3 * gray;
  static constexpr double microgray = 1.e-6 * gray;

  //
  // Luminous intensity [I]
  //
  static constexpr double candela = 1.;

  //
  // Luminous flux [I]
  //
  static constexpr double lumen = candela * steradian;

  //
  // Illuminance [I][L^-2]
  //
  static constexpr double lux = lumen / meter2;

  //
  // Miscellaneous
  //
  static constexpr double perCent = 0.01;
  static constexpr double perThousand = 0.001;
  static constexpr double perMillion = 0.000001;
} // namespace TGeoUnit

#endif /* TGEANT4_SYSTEM_OF_UNITS_H */
