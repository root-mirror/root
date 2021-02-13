/* @(#)root/gdml:$Id$ */
// Author: Ben Lloyd 09/11/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGDMLParse
\ingroup Geometry_gdml

 This class contains the implementation of the GDML  parser associated to
 all the supported GDML elements. User should never need to explicitly
 instaciate this class. It is internally used by the TGeoManager.

 Each element process has a 'Binding' to ROOT. The 'binding' is specific
 mapping of GDML elements (materials, solids, etc) to specific objects which
 should be instanciated by the converted. In the present case (ROOT) the
 binding is implemented at the near the end of each process function. Most
 bindings follow similar format, dependent on what is being added to the
 geometry.

 This file also contains the implementation of the TGDMLRefl class. This is
 just a small helper class used internally by the 'reflection' method (for
 reflected solids).

 The presently supported list of TGeo classes is the following:

#### Materials:
  - TGeoElement
  - TGeoMaterial
  - TGeoMixture

#### Solids:
  - TGeoBBox
  - TGeoArb8
  - TGeoTubeSeg
  - TGeoConeSeg
  - TGeoCtub
  - TGeoPcon
  - TGeoTrap
  - TGeoGtra
  - TGeoTrd2
  - TGeoSphere
  - TGeoPara
  - TGeoTorus
  - TGeoHype
  - TGeoPgon
  - TGeoXtru
  - TGeoEltu
  - TGeoParaboloid
  - TGeoCompositeShape (subtraction, union, intersection)

#### Approximated Solids:
  - Ellipsoid (approximated to a TGeoBBox)
  - Elliptical cone (approximated to a TGeoCone)

#### Geometry:
  - TGeoVolume
  - TGeoVolumeAssembly
  - divisions
  - reflection

When most solids or volumes are added to the geometry they


 Whenever a new element is added to GDML schema, this class needs to be extended.
 The appropriate method (process) needs to be implemented, as well as the new
 element process then needs to be linked thru the function TGDMLParse

 For any question or remarks concerning this code, please send an email to
 ben.lloyd@cern.ch

*/

#include "TGDMLParse.h"
#include "TGDMLMatrix.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TXMLEngine.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoParaboloid.h"
#include "TGeoArb8.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoTrd2.h"
#include "TGeoPcon.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TGeoTorus.h"
#include "TGeoPara.h"
#include "TGeoHype.h"
#include "TGeoEltu.h"
#include "TGeoXtru.h"
#include "TGeoScaledShape.h"
#include "TGeoTessellated.h"
#include "TMath.h"
#include "TMap.h"
#include "TObjString.h"
#include "TGeoExtension.h"
#include "TGeoMaterial.h"
#include "TGeoBoolNode.h"
#include "TGeoMedium.h"
#include "TGeoElement.h"
#include "TGeoShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoRegion.h"
#include "TGeoOpticalSurface.h"
#include "TGeoSystemOfUnits.h"
#include "TGeant4SystemOfUnits.h"

#include <cstdlib>
#include <string>
#include <sstream>
#include <locale>

ClassImp(TGDMLParse);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGDMLParse::TGDMLParse()
{
   fWorldName = "";
   fWorld = 0;
   fVolID = 0;
   fFILENO = 0;
   for (Int_t i = 0; i < 20; i++)
      fFileEngine[i] = 0;
   fStartFile = 0;
   fCurrentFile = 0;
   auto def_units = gGeoManager->GetDefaultUnits();
   switch (def_units) {
   case TGeoManager::kG4Units:
      fDefault_lunit = "mm";
      fDefault_aunit = "rad";
      break;
   case TGeoManager::kRootUnits:
      fDefault_lunit = "cm";
      fDefault_aunit = "deg";
      break;
   default: // G4 units
      fDefault_lunit = "mm";
      fDefault_aunit = "rad";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the new instance of the XMLEngine called 'gdml', using the filename >>
/// then parses the file and creates the DOM tree. Then passes the DOM to the
/// next function to translate it.

TGeoVolume *TGDMLParse::GDMLReadFile(const char *filename)
{
   // First create engine
   TXMLEngine *gdml = new TXMLEngine;
   gdml->SetSkipComments(kTRUE);

   // Now try to parse xml file
   XMLDocPointer_t gdmldoc = gdml->ParseFile(filename);
   if (gdmldoc == 0) {
      delete gdml;
      return 0;
   } else {

      // take access to main node
      XMLNodePointer_t mainnode = gdml->DocGetRootElement(gdmldoc);

      fFileEngine[fFILENO] = gdml;
      fStartFile = filename;
      fCurrentFile = filename;

      // display recursively all nodes and subnodes
      ParseGDML(gdml, mainnode);

      // Release memory before exit
      gdml->FreeDoc(gdmldoc);
      delete gdml;
   }
   return fWorld;
}

////////////////////////////////////////////////////////////////////////////////
/// This function recursively moves thru the DOM tree of the GDML file. It checks for
/// key words along the way and if a key word is found it calls the corresponding
/// function to interpret the node.

const char *TGDMLParse::ParseGDML(TXMLEngine *gdml, XMLNodePointer_t node)
{
   DefineConstants();
   XMLAttrPointer_t attr = gdml->GetFirstAttr(node);
   const char *name = gdml->GetNodeName(node);
   XMLNodePointer_t parentn = gdml->GetParent(node);
   const char *parent = gdml->GetNodeName(parentn);
   XMLNodePointer_t childtmp = 0;

   const char *posistr = "position";
   const char *setustr = "setup";
   const char *consstr = "constant";
   const char *varistr = "variable";
   const char *quanstr = "quantity";
   const char *matrstr = "matrix";
   const char *rotastr = "rotation";
   const char *scalstr = "scale";
   const char *elemstr = "element";
   const char *istpstr = "isotope";
   const char *matestr = "material";
   const char *volustr = "volume";
   const char *assestr = "assembly";
   const char *twtrstr = "twistedtrap";
   const char *cutTstr = "cutTube";
   const char *bboxstr = "box";
   const char *xtrustr = "xtru";
   const char *arb8str = "arb8";
   const char *tubestr = "tube";
   const char *conestr = "cone";
   const char *polystr = "polycone";
   const char *hypestr = "hype";
   const char *trapstr = "trap";
   const char *trdstr = "trd";
   const char *sphestr = "sphere";
   const char *orbstr = "orb";
   const char *parastr = "para";
   const char *torustr = "torus";
   const char *hedrstr = "polyhedra";
   const char *eltustr = "eltube";
   const char *subtstr = "subtraction";
   const char *uniostr = "union";
   const char *parbstr = "paraboloid";
   const char *intestr = "intersection";
   const char *reflstr = "reflectedSolid";
   const char *ellistr = "ellipsoid";
   const char *elcnstr = "elcone";
   const char *optsstr = "opticalsurface";
   const char *skinstr = "skinsurface";
   const char *bordstr = "bordersurface";
   const char *usrstr = "userinfo";
   const char *tslstr = "tessellated";
   Bool_t hasIsotopes;
   Bool_t hasIsotopesExtended;

   if ((strcmp(name, posistr)) == 0) {
      node = PosProcess(gdml, node, attr);
   } else if ((strcmp(name, rotastr)) == 0) {
      node = RotProcess(gdml, node, attr);
   } else if ((strcmp(name, scalstr)) == 0) {
      node = SclProcess(gdml, node, attr);
   } else if ((strcmp(name, setustr)) == 0) {
      node = TopProcess(gdml, node);
   } else if ((strcmp(name, consstr)) == 0) {
      node = ConProcess(gdml, node, attr);
   } else if ((strcmp(name, varistr)) == 0) {
      node = ConProcess(gdml, node, attr);
   } else if ((strcmp(name, quanstr)) == 0) {
      node = QuantityProcess(gdml, node, attr);
   } else if ((strcmp(name, matrstr)) == 0) {
      node = MatrixProcess(gdml, node, attr);
   } else if ((strcmp(name, optsstr)) == 0) {
      node = OpticalSurfaceProcess(gdml, node, attr);
   } else if ((strcmp(name, skinstr)) == 0) {
      node = SkinSurfaceProcess(gdml, node, attr);
   } else if ((strcmp(name, bordstr)) == 0) {
      node = BorderSurfaceProcess(gdml, node, attr);
   }
   //*************eleprocess********************************

   else if (((strcmp(name, "atom")) == 0) && ((strcmp(parent, elemstr)) == 0)) {
      hasIsotopes = kFALSE;
      hasIsotopesExtended = kFALSE;
      node = EleProcess(gdml, node, parentn, hasIsotopes, hasIsotopesExtended);
   } else if ((strcmp(name, elemstr) == 0) && !gdml->HasAttr(node, "Z")) {
      hasIsotopes = kTRUE;
      hasIsotopesExtended = kFALSE;
      node = EleProcess(gdml, node, parentn, hasIsotopes, hasIsotopesExtended);
   }

   else if ((strcmp(name, elemstr) == 0) && gdml->HasAttr(node, "Z")) {
      childtmp = gdml->GetChild(node);
      if ((strcmp(gdml->GetNodeName(childtmp), "fraction") == 0)) {
         hasIsotopes = kFALSE;
         hasIsotopesExtended = kTRUE;
         node = EleProcess(gdml, node, parentn, hasIsotopes, hasIsotopesExtended);
      }
   }

   //********isoprocess******************************

   else if (((strcmp(name, "atom")) == 0) && ((strcmp(parent, istpstr)) == 0)) {
      node = IsoProcess(gdml, node, parentn);
   }

   //********matprocess***********************************
   else if ((strcmp(name, matestr)) == 0 && gdml->HasAttr(node, "Z")) {
      childtmp = gdml->GetChild(node);
      //     if ((strcmp(gdml->GetNodeName(childtmp), "fraction") == 0) || (strcmp(gdml->GetNodeName(childtmp), "D") ==
      //     0)){
      // Bool_t frac = kFALSE;
      Bool_t atom = kFALSE;
      while (childtmp) {
         // frac = strcmp(gdml->GetNodeName(childtmp),"fraction")==0;
         atom = strcmp(gdml->GetNodeName(childtmp), "atom") == 0;
         gdml->ShiftToNext(childtmp);
      }
      int z = (atom) ? 1 : 0;
      node = MatProcess(gdml, node, attr, z);
   } else if ((strcmp(name, matestr)) == 0 && !gdml->HasAttr(node, "Z")) {
      int z = 0;
      node = MatProcess(gdml, node, attr, z);
   }

   //*********************************************
   else if ((strcmp(name, volustr)) == 0) {
      node = VolProcess(gdml, node);
   } else if ((strcmp(name, bboxstr)) == 0) {
      node = Box(gdml, node, attr);
   } else if ((strcmp(name, ellistr)) == 0) {
      node = Ellipsoid(gdml, node, attr);
   } else if ((strcmp(name, elcnstr)) == 0) {
      node = ElCone(gdml, node, attr);
   } else if ((strcmp(name, cutTstr)) == 0) {
      node = CutTube(gdml, node, attr);
   } else if ((strcmp(name, arb8str)) == 0) {
      node = Arb8(gdml, node, attr);
   } else if ((strcmp(name, tubestr)) == 0) {
      node = Tube(gdml, node, attr);
   } else if ((strcmp(name, conestr)) == 0) {
      node = Cone(gdml, node, attr);
   } else if ((strcmp(name, polystr)) == 0) {
      node = Polycone(gdml, node, attr);
   } else if ((strcmp(name, trapstr)) == 0) {
      node = Trap(gdml, node, attr);
   } else if ((strcmp(name, trdstr)) == 0) {
      node = Trd(gdml, node, attr);
   } else if ((strcmp(name, sphestr)) == 0) {
      node = Sphere(gdml, node, attr);
   } else if ((strcmp(name, xtrustr)) == 0) {
      node = Xtru(gdml, node, attr);
   } else if ((strcmp(name, twtrstr)) == 0) {
      node = TwistTrap(gdml, node, attr);
   } else if ((strcmp(name, hypestr)) == 0) {
      node = Hype(gdml, node, attr);
   } else if ((strcmp(name, orbstr)) == 0) {
      node = Orb(gdml, node, attr);
   } else if ((strcmp(name, parastr)) == 0) {
      node = Para(gdml, node, attr);
   } else if ((strcmp(name, torustr)) == 0) {
      node = Torus(gdml, node, attr);
   } else if ((strcmp(name, eltustr)) == 0) {
      node = ElTube(gdml, node, attr);
   } else if ((strcmp(name, hedrstr)) == 0) {
      node = Polyhedra(gdml, node, attr);
   } else if ((strcmp(name, tslstr)) == 0) {
      node = Tessellated(gdml, node, attr);
   } else if ((strcmp(name, parbstr)) == 0) {
      node = Paraboloid(gdml, node, attr);
   } else if ((strcmp(name, subtstr)) == 0) {
      node = BooSolid(gdml, node, attr, 1);
   } else if ((strcmp(name, intestr)) == 0) {
      node = BooSolid(gdml, node, attr, 2);
   } else if ((strcmp(name, uniostr)) == 0) {
      node = BooSolid(gdml, node, attr, 3);
   } else if ((strcmp(name, reflstr)) == 0) {
      node = Reflection(gdml, node, attr);
   } else if ((strcmp(name, assestr)) == 0) {
      node = AssProcess(gdml, node);
   } else if ((strcmp(name, usrstr)) == 0) {
      node = UsrProcess(gdml, node);
      // CHECK FOR TAGS NOT SUPPORTED
   } else if (((strcmp(name, "gdml")) != 0) && ((strcmp(name, "define")) != 0) && ((strcmp(name, "element")) != 0) &&
              ((strcmp(name, "materials")) != 0) && ((strcmp(name, "solids")) != 0) &&
              ((strcmp(name, "structure")) != 0) && ((strcmp(name, "zplane")) != 0) && ((strcmp(name, "first")) != 0) &&
              ((strcmp(name, "second")) != 0) && ((strcmp(name, "twoDimVertex")) != 0) &&
              ((strcmp(name, "firstposition")) != 0) && ((strcmp(name, "firstpositionref")) != 0) &&
              ((strcmp(name, "firstrotation")) != 0) && ((strcmp(name, "firstrotationref")) != 0) &&
              ((strcmp(name, "section")) != 0) && ((strcmp(name, "world")) != 0) && ((strcmp(name, "isotope")) != 0) &&
              ((strcmp(name, "triangular")) != 0) && ((strcmp(name, "quadrangular")) != 0)) {
      std::cout << "Error: Unsupported GDML Tag Used :" << name << ". Please Check Geometry/Schema." << std::endl;
   }

   // Check for Child node - if present call this funct. recursively until no more

   XMLNodePointer_t child = gdml->GetChild(node);
   while (child != 0) {
      ParseGDML(gdml, child);
      child = gdml->GetNext(child);
   }

   return fWorldName;
}

////////////////////////////////////////////////////////////////////////////////
/// Takes a string containing a mathematical expression and returns the value of
/// the expression

double TGDMLParse::Evaluate(const char *evalline)
{

   return TFormula("TFormula", evalline).Eval(0);
}

////////////////////////////////////////////////////////////////////////////////
/// When using the 'divide' process in the geometry this function
/// sets the variable 'axis' depending on what is specified.

Int_t TGDMLParse::SetAxis(const char *axisString)
{
   Int_t axis = 0;

   if ((strcmp(axisString, "kXAxis")) == 0) {
      axis = 1;
   } else if ((strcmp(axisString, "kYAxis")) == 0) {
      axis = 2;
   } else if ((strcmp(axisString, "kZAxis")) == 0) {
      axis = 3;
   } else if ((strcmp(axisString, "kRho")) == 0) {
      axis = 1;
   } else if ((strcmp(axisString, "kPhi")) == 0) {
      axis = 2;
   }

   return axis;
}

////////////////////////////////////////////////////////////////////////////////
/// This function looks thru a string for the chars '0x' next to
/// each other, when it finds this, it calls another function to strip
/// the hex address.   It does this recursively until the end of the
/// string is reached, returning a string without any hex addresses.

const char *TGDMLParse::NameShort(const char *name)
{
   static TString stripped;
   stripped = name;
   Int_t index = stripped.Index("0x");
   if (index >= 0)
      stripped = stripped(0, index);
   return stripped.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, constants can be declared.
/// when the constant keyword is found, this function is called, and the
/// name and value of the constant is stored in the "fformvec" vector as
/// a TFormula class, representing a constant function

XMLNodePointer_t TGDMLParse::ConProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name = "";
   TString value = "";
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "value") {
         value = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   // if ((strcmp(fCurrentFile, fStartFile)) != 0) {
   // name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   //}

   Double_t val = Value(value);
   fconsts[name.Data()] = val;
   gGeoManager->AddProperty(name.Data(), val);

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// Define constant expressions used.
void TGDMLParse::DefineConstants()
{
   auto def_units = gGeoManager->GetDefaultUnits();

   // Units used in TGeo. Note that they are based on cm/degree/GeV and they are different from Geant4
   fconsts["mm"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::mm : TGeant4Unit::mm;
   fconsts["millimeter"] = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::mm : TGeant4Unit::mm;
   fconsts["cm"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::cm : TGeant4Unit::cm;
   fconsts["centimeter"] = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::cm : TGeant4Unit::cm;
   fconsts["m"]          = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::m  : TGeant4Unit::m;
   fconsts["meter"]      = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::m  : TGeant4Unit::m;
   fconsts["km"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::km : TGeant4Unit::km;
   fconsts["kilometer"]  = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::km : TGeant4Unit::km;
   fconsts["rad"]        = TGeoUnit::rad;
   fconsts["radian"]     = TGeoUnit::rad;
   fconsts["deg"]        = TGeoUnit::deg;
   fconsts["degree"]     = TGeoUnit::deg;
   fconsts["pi"]         = TGeoUnit::pi;
   fconsts["twopi"]      = TGeoUnit::twopi;
   fconsts["avogadro"]   = TMath::Na();
   fconsts["gev"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::GeV : TGeant4Unit::GeV;
   fconsts["GeV"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::GeV : TGeant4Unit::GeV;
   fconsts["mev"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::MeV : TGeant4Unit::MeV;
   fconsts["MeV"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::MeV : TGeant4Unit::MeV;
   fconsts["kev"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::keV : TGeant4Unit::keV;
   fconsts["keV"]        = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::keV : TGeant4Unit::keV;
   fconsts["ev"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::eV  : TGeant4Unit::eV;
   fconsts["eV"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::eV  : TGeant4Unit::eV;
   fconsts["s"]          = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::s   : TGeant4Unit::s;
   fconsts["ms"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::ms  : TGeant4Unit::ms;
   fconsts["ns"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::ns  : TGeant4Unit::ns;
   fconsts["us"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::us  : TGeant4Unit::us;
   fconsts["kg"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::kg  : TGeant4Unit::kg;
   fconsts["g"]          = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::g   : TGeant4Unit::g;
   fconsts["mg"]         = (def_units == TGeoManager::kRootUnits) ? TGeoUnit::mg  : TGeant4Unit::mg;
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, quantities can be declared.
/// These are treated the same as constants, but the unit has to be multiplied

XMLNodePointer_t TGDMLParse::QuantityProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name = "";
   TString value = "";
   TString unit = "1.0";
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "value") {
         value = gdml->GetAttrValue(attr);
      }
      if (tempattr == "unit") {
         unit = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   fconsts[name.Data()] = GetScaleVal(unit) * Value(value);

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, matrices
/// These are referenced by other GDML tags, such as optical surfaces
XMLNodePointer_t TGDMLParse::MatrixProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name = "";
   Int_t coldim = 1;
   std::string values;
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "coldim") {
         coldim = (Int_t)Value(gdml->GetAttrValue(attr));
      }
      if (tempattr == "values") {
         values = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   // Parse the values and create the matrix
   std::stringstream valueStream(values);
   std::vector<Double_t> valueList;
   while (!valueStream.eof()) {
      std::string matrixValue;
      valueStream >> matrixValue;
      // protect against trailing '\n' and other white spaces
      if (matrixValue.empty())
         continue;
      valueList.push_back(Value(matrixValue.c_str()));
   }

   TGDMLMatrix *matrix = new TGDMLMatrix(name, valueList.size() / coldim, coldim);
   matrix->SetMatrixAsString(values.c_str());
   for (size_t i = 0; i < valueList.size(); ++i)
      matrix->Set(i / coldim, i % coldim, valueList[i]);

   gGeoManager->AddGDMLMatrix(matrix);
   fmatrices[name.Data()] = matrix;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, optical surfaces can be defined
///
XMLNodePointer_t TGDMLParse::OpticalSurfaceProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name, propname, ref;
   TGeoOpticalSurface::ESurfaceModel model = TGeoOpticalSurface::kMglisur;
   TGeoOpticalSurface::ESurfaceFinish finish = TGeoOpticalSurface::kFpolished;
   TGeoOpticalSurface::ESurfaceType type = TGeoOpticalSurface::kTdielectric_metal;
   Double_t value = 0;
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "model") {
         model = TGeoOpticalSurface::StringToModel(gdml->GetAttrValue(attr));
      }
      if (tempattr == "finish") {
         finish = TGeoOpticalSurface::StringToFinish(gdml->GetAttrValue(attr));
      }
      if (tempattr == "type") {
         type = TGeoOpticalSurface::StringToType(gdml->GetAttrValue(attr));
      }
      if (tempattr == "value") {
         value = Value(gdml->GetAttrValue(attr));
      }
      attr = gdml->GetNextAttr(attr);
   }

   TGeoOpticalSurface *surf = new TGeoOpticalSurface(name, model, finish, type, value);

   XMLNodePointer_t child = gdml->GetChild(node);
   while (child != 0) {
      attr = gdml->GetFirstAttr(child);
      if ((strcmp(gdml->GetNodeName(child), "property")) == 0) {
         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();
            if (tempattr == "name") {
               propname = gdml->GetAttrValue(attr);
            } else if (tempattr == "ref") {
               ref = gdml->GetAttrValue(attr);
               TGDMLMatrix *matrix = fmatrices[ref.Data()];
               if (!matrix)
                  Error("OpticalSurfaceProcess", "Reference matrix %s for optical surface %s not found", ref.Data(),
                        name.Data());
               surf->AddProperty(propname, ref);
            }
            attr = gdml->GetNextAttr(attr);
         }
      } // loop on child attributes
      child = gdml->GetNext(child);
   } // loop on children
   gGeoManager->AddOpticalSurface(surf);
   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// Throughout the GDML file, a unit can de specified.   Whether it be
/// angular or linear, values can be used as well as abbreviations such as
/// 'mm' or 'deg'. This function is passed the specified unit and if it is
/// found, replaces it with the appropriate value.

TString TGDMLParse::GetScale(const char *unit)
{
   TString retunit = "";

   if (strcmp(unit, "mm") == 0) {
      retunit = "0.1";
   } else if (strcmp(unit, "milimeter") == 0) {
      retunit = "0.1";
   } else if (strcmp(unit, "cm") == 0) {
      retunit = "1.0";
   } else if (strcmp(unit, "centimeter") == 0) {
      retunit = "1.0";
   } else if (strcmp(unit, "m") == 0) {
      retunit = "100.0";
   } else if (strcmp(unit, "meter") == 0) {
      retunit = "100.0";
   } else if (strcmp(unit, "km") == 0) {
      retunit = "100000.0";
   } else if (strcmp(unit, "kilometer") == 0) {
      retunit = "100000.0";
   } else if (strcmp(unit, "rad") == 0) {
      retunit = TString::Format("%.12f", TMath::RadToDeg());
   } else if (strcmp(unit, "radian") == 0) {
      retunit = TString::Format("%.12f", TMath::RadToDeg());
   } else if (strcmp(unit, "deg") == 0) {
      retunit = "1.0";
   } else if (strcmp(unit, "degree") == 0) {
      retunit = "1.0";
   } else if (strcmp(unit, "pi") == 0) {
      retunit = "pi";
   } else if (strcmp(unit, "avogadro") == 0) {
      retunit = TString::Format("%.12g", TMath::Na());
   } else {
      Fatal("GetScale", "Unit <%s> not known", unit);
      retunit = "0";
   }
   return retunit;
}

////////////////////////////////////////////////////////////////////////////////
/// Throughout the GDML file, a unit can de specified.   Whether it be
/// angular or linear, values can be used as well as abbreviations such as
/// 'mm' or 'deg'. This function is passed the specified unit and if it is
/// found, replaces it with the appropriate value.

Double_t TGDMLParse::GetScaleVal(const char *sunit)
{
   auto def_units = gGeoManager->GetDefaultUnits();
   Double_t retunit = 0.;
   TString unit(sunit);
   unit.ToLower();

   if ((unit == "mm") || (unit == "milimeter")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 0.1 : 1.0;
   } else if ((unit == "cm") || (unit == "centimeter")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 1.0 : 10.0;
   } else if ((unit == "m") || (unit == "meter")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 100.0 : 1e3;
   } else if ((unit == "km") || (unit == "kilometer")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 100000.0 : 1e6;
   } else if ((unit == "rad") || (unit == "radian")) {
     retunit = TMath::RadToDeg();
   } else if ((unit == "deg") || (unit == "degree")) {
     retunit = 1.0;
   } else if ((unit == "ev") || (unit == "electronvolt")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 0.000000001 : 1e-6;
   } else if ((unit == "kev") || (unit == "kiloelectronvolt")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 0.000001 : 1e-3;
   } else if ((unit == "mev") || (unit == "megaelectronvolt")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 0.001 : 1.0;
   } else if ((unit == "gev") || (unit == "gigaelectronvolt")) {
     retunit = (def_units == TGeoManager::kRootUnits) ? 1.0 : 1000.0;
   } else if (unit == "pi") {
     retunit = TMath::Pi();
   } else if (unit == "avogadro") {
     retunit = TMath::Na();
   } else {
     Fatal("GetScaleVal", "Unit <%s> not known", sunit);
     retunit = 0;
   }
   return retunit;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert number in string format to double value.

Double_t TGDMLParse::Value(const char *svalue) const
{
   char *end;
   double val = strtod(svalue, &end);

   // ignore white spaces.
   while (*end != 0 && isspace(*end))
      ++end;

   // Successfully parsed all the characters up to the ending NULL, so svalue
   // was a simple number.
   if (*end == 0)
      return val;

   // Otherwise we'll use TFormula to evaluate the string, having first found
   // all the GDML variable names in it and marked them with [] so that
   // TFormula will recognize them as parameters.

   std::string expanded;
   expanded.reserve(strlen(svalue) * 2);

   // Be careful about locale so we always mean the same thing by
   // "alphanumeric"
   const std::locale &loc = std::locale::classic(); // "C" locale

   // Walk through the string inserting '[' and ']' where necessary
   const char *p = svalue;
   while (*p) {
      // Find a site for a '['. Just before the first alphabetic character
      for (; *p != 0; ++p) {
         if (std::isalpha(*p, loc) || *p == '_') {
            const char *pe = p + 1;
            // Now look for the position of the following ']'. Straight before the
            // first non-alphanumeric character
            for (; *pe != 0; ++pe) {
               if (!isalnum(*pe, loc) && *pe != '_') {
                  if (*pe == '(') {
                     // The string represents a function, so no brackets needed: copy chars and advance
                     for (; p < pe; ++p)
                        expanded += *p;
                     break;
                  } else {
                     expanded += '[';
                     for (; p < pe; ++p)
                        expanded += *p;
                     expanded += ']';
                     break;
                  }
               }
            }
            if (*pe == 0) {
               expanded += '[';
               for (; p < pe; ++p)
                  expanded += *p;
               expanded += ']';
            }
         }
         expanded += *p;
      }
   } // end loop over svalue

   TFormula f("TFormula", expanded.c_str());

   // Tell the TFormula about every parameter we know about
   for (auto it : fconsts)
      f.SetParameter(it.first.c_str(), it.second);

   val = f.Eval(0);

   if (std::isnan(val) || std::isinf(val)) {
      Fatal("Value", "Got bad value %lf from string '%s'", val, svalue);
   }

   return val;
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, positions can be declared.
/// when the position keyword is found, this function is called, and the
/// name and values of the position are converted into type TGeoPosition
/// and stored in fposmap map using the name as its key. This function
/// can also be called when declaring solids.

XMLNodePointer_t TGDMLParse::PosProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString xpos = "0";
   TString ypos = "0";
   TString zpos = "0";
   TString name = "0";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x") {
         xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "y") {
         ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         zpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "unit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);
   Double_t xline = Value(xpos) * retunit;
   Double_t yline = Value(ypos) * retunit;
   Double_t zline = Value(zpos) * retunit;

   TGeoTranslation *pos = new TGeoTranslation(xline, yline, zline);

   fposmap[name.Data()] = pos;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, rotations can be declared.
/// when the rotation keyword is found, this function is called, and the
/// name and values of the rotation are converted into type TGeoRotation
/// and stored in frotmap map using the name as its key. This function
/// can also be called when declaring solids.

XMLNodePointer_t TGDMLParse::RotProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString aunit = fDefault_aunit.c_str();
   TString xpos = "0";
   TString ypos = "0";
   TString zpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x") {
         xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "y") {
         ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         zpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "unit") {
         aunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(aunit);

   Double_t xline = Value(xpos) * retunit;
   Double_t yline = Value(ypos) * retunit;
   Double_t zline = Value(zpos) * retunit;

   TGeoRotation *rot = new TGeoRotation();

   rot->RotateZ(-zline);
   rot->RotateY(-yline);
   rot->RotateX(-xline);

   frotmap[name.Data()] = rot;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the define section of the GDML file, rotations can be declared.
/// when the scale keyword is found, this function is called, and the
/// name and values of the scale are converted into type TGeoScale
/// and stored in fsclmap map using the name as its key. This function
/// can also be called when declaring solids.

XMLNodePointer_t TGDMLParse::SclProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString xpos = "0";
   TString ypos = "0";
   TString zpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x") {
         xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "y") {
         ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         zpos = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   TGeoScale *scl = new TGeoScale(Value(xpos), Value(ypos), Value(zpos));

   fsclmap[name.Data()] = scl;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the material section of the GDML file, an isotope may be declared.
/// when the isotope keyword is found, this function is called, and the
/// required parameters are taken and stored, these are then bound and
/// converted to type TGeoIsotope and stored in fisomap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::IsoProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t parentn)
{
   TString z = "0";
   TString name = "";
   TString n = "0";
   TString atom = "0";
   TString tempattr;

   // obtain attributes for the element

   XMLAttrPointer_t attr = gdml->GetFirstAttr(parentn);

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "n") {
         n = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   // get the atom value for the element

   attr = gdml->GetFirstAttr(node);

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);

      if (tempattr == "value") {
         atom = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Int_t z2 = (Int_t)Value(z);
   Int_t n2 = (Int_t)Value(n);
   Double_t atom2 = Value(atom);

   TGeoManager *mgr = gGeoManager;
   TString iso_name = NameShort(name);
   TGeoElementTable *tab = mgr->GetElementTable();
   TGeoIsotope *iso = tab->FindIsotope(iso_name);
   if (!iso) {
      iso = new TGeoIsotope(iso_name, z2, n2, atom2);
   } else if (gDebug >= 2) {
      Info("TGDMLParse", "Re-use existing isotope: %s", iso->GetName());
   }
   fisomap[name.Data()] = iso;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// When the element keyword is found, this function is called, and the
/// name and values of the element are converted into type TGeoElement and
/// stored in felemap map using the name as its key.

XMLNodePointer_t TGDMLParse::EleProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLNodePointer_t parentn,
                                        Bool_t hasIsotopes, Bool_t hasIsotopesExtended)

{
   TString z = "0";
   TString name = "";
   TString formula = "";
   TString atom = "0";
   TString tempattr;
   Int_t ncompo = 0;
   TGeoManager *mgr = gGeoManager;
   TGeoElementTable *tab = mgr->GetElementTable();
   typedef FracMap::iterator fractions;
   FracMap fracmap;

   XMLNodePointer_t child = 0;

   // obtain attributes for the element

   XMLAttrPointer_t attr = gdml->GetFirstAttr(node);

   if (hasIsotopes) {

      // Get the name of the element
      while (attr != 0) {
         tempattr = gdml->GetAttrName(attr);
         if (tempattr == "name") {
            name = gdml->GetAttrValue(attr);

            if ((strcmp(fCurrentFile, fStartFile)) != 0) {
               name = TString::Format("%s_%s", name.Data(), fCurrentFile);
            }
            break;
         }
         attr = gdml->GetNextAttr(attr);
      }
      // Get component isotopes. Loop all children.
      child = gdml->GetChild(node);
      while (child != 0) {

         // Check for fraction node name
         if ((strcmp(gdml->GetNodeName(child), "fraction")) == 0) {
            Double_t n = 0;
            TString ref = "";
            ncompo = ncompo + 1;
            attr = gdml->GetFirstAttr(child);
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();
               if (tempattr == "n") {
                  n = Value(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  ref = gdml->GetAttrValue(attr);
                  if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                     ref = TString::Format("%s_%s", ref.Data(), fCurrentFile);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            } // loop on child attributes
            fracmap[ref.Data()] = n;
         }
         child = gdml->GetNext(child);
      } // loop on children
        // Create TGeoElement - note: Object(name, title) corresponds to Element(formula, name)
      TGeoElement *ele = tab->FindElement(NameShort(name));
      // We cannot use elements with Z = 0, so we expect a user definition
      if (ele && ele->Z() == 0)
         ele = nullptr;
      if (!ele) {
         ele = new TGeoElement(NameShort(name), NameShort(name), ncompo);
         for (fractions f = fracmap.begin(); f != fracmap.end(); ++f) {
            if (fisomap.find(f->first) != fisomap.end()) {
               ele->AddIsotope((TGeoIsotope *)fisomap[f->first], f->second);
            }
         }
      } else if (gDebug >= 2) {
         Info("TGDMLParse", "Re-use existing element: %s", ele->GetName());
      }
      felemap[name.Data()] = ele;
      return child;
   } // hasisotopes end loop

   //*************************

   if (hasIsotopesExtended) {

      while (attr != 0) {
         tempattr = gdml->GetAttrName(attr);

         if (tempattr == "name") {
            name = gdml->GetAttrValue(attr);

            if ((strcmp(fCurrentFile, fStartFile)) != 0) {
               name = TString::Format("%s_%s", name.Data(), fCurrentFile);
            }
            break;
         }
         attr = gdml->GetNextAttr(attr);
      }
      // Get component isotopes. Loop all children.
      child = gdml->GetChild(node);
      while (child != 0) {

         // Check for fraction node name
         if ((strcmp(gdml->GetNodeName(child), "fraction")) == 0) {
            Double_t n = 0;
            TString ref = "";
            ncompo = ncompo + 1;
            attr = gdml->GetFirstAttr(child);
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();
               if (tempattr == "n") {
                  n = Value(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  ref = gdml->GetAttrValue(attr);
                  if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                     ref = TString::Format("%s_%s", ref.Data(), fCurrentFile);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            } // loop on child attributes
            fracmap[ref.Data()] = n;
         }
         child = gdml->GetNext(child);
      } // loop on children
        // Create TGeoElement - note: Object(name, title) corresponds to Element(formula, name)
      TGeoElement *ele = tab->FindElement(NameShort(name));
      // We cannot use elements with Z = 0, so we expect a user definition
      if (ele && ele->Z() == 0)
         ele = nullptr;
      if (!ele) {
         ele = new TGeoElement(NameShort(name), NameShort(name), ncompo);
         for (fractions f = fracmap.begin(); f != fracmap.end(); ++f) {
            if (fisomap.find(f->first) != fisomap.end()) {
               ele->AddIsotope((TGeoIsotope *)fisomap[f->first], f->second);
            }
         }
      } else if (gDebug >= 2) {
         Info("TGDMLParse", "Re-use existing element: %s", ele->GetName());
      }
      felemap[name.Data()] = ele;
      return child;
   } // hasisotopesExtended end loop

   //***************************

   attr = gdml->GetFirstAttr(parentn);
   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);

      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "formula") {
         formula = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   // get the atom value for the element

   attr = gdml->GetFirstAttr(node);

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "value") {
         atom = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Int_t z2 = (Int_t)Value(z);
   Double_t atom2 = Value(atom);
   TGeoElement *ele = tab->FindElement(formula);
   // We cannot use elements with Z = 0, so we expect a user definition
   if (ele && ele->Z() == 0)
      ele = nullptr;

   if (!ele) {
      ele = new TGeoElement(formula, NameShort(name), z2, atom2);
   } else if (gDebug >= 2) {
      Info("TGDMLParse", "Re-use existing element: %s", ele->GetName());
   }
   felemap[name.Data()] = ele;
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the materials section of the GDML file, materials can be declared.
/// when the material keyword is found, this function is called, and the
/// name and values of the material are converted into type TGeoMaterial
/// and stored in fmatmap map using the name as its key. Mixtures can also
/// be declared, and they are converted to TGeoMixture and stored in
/// fmixmap.   These mixtures and materials are then all converted into one
/// common type - TGeoMedium.   The map fmedmap is then built up of all the
/// mixtures and materials.

XMLNodePointer_t TGDMLParse::MatProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr, int z)
{
   //! Map to hold fractions while being processed
   typedef FracMap::iterator fractions;
   //  typedef FracMap::iterator i;
   FracMap fracmap;

   TGeoManager *mgr = gGeoManager;
   TGeoElementTable *tab_ele = mgr->GetElementTable();
   TList properties, constproperties;
   properties.SetOwner();
   constproperties.SetOwner();
   // We have to assume the media are monotonic increasing starting with 1
   static int medid = mgr->GetListOfMedia()->GetSize() + 1;
   XMLNodePointer_t child = gdml->GetChild(node);
   TString tempattr = "";
   Int_t ncompo = 0, mixflag = 2;
   Double_t density = 0;
   TString name = "";
   TGeoMixture *mix = 0;
   TGeoMaterial *mat = 0;
   TString tempconst = "";
   TString matname;
   Bool_t composite = kFALSE;

   if (z == 1) {
      Double_t a = 0;
      Double_t d = 0;

      name = gdml->GetAttr(node, "name");
      if ((strcmp(fCurrentFile, fStartFile)) != 0) {
         name = TString::Format("%s_%s", name.Data(), fCurrentFile);
      }

      while (child != 0) {
         attr = gdml->GetFirstAttr(child);

         if ((strcmp(gdml->GetNodeName(child), "property")) == 0) {
            TNamed *property = new TNamed();
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "name") {
                  property->SetName(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  property->SetTitle(gdml->GetAttrValue(attr));
                  TGDMLMatrix *matrix = fmatrices[property->GetTitle()];
                  if (matrix)
                     properties.Add(property);
                  else {
                     Bool_t error = 0;
                     gGeoManager->GetProperty(property->GetTitle(), &error);
                     if (error)
                        Error("MatProcess", "Reference %s for material %s not found", property->GetTitle(),
                              name.Data());
                     else
                        constproperties.Add(property);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            }
         }

         if ((strcmp(gdml->GetNodeName(child), "atom")) == 0) {
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "value") {
                  a = Value(gdml->GetAttrValue(attr));
               }
               attr = gdml->GetNextAttr(attr);
            }
         }

         if ((strcmp(gdml->GetNodeName(child), "D")) == 0) {
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "value") {
                  d = Value(gdml->GetAttrValue(attr));
               }
               attr = gdml->GetNextAttr(attr);
            }
         }
         child = gdml->GetNext(child);
      }
      // still in the is Z else...but not in the while..
      // CHECK FOR CONSTANTS
      tempconst = gdml->GetAttr(node, "Z");

      Double_t valZ = Value(tempconst);

      TString tmpname = name;
      // deal with special case - Z of vacuum is always 0
      tmpname.ToLower();
      if (tmpname == "vacuum") {
         valZ = 0;
      }
      TString mat_name = NameShort(name);
      mat = mgr->GetMaterial(mat_name);
      if (!mat) {
         mat = new TGeoMaterial(mat_name, a, valZ, d);
      } else {
         Info("TGDMLParse", "Re-use existing material: %s", mat->GetName());
      }
      if (properties.GetSize()) {
         TNamed *property;
         TIter next(&properties);
         while ((property = (TNamed *)next()))
            mat->AddProperty(property->GetName(), property->GetTitle());
      }
      if (constproperties.GetSize()) {
         TNamed *property;
         TIter next(&constproperties);
         while ((property = (TNamed *)next()))
            mat->AddConstProperty(property->GetName(), property->GetTitle());
      }
      mixflag = 0;
      // Note: Object(name, title) corresponds to Element(formula, name)
      TGeoElement *mat_ele = tab_ele->FindElement(mat_name);
      // We cannot use elements with Z = 0, so we expect a user definition
      if (mat_ele && mat_ele->Z() == 0)
         mat_ele = nullptr;

      if (!mat_ele) {
         mat_ele = new TGeoElement(mat_name, mat_name, atoi(tempconst), a);
      } else if (gDebug >= 2) {
         Info("TGDMLParse", "Re-use existing material-element: %s", mat_ele->GetName());
      }
      felemap[name.Data()] = mat_ele;
   }

   else if (z == 0) {
      while (child != 0) {
         attr = gdml->GetFirstAttr(child);

         if ((strcmp(gdml->GetNodeName(child), "property")) == 0) {
            TNamed *property = new TNamed();
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "name") {
                  property->SetName(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  property->SetTitle(gdml->GetAttrValue(attr));
                  TGDMLMatrix *matrix = fmatrices[property->GetTitle()];
                  if (matrix)
                     properties.Add(property);
                  else {
                     Bool_t error = 0;
                     gGeoManager->GetProperty(property->GetTitle(), &error);
                     if (error)
                        Error("MatProcess", "Reference %s for material %s not found", property->GetTitle(),
                              name.Data());
                     else
                        constproperties.Add(property);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            }
         }
         if ((strcmp(gdml->GetNodeName(child), "fraction")) == 0) {
            Double_t n = 0;
            TString ref = "";
            ncompo = ncompo + 1;

            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "n") {
                  n = Value(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  ref = gdml->GetAttrValue(attr);
                  if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                     ref = TString::Format("%s_%s", ref.Data(), fCurrentFile);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            }
            fracmap[ref.Data()] = n;
         }

         else if ((strcmp(gdml->GetNodeName(child), "composite")) == 0) {
            composite = kTRUE;
            Double_t n = 0;
            TString ref = "";
            ncompo = ncompo + 1;

            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();
               if (tempattr == "n") {
                  n = Value(gdml->GetAttrValue(attr));
               } else if (tempattr == "ref") {
                  ref = gdml->GetAttrValue(attr);
                  if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                     ref = TString::Format("%s_%s", ref.Data(), fCurrentFile);
                  }
               }
               attr = gdml->GetNextAttr(attr);
            }
            fracmap[ref.Data()] = n;
         } else if ((strcmp(gdml->GetNodeName(child), "D")) == 0) {
            while (attr != 0) {
               tempattr = gdml->GetAttrName(attr);
               tempattr.ToLower();

               if (tempattr == "value") {
                  density = Value(gdml->GetAttrValue(attr));
               }
               attr = gdml->GetNextAttr(attr);
            }
         }
         child = gdml->GetNext(child);
      }
      // still in the not Z else...but not in the while..

      name = gdml->GetAttr(node, "name");
      if ((strcmp(fCurrentFile, fStartFile)) != 0) {
         name = TString::Format("%s_%s", name.Data(), fCurrentFile);
      }
      // mix = new TGeoMixture(NameShort(name), 0 /*ncompo*/, density);
      mixflag = 1;
      TString mat_name = NameShort(name);
      mat = mgr->GetMaterial(mat_name);
      if (!mat) {
         mix = new TGeoMixture(mat_name, ncompo, density);
      } else if (mat->IsMixture()) {
         mix = (TGeoMixture *)mat;
         if (gDebug >= 2)
            Info("TGDMLParse", "Re-use existing material-mixture: %s", mix->GetName());
      } else {
         Fatal("TGDMLParse", "WARNING! Inconsistent material definitions between GDML and TGeoManager");
         return child;
      }
      if (properties.GetSize()) {
         TNamed *property;
         TIter next(&properties);
         while ((property = (TNamed *)next()))
            mix->AddProperty(property->GetName(), property->GetTitle());
      }
      if (constproperties.GetSize()) {
         TNamed *property;
         TIter next(&constproperties);
         while ((property = (TNamed *)next()))
            mix->AddConstProperty(property->GetName(), property->GetTitle());
      }
      Int_t natoms;
      Double_t weight;

      for (fractions f = fracmap.begin(); f != fracmap.end(); ++f) {
         matname = f->first;
         matname = NameShort(matname);

         TGeoMaterial *mattmp = (TGeoMaterial *)gGeoManager->GetListOfMaterials()->FindObject(matname);

         if (mattmp || (felemap.find(f->first) != felemap.end())) {
            if (composite) {
               natoms = (Int_t)f->second;

               mix->AddElement(felemap[f->first], natoms);

            }

            else {
               weight = f->second;
               if (mattmp) {
                  mix->AddElement(mattmp, weight);
               } else {
                  mix->AddElement(felemap[f->first], weight);
               }
            }
         }
      }
   } // end of not Z else

   medid = medid + 1;

   TGeoMedium *med = mgr->GetMedium(NameShort(name));
   if (!med) {
      if (mixflag == 1) {
         fmixmap[name.Data()] = mix;
         med = new TGeoMedium(NameShort(name), medid, mix);
      } else if (mixflag == 0) {
         fmatmap[name.Data()] = mat;
         med = new TGeoMedium(NameShort(name), medid, mat);
      }
   } else if (gDebug >= 2) {
      Info("TGDMLParse", "Re-use existing medium: %s", med->GetName());
   }
   fmedmap[name.Data()] = med;

   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the structure section of the GDML file, skin surfaces can be declared.

XMLNodePointer_t TGDMLParse::SkinSurfaceProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name, surfname, volname;
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "surfaceproperty") {
         surfname = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   XMLNodePointer_t child = gdml->GetChild(node);
   while (child != 0) {
      attr = gdml->GetFirstAttr(child);
      if ((strcmp(gdml->GetNodeName(child), "volumeref")) == 0) {
         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();
            if (tempattr == "ref") {
               volname = gdml->GetAttrValue(attr);
            }
            attr = gdml->GetNextAttr(attr);
         }
      } // loop on child attributes
      child = gdml->GetNext(child);
   } // loop on children
   TGeoOpticalSurface *surf = gGeoManager->GetOpticalSurface(surfname);
   if (!surf)
      Fatal("SkinSurfaceProcess", "Skin surface %s: referenced optical surface %s not defined", name.Data(),
            surfname.Data());
   TGeoVolume *vol = fvolmap[volname.Data()];
   TGeoSkinSurface *skin = new TGeoSkinSurface(name, surfname, surf, vol);
   gGeoManager->AddSkinSurface(skin);
   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the structure section of the GDML file, border surfaces can be declared.

XMLNodePointer_t TGDMLParse::BorderSurfaceProcess(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name, surfname, nodename[2];
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      if (tempattr == "surfaceproperty") {
         surfname = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   XMLNodePointer_t child = gdml->GetChild(node);
   Int_t inode = 0;
   while (child != 0) {
      attr = gdml->GetFirstAttr(child);
      if ((strcmp(gdml->GetNodeName(child), "physvolref")) == 0) {
         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();
            if (tempattr == "ref") {
               nodename[inode++] = gdml->GetAttrValue(attr);
            }
            attr = gdml->GetNextAttr(attr);
         }
      } // loop on child attributes
      child = gdml->GetNext(child);
   } // loop on children
   if (inode != 2)
      Fatal("BorderSurfaceProcess", "Border surface %s not referencing two nodes", name.Data());
   TGeoOpticalSurface *surf = gGeoManager->GetOpticalSurface(surfname);
   if (!surf)
      Fatal("BorderSurfaceProcess", "Border surface %s: referenced optical surface %s not defined", name.Data(),
            surfname.Data());
   TGeoNode *node1 = fpvolmap[nodename[0].Data()];
   TGeoNode *node2 = fpvolmap[nodename[1].Data()];
   if (!node1 || !node2)
      Fatal("BorderSurfaceProcess", "Border surface %s: not found nodes %s [%s] or %s [%s]", name.Data(),
            nodename[0].Data(), node1 ? "present" : "missing", nodename[1].Data(), node2 ? "present" : "missing");

   TGeoBorderSurface *border = new TGeoBorderSurface(name, surfname, surf, node1, node2);
   gGeoManager->AddBorderSurface(border);
   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the structure section of the GDML file, volumes can be declared.
/// when the volume keyword is found, this function is called, and the
/// name and values of the volume are converted into type TGeoVolume and
/// stored in fvolmap map using the name as its key. Volumes reference to
/// a solid declared higher up in the solids section of the GDML file.
/// Some volumes reference to other physical volumes to contain inside
/// that volume, declaring positions and rotations within that volume.
/// when each 'physvol' is declared, a matrix for its rotation and
/// translation is built and the 'physvol node' is added to the original
/// volume using TGeoVolume->AddNode.
/// volume division is also declared within the volume node, and once the
/// values for the division have been collected, using TGeoVolume->divide,
/// the division can be applied.

XMLNodePointer_t TGDMLParse::VolProcess(TXMLEngine *gdml, XMLNodePointer_t node)
{
   XMLAttrPointer_t attr;
   XMLNodePointer_t subchild;
   XMLNodePointer_t subsubchild;

   XMLNodePointer_t child = gdml->GetChild(node);
   TString name;
   TString solidname = "";
   TString tempattr = "";
   TGeoShape *solid = 0;
   TGeoMedium *medium = 0;
   TGeoVolume *vol = 0;
   TGeoVolume *lv = 0;
   TGeoShape *reflex = 0;
   const Double_t *parentrot = 0;
   int yesrefl = 0;
   TString reftemp = "";
   TMap *auxmap = 0;

   while (child != 0) {
      if ((strcmp(gdml->GetNodeName(child), "solidref")) == 0) {

         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fsolmap.find(reftemp.Data()) != fsolmap.end()) {
            solid = fsolmap[reftemp.Data()];
         } else if (freflectmap.find(reftemp.Data()) != freflectmap.end()) {
            solidname = reftemp;
            reflex = fsolmap[freflectmap[reftemp.Data()]];
         } else {
            printf("Solid: %s, Not Yet Defined!\n", reftemp.Data());
         }
      }

      if ((strcmp(gdml->GetNodeName(child), "materialref")) == 0) {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fmedmap.find(reftemp.Data()) != fmedmap.end()) {
            medium = fmedmap[reftemp.Data()];
         } else {
            printf("Medium: %s, Not Yet Defined!\n", gdml->GetAttr(child, "ref"));
         }
      }

      child = gdml->GetNext(child);
   }

   name = gdml->GetAttr(node, "name");

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   if (reflex == 0) {
      vol = new TGeoVolume(NameShort(name), solid, medium);
   } else {
      vol = new TGeoVolume(NameShort(name), reflex, medium);
      freflvolmap[name.Data()] = solidname;
      TGDMLRefl *parentrefl = freflsolidmap[solidname.Data()];
      parentrot = parentrefl->GetMatrix()->GetRotationMatrix();
      yesrefl = 1;
   }

   fvolmap[name.Data()] = vol;

   // PHYSVOL - run through child nodes of VOLUME again..

   child = gdml->GetChild(node);

   while (child != 0) {
      if ((strcmp(gdml->GetNodeName(child), "physvol")) == 0) {

         TString volref = "";

         TGeoTranslation *pos = 0;
         TGeoRotation *rot = 0;
         TGeoScale *scl = 0;
         TString pnodename = gdml->GetAttr(child, "name");
         TString scopynum = gdml->GetAttr(child, "copynumber");
         Int_t copynum = (scopynum.IsNull()) ? 0 : (Int_t)Value(scopynum);

         subchild = gdml->GetChild(child);

         while (subchild != 0) {
            tempattr = gdml->GetNodeName(subchild);
            tempattr.ToLower();

            if (tempattr == "volumeref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               lv = fvolmap[reftemp.Data()];
               volref = reftemp;
            } else if (tempattr == "file") {
               const char *filevol;
               const char *prevfile = fCurrentFile;

               fCurrentFile = gdml->GetAttr(subchild, "name");
               filevol = gdml->GetAttr(subchild, "volname");

               TXMLEngine *gdml2 = new TXMLEngine;
               gdml2->SetSkipComments(kTRUE);
               XMLDocPointer_t filedoc1 = gdml2->ParseFile(fCurrentFile);
               if (filedoc1 == 0) {
                  Fatal("VolProcess", "Bad filename given %s", fCurrentFile);
               }
               // take access to main node
               XMLNodePointer_t mainnode2 = gdml2->DocGetRootElement(filedoc1);
               // increase depth counter + add DOM pointer
               fFILENO = fFILENO + 1;
               fFileEngine[fFILENO] = gdml2;

               if (ffilemap.find(fCurrentFile) != ffilemap.end()) {
                  volref = ffilemap[fCurrentFile];
               } else {
                  volref = ParseGDML(gdml2, mainnode2);
                  ffilemap[fCurrentFile] = volref;
               }

               if (filevol) {
                  volref = filevol;
                  if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                     volref = TString::Format("%s_%s", volref.Data(), fCurrentFile);
                  }
               }

               fFILENO = fFILENO - 1;
               gdml = fFileEngine[fFILENO];
               fCurrentFile = prevfile;

               lv = fvolmap[volref.Data()];
               // File tree complete - Release memory before exit

               gdml->FreeDoc(filedoc1);
               delete gdml2;
            } else if (tempattr == "position") {
               attr = gdml->GetFirstAttr(subchild);
               PosProcess(gdml, subchild, attr);
               reftemp = gdml->GetAttr(subchild, "name");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               pos = fposmap[reftemp.Data()];
            } else if (tempattr == "positionref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               if (fposmap.find(reftemp.Data()) != fposmap.end())
                  pos = fposmap[reftemp.Data()];
               else
                  std::cout << "ERROR! Physvol's position " << reftemp << " not found!" << std::endl;
            } else if (tempattr == "rotation") {
               attr = gdml->GetFirstAttr(subchild);
               RotProcess(gdml, subchild, attr);
               reftemp = gdml->GetAttr(subchild, "name");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               rot = frotmap[reftemp.Data()];
            } else if (tempattr == "rotationref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               if (frotmap.find(reftemp.Data()) != frotmap.end())
                  rot = frotmap[reftemp.Data()];
               else
                  std::cout << "ERROR! Physvol's rotation " << reftemp << " not found!" << std::endl;
            } else if (tempattr == "scale") {
               attr = gdml->GetFirstAttr(subchild);
               SclProcess(gdml, subchild, attr);
               reftemp = gdml->GetAttr(subchild, "name");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               scl = fsclmap[reftemp.Data()];
            } else if (tempattr == "scaleref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               if (fsclmap.find(reftemp.Data()) != fsclmap.end())
                  scl = fsclmap[reftemp.Data()];
               else
                  std::cout << "ERROR! Physvol's scale " << reftemp << " not found!" << std::endl;
            }

            subchild = gdml->GetNext(subchild);
         }

         // ADD PHYSVOL TO GEOMETRY
         fVolID = fVolID + 1;

         TGeoHMatrix *transform = new TGeoHMatrix();

         if (pos != 0)
            transform->SetTranslation(pos->GetTranslation());
         if (rot != 0)
            transform->SetRotation(rot->GetRotationMatrix());

         if (scl != 0) { // Scaling must be added to the rotation matrix!

            Double_t scale3x3[9];
            memset(scale3x3, 0, 9 * sizeof(Double_t));
            const Double_t *diagonal = scl->GetScale();

            scale3x3[0] = diagonal[0];
            scale3x3[4] = diagonal[1];
            scale3x3[8] = diagonal[2];

            TGeoRotation scaleMatrix;
            scaleMatrix.SetMatrix(scale3x3);
            transform->Multiply(&scaleMatrix);
         }

         // BEGIN: reflectedSolid. Remove lines between if reflectedSolid will be removed from GDML!!!

         if (freflvolmap.find(volref.Data()) != freflvolmap.end()) {
            // if the volume is a reflected volume the matrix needs to be CHANGED
            TGDMLRefl *temprefl = freflsolidmap[freflvolmap[volref.Data()]];
            transform->Multiply(temprefl->GetMatrix());
         }

         if (yesrefl == 1) {
            // reflection is done per solid so that we cancel it if exists in mother volume!!!
            TGeoRotation prot;
            prot.SetMatrix(parentrot);
            transform->MultiplyLeft(&prot);
         }

         // END: reflectedSolid

         vol->AddNode(lv, copynum, transform);
         TGeoNode *lastnode = (TGeoNode *)vol->GetNodes()->Last();
         if (!pnodename.IsNull())
            lastnode->SetName(pnodename);
         fpvolmap[lastnode->GetName()] = lastnode;
      } else if ((strcmp(gdml->GetNodeName(child), "divisionvol")) == 0) {

         TString divVolref = "";
         Int_t axis = 0;
         TString number = "";
         TString width = "";
         TString offset = "";
         TString lunit = fDefault_lunit.c_str();

         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {

            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();

            if (tempattr == "axis") {
               axis = SetAxis(gdml->GetAttrValue(attr));
            } else if (tempattr == "number") {
               number = gdml->GetAttrValue(attr);
            } else if (tempattr == "width") {
               width = gdml->GetAttrValue(attr);
            } else if (tempattr == "offset") {
               offset = gdml->GetAttrValue(attr);
            } else if (tempattr == "unit") {
               lunit = gdml->GetAttrValue(attr);
            }

            attr = gdml->GetNextAttr(attr);
         }

         subchild = gdml->GetChild(child);

         while (subchild != 0) {
            tempattr = gdml->GetNodeName(subchild);
            tempattr.ToLower();

            if (tempattr == "volumeref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               divVolref = reftemp;
            }

            subchild = gdml->GetNext(subchild);
         }

         Double_t numberline = Value(number);
         Double_t retunit = GetScaleVal(lunit);
         Double_t step = Value(width) * retunit;
         Double_t offsetline = Value(offset) * retunit;

         fVolID = fVolID + 1;
         Double_t xlo, xhi;
         vol->GetShape()->GetAxisRange(axis, xlo, xhi);

         Int_t ndiv = (Int_t)numberline;
         Double_t start = xlo + offsetline;

         Int_t numed = 0;
         TGeoVolume *old = fvolmap[NameShort(reftemp)];
         if (old) {
            // We need to recreate the content of the divided volume
            old = fvolmap[NameShort(reftemp)];
            // medium id
            numed = old->GetMedium()->GetId();
         }
         TGeoVolume *divvol = vol->Divide(NameShort(reftemp), axis, ndiv, start, step, numed);
         if (!divvol) {
            Fatal("VolProcess", "Cannot divide volume %s", vol->GetName());
            return child;
         }
         if (old && old->GetNdaughters()) {
            divvol->ReplayCreation(old);
         }
         fvolmap[NameShort(reftemp)] = divvol;

      } // end of Division else if

      else if ((strcmp(gdml->GetNodeName(child), "replicavol")) == 0) {

         TString divVolref = "";
         Int_t axis = 0;
         TString number = "";
         TString width = "";
         TString offset = "";
         TString wunit = fDefault_lunit.c_str();
         TString ounit = fDefault_lunit.c_str();
         Double_t wvalue = 0;
         Double_t ovalue = 0;

         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {

            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();

            if (tempattr == "number") {
               number = gdml->GetAttrValue(attr);
            }
            attr = gdml->GetNextAttr(attr);
         }

         subchild = gdml->GetChild(child);

         while (subchild != 0) {
            tempattr = gdml->GetNodeName(subchild);
            tempattr.ToLower();

            if (tempattr == "volumeref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               divVolref = reftemp;
            }

            if (tempattr == "replicate_along_axis") {
               subsubchild = gdml->GetChild(subchild);

               while (subsubchild != 0) {
                  if ((strcmp(gdml->GetNodeName(subsubchild), "width")) == 0) {
                     attr = gdml->GetFirstAttr(subsubchild);
                     while (attr != 0) {
                        tempattr = gdml->GetAttrName(attr);
                        tempattr.ToLower();
                        if (tempattr == "value") {
                           wvalue = Value(gdml->GetAttrValue(attr));
                        } else if (tempattr == "unit") {
                           wunit = gdml->GetAttrValue(attr);
                        }

                        attr = gdml->GetNextAttr(attr);
                     }
                  } else if ((strcmp(gdml->GetNodeName(subsubchild), "offset")) == 0) {
                     attr = gdml->GetFirstAttr(subsubchild);
                     while (attr != 0) {
                        tempattr = gdml->GetAttrName(attr);
                        tempattr.ToLower();
                        if (tempattr == "value") {
                           ovalue = Value(gdml->GetAttrValue(attr));
                        } else if (tempattr == "unit") {
                           ounit = gdml->GetAttrValue(attr);
                        }
                        attr = gdml->GetNextAttr(attr);
                     }
                  } else if ((strcmp(gdml->GetNodeName(subsubchild), "direction")) == 0) {
                     attr = gdml->GetFirstAttr(subsubchild);
                     while (attr != 0) {
                        tempattr = gdml->GetAttrName(attr);
                        tempattr.ToLower();
                        if (tempattr == "x") {
                           axis = 1;
                        } else if (tempattr == "y") {
                           axis = 2;
                        } else if (tempattr == "z") {
                           axis = 3;
                        } else if (tempattr == "rho") {
                           axis = 1;
                        } else if (tempattr == "phi") {
                           axis = 2;
                        }

                        attr = gdml->GetNextAttr(attr);
                     }
                  }

                  subsubchild = gdml->GetNext(subsubchild);
               }
            }

            subchild = gdml->GetNext(subchild);
         }

         Double_t retwunit = GetScaleVal(wunit);
         Double_t retounit = GetScaleVal(ounit);

         Double_t numberline = Value(number);
         Double_t widthline = wvalue * retwunit;
         Double_t offsetline = ovalue * retounit;

         fVolID = fVolID + 1;
         Double_t xlo, xhi;
         vol->GetShape()->GetAxisRange(axis, xlo, xhi);

         Int_t ndiv = (Int_t)numberline;
         Double_t start = xlo + offsetline;

         Double_t step = widthline;
         Int_t numed = 0;
         TGeoVolume *old = fvolmap[NameShort(reftemp)];
         if (old) {
            // We need to recreate the content of the divided volume
            old = fvolmap[NameShort(reftemp)];
            // medium id
            numed = old->GetMedium()->GetId();
         }
         TGeoVolume *divvol = vol->Divide(NameShort(reftemp), axis, ndiv, start, step, numed);
         if (!divvol) {
            Fatal("VolProcess", "Cannot divide volume %s", vol->GetName());
            return child;
         }
         if (old && old->GetNdaughters()) {
            divvol->ReplayCreation(old);
         }
         fvolmap[NameShort(reftemp)] = divvol;

      } // End of replicavol
      else if (strcmp(gdml->GetNodeName(child), "auxiliary") == 0) {
         TString auxType, auxUnit, auxValue;
         if (!auxmap) {
            // printf("Auxiliary values for volume %s\n",vol->GetName());
            auxmap = new TMap();
            vol->SetUserExtension(new TGeoRCExtension(auxmap));
         }
         attr = gdml->GetFirstAttr(child);
         while (attr) {
            if (!strcmp(gdml->GetAttrName(attr), "auxtype"))
               auxType = gdml->GetAttrValue(attr);
            else if (!strcmp(gdml->GetAttrName(attr), "auxvalue"))
               auxValue = gdml->GetAttrValue(attr);
            else if (!strcmp(gdml->GetAttrName(attr), "auxunit"))
               auxUnit = gdml->GetAttrValue(attr);
            attr = gdml->GetNextAttr(attr);
         }
         if (!auxUnit.IsNull())
            auxValue = TString::Format("%s*%s", auxValue.Data(), auxUnit.Data());
         auxmap->Add(new TObjString(auxType), new TObjString(auxValue));
         // printf("  %s: %s\n", auxType.Data(), auxValue.Data());
      }

      child = gdml->GetNext(child);
   }

   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solid section of the GDML file, boolean solids can be
/// declared. when the subtraction, intersection or union   keyword
/// is found, this function is called, and the values (rotation and
/// translation) of the solid are converted into type TGeoCompositeShape
/// and stored in fsolmap map using the name as its key.
///
///  - 1 = SUBTRACTION
///  - 2 = INTERSECTION
///  - 3 = UNION

XMLNodePointer_t TGDMLParse::BooSolid(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr, int num)
{
   TString reftemp = "";
   TString tempattr = "";
   XMLNodePointer_t child = gdml->GetChild(node);

   TGeoShape *first = 0;
   TGeoShape *second = 0;

   TGeoTranslation *firstPos = new TGeoTranslation(0, 0, 0);
   TGeoTranslation *secondPos = new TGeoTranslation(0, 0, 0);

   TGeoRotation *firstRot = new TGeoRotation();
   TGeoRotation *secondRot = new TGeoRotation();

   firstRot->RotateZ(0);
   firstRot->RotateY(0);
   firstRot->RotateX(0);

   secondRot->RotateZ(0);
   secondRot->RotateY(0);
   secondRot->RotateX(0);

   TString name = gdml->GetAttr(node, "name");

   if ((strcmp(fCurrentFile, fStartFile)) != 0)
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);

   while (child != 0) {
      tempattr = gdml->GetNodeName(child);
      tempattr.ToLower();

      if (tempattr == "first") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fsolmap.find(reftemp.Data()) != fsolmap.end()) {
            first = fsolmap[reftemp.Data()];
         }
      } else if (tempattr == "second") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fsolmap.find(reftemp.Data()) != fsolmap.end()) {
            second = fsolmap[reftemp.Data()];
         }
      } else if (tempattr == "position") {
         attr = gdml->GetFirstAttr(child);
         PosProcess(gdml, child, attr);
         reftemp = gdml->GetAttr(child, "name");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         secondPos = fposmap[reftemp.Data()];
      } else if (tempattr == "positionref") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fposmap.find(reftemp.Data()) != fposmap.end()) {
            secondPos = fposmap[reftemp.Data()];
         }
      } else if (tempattr == "rotation") {
         attr = gdml->GetFirstAttr(child);
         RotProcess(gdml, child, attr);
         reftemp = gdml->GetAttr(child, "name");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         secondRot = frotmap[reftemp.Data()];
      } else if (tempattr == "rotationref") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (frotmap.find(reftemp.Data()) != frotmap.end()) {
            secondRot = frotmap[reftemp.Data()];
         }
      } else if (tempattr == "firstposition") {
         attr = gdml->GetFirstAttr(child);
         PosProcess(gdml, child, attr);
         reftemp = gdml->GetAttr(child, "name");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         firstPos = fposmap[reftemp.Data()];
      } else if (tempattr == "firstpositionref") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (fposmap.find(reftemp.Data()) != fposmap.end()) {
            firstPos = fposmap[reftemp.Data()];
         }
      } else if (tempattr == "firstrotation") {
         attr = gdml->GetFirstAttr(child);
         RotProcess(gdml, child, attr);
         reftemp = gdml->GetAttr(child, "name");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         firstRot = frotmap[reftemp.Data()];
      } else if (tempattr == "firstrotationref") {
         reftemp = gdml->GetAttr(child, "ref");
         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         if (frotmap.find(reftemp.Data()) != frotmap.end()) {
            firstRot = frotmap[reftemp.Data()];
         }
      }
      child = gdml->GetNext(child);
   }

   TGeoMatrix *firstMatrix = new TGeoCombiTrans(*firstPos, firstRot->Inverse());
   TGeoMatrix *secondMatrix = new TGeoCombiTrans(*secondPos, secondRot->Inverse());

   TGeoCompositeShape *boolean = 0;
   if (!first || !second) {
      Fatal("BooSolid", "Incomplete solid %s, missing shape components", name.Data());
      return child;
   }
   switch (num) {
   case 1:
      boolean = new TGeoCompositeShape(NameShort(name), new TGeoSubtraction(first, second, firstMatrix, secondMatrix));
      break; // SUBTRACTION
   case 2:
      boolean = new TGeoCompositeShape(NameShort(name), new TGeoIntersection(first, second, firstMatrix, secondMatrix));
      break; // INTERSECTION
   case 3:
      boolean = new TGeoCompositeShape(NameShort(name), new TGeoUnion(first, second, firstMatrix, secondMatrix));
      break; // UNION
   default: break;
   }

   fsolmap[name.Data()] = boolean;

   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// User data to be processed.

XMLNodePointer_t TGDMLParse::UsrProcess(TXMLEngine *gdml, XMLNodePointer_t node)
{
   XMLNodePointer_t child = gdml->GetChild(node);
   TString nodename, auxtype, auxtypec, auxvalue, auxvaluec, auxunit, auxunitc;
   double value = 0.;
   TGeoRegion *region;
   while (child) {
      region = nullptr;
      nodename = gdml->GetNodeName(child);
      if (nodename == "auxiliary") {
         auxtype = gdml->GetAttr(child, "auxtype");
         auxvalue = gdml->GetAttr(child, "auxvalue");
         if (auxtype == "Region") {
            auxvalue = NameShort(auxvalue);
            region = new TGeoRegion(auxvalue);
         }
      }
      XMLNodePointer_t subchild = gdml->GetChild(child);
      while (subchild) {
         auxtypec = gdml->GetAttr(subchild, "auxtype");
         auxvaluec = gdml->GetAttr(subchild, "auxvalue");
         auxunitc = gdml->GetAttr(subchild, "auxunit");
         if (auxtypec == "volume") {
            auxvaluec = NameShort(auxvaluec);
            if (region)
               region->AddVolume(auxvaluec);
         }
         if (auxtypec.Contains("cut")) {
            value = Value(auxvaluec) * GetScaleVal(auxunitc);
            if (region)
               region->AddCut(auxtypec, value);
         }
         subchild = gdml->GetNext(subchild);
      }
      if (region) {
         gGeoManager->AddRegion(region);
         // region->Print();
      }
      child = gdml->GetNext(child);
   }
   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the structure section of the GDML file, assembly volumes can be
/// declared. when the assembly keyword is found, this function is called,
/// and the name is converted into type TGeoVolumeAssembly and
/// stored in fvolmap map using the name as its key. Some assembly volumes
/// reference to other physical volumes to contain inside that assembly,
/// declaring positions and rotations within that volume. When each 'physvol'
/// is declared, a matrix for its rotation and translation is built and the
/// 'physvol node' is added to the original assembly using TGeoVolume->AddNode.

XMLNodePointer_t TGDMLParse::AssProcess(TXMLEngine *gdml, XMLNodePointer_t node)
{
   TString name = gdml->GetAttr(node, "name");
   TString reftemp = "";

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   XMLAttrPointer_t attr;
   XMLNodePointer_t subchild;
   XMLNodePointer_t child = gdml->GetChild(node);
   TString tempattr = "";
   TGeoVolume *lv = 0;
   TGeoTranslation *pos = 0;
   TGeoRotation *rot = 0;
   TGeoCombiTrans *matr;

   TGeoVolumeAssembly *assem = new TGeoVolumeAssembly(NameShort(name));

   // PHYSVOL - run through child nodes of VOLUME again..

   //   child = gdml->GetChild(node);

   while (child != 0) {
      if ((strcmp(gdml->GetNodeName(child), "physvol")) == 0) {
         TString pnodename = gdml->GetAttr(child, "name");
         TString scopynum = gdml->GetAttr(child, "copynumber");
         Int_t copynum = (scopynum.IsNull()) ? 0 : (Int_t)Value(scopynum);

         subchild = gdml->GetChild(child);
         pos = new TGeoTranslation(0, 0, 0);
         rot = new TGeoRotation();

         while (subchild != 0) {
            tempattr = gdml->GetNodeName(subchild);
            tempattr.ToLower();

            if (tempattr == "volumeref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               lv = fvolmap[reftemp.Data()];
            } else if (tempattr == "positionref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               if (fposmap.find(reftemp.Data()) != fposmap.end()) {
                  pos = fposmap[reftemp.Data()];
               }
            } else if (tempattr == "position") {
               attr = gdml->GetFirstAttr(subchild);
               PosProcess(gdml, subchild, attr);
               reftemp = gdml->GetAttr(subchild, "name");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               pos = fposmap[reftemp.Data()];
            } else if (tempattr == "rotationref") {
               reftemp = gdml->GetAttr(subchild, "ref");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               if (frotmap.find(reftemp.Data()) != frotmap.end()) {
                  rot = frotmap[reftemp.Data()];
               }
            } else if (tempattr == "rotation") {
               attr = gdml->GetFirstAttr(subchild);
               RotProcess(gdml, subchild, attr);
               reftemp = gdml->GetAttr(subchild, "name");
               if ((strcmp(fCurrentFile, fStartFile)) != 0) {
                  reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
               }
               rot = frotmap[reftemp.Data()];
            }

            subchild = gdml->GetNext(subchild);
         }

         // ADD PHYSVOL TO GEOMETRY
         fVolID = fVolID + 1;
         matr = new TGeoCombiTrans(*pos, *rot);
         assem->AddNode(lv, copynum, matr);
         TGeoNode *lastnode = (TGeoNode *)assem->GetNodes()->Last();
         if (!pnodename.IsNull())
            lastnode->SetName(pnodename);
         fpvolmap[lastnode->GetName()] = lastnode;
      }
      child = gdml->GetNext(child);
   }

   fvolmap[name.Data()] = assem;
   return child;
}

////////////////////////////////////////////////////////////////////////////////
/// In the setup section of the GDML file, the top volume need to be
/// declared. when the setup keyword is found, this function is called,
/// and the top volume ref is taken and 'world' is set

XMLNodePointer_t TGDMLParse::TopProcess(TXMLEngine *gdml, XMLNodePointer_t node)
{
   const char *name = gdml->GetAttr(node, "name");
   gGeoManager->SetName(name);
   XMLNodePointer_t child = gdml->GetChild(node);
   TString reftemp = "";

   while (child != 0) {

      if ((strcmp(gdml->GetNodeName(child), "world") == 0)) {
         // const char* reftemp;
         // TString reftemp = "";
         reftemp = gdml->GetAttr(child, "ref");

         if ((strcmp(fCurrentFile, fStartFile)) != 0) {
            reftemp = TString::Format("%s_%s", reftemp.Data(), fCurrentFile);
         }
         fWorld = fvolmap[reftemp.Data()];
         fWorldName = reftemp.Data();
      }
      child = gdml->GetNext(child);
   }
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a box may be declared.
/// when the box keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoBBox and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Box(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString xpos = "0";
   TString ypos = "0";
   TString zpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x") {
         xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "y") {
         ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         zpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t xline = 0.5 * Value(xpos) * retunit;
   Double_t yline = 0.5 * Value(ypos) * retunit;
   Double_t zline = 0.5 * Value(zpos) * retunit;

   TGeoBBox *box = new TGeoBBox(NameShort(name), xline, yline, zline);

   fsolmap[name.Data()] = box;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, an ellipsoid may be declared.
/// Unfortunately, the ellipsoid is not supported under ROOT so,
/// when the ellipsoid keyword is found, this function is called
/// to convert it to a simple box with similar dimensions, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoBBox and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Ellipsoid(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString ax = "0";
   TString by = "0";
   TString cz = "0";
   // initialization to empty string
   TString zcut1 = "";
   TString zcut2 = "";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "ax") {
         ax = gdml->GetAttrValue(attr);
      } else if (tempattr == "by") {
         by = gdml->GetAttrValue(attr);
      } else if (tempattr == "cz") {
         cz = gdml->GetAttrValue(attr);
      } else if (tempattr == "zcut1") {
         zcut1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "zcut2") {
         zcut2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t dx = Value(ax) * retunit;
   Double_t dy = Value(by) * retunit;
   Double_t radius = Value(cz) * retunit;
   Double_t sx = dx / radius;
   Double_t sy = dy / radius;
   Double_t sz = 1.;
   Double_t z1, z2;
   // Initialization of cutting
   if (zcut1 == "") {
      z1 = -radius;
   } else {
      z1 = Value(zcut1) * retunit;
   }
   if (zcut2 == "") {
      z2 = radius;
   } else {
      z2 = Value(zcut2) * retunit;
   }

   TGeoSphere *sph = new TGeoSphere(0, radius);
   TGeoScale *scl = new TGeoScale("", sx, sy, sz);
   TGeoScaledShape *shape = new TGeoScaledShape(NameShort(name), sph, scl);

   Double_t origin[3] = {0., 0., 0.};
   origin[2] = 0.5 * (z1 + z2);
   Double_t dz = 0.5 * (z2 - z1);
   TGeoBBox *pCutBox = new TGeoBBox("cutBox", dx, dy, dz, origin);
   TGeoBoolNode *pBoolNode = new TGeoIntersection(shape, pCutBox, 0, 0);
   TGeoCompositeShape *cs = new TGeoCompositeShape(NameShort(name), pBoolNode);
   fsolmap[name.Data()] = cs;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, an elliptical cone may be declared.
/// Unfortunately, the elliptical cone is not supported under ROOT so,
/// when the elcone keyword is found, this function is called
/// to convert it to a simple box with similar dimensions, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoBBox and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::ElCone(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString dx = "0";
   TString dy = "0";
   TString zmax = "0";
   TString zcut = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "dx") {
         dx = gdml->GetAttrValue(attr);
      } else if (tempattr == "dy") {
         dy = gdml->GetAttrValue(attr);
      } else if (tempattr == "zmax") {
         zmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "zcut") {
         zcut = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   // semiaxises of elliptical cone (elcone) are different then ellipsoid

   Double_t retunit = GetScaleVal(lunit);

   // dxline and dyline are without units because they are as a ration
   Double_t dxratio = Value(dx);
   Double_t dyratio = Value(dy);
   Double_t z = Value(zmax) * retunit;
   Double_t z1 = Value(zcut) * retunit;

   if (z1 <= 0) {
      Info("ElCone", "ERROR! Parameter zcut = %.12g is not set properly, elcone will not be imported.", z1);
      return node;
   }
   if (z1 > z) {
      z1 = z;
   }
   Double_t rx1 = (z + z1) * dxratio;
   Double_t ry1 = (z + z1) * dyratio;
   Double_t rx2 = (z - z1) * dxratio;
   Double_t sx = 1.;
   Double_t sy = ry1 / rx1;
   Double_t sz = 1.;

   TGeoCone *con = new TGeoCone(z1, 0, rx1, 0, rx2);
   TGeoScale *scl = new TGeoScale("", sx, sy, sz);
   TGeoScaledShape *shape = new TGeoScaledShape(NameShort(name), con, scl);

   fsolmap[name.Data()] = shape;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Paraboloid may be declared.
/// when the paraboloid keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoParaboloid and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Paraboloid(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString rlopos = "0";
   TString rhipos = "0";
   TString dzpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rlo") {
         rlopos = gdml->GetAttrValue(attr);
      } else if (tempattr == "rhi") {
         rhipos = gdml->GetAttrValue(attr);
      } else if (tempattr == "dz") {
         dzpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t rlo = Value(rlopos) * retunit;
   Double_t rhi = Value(rhipos) * retunit;
   Double_t dz = Value(dzpos) * retunit;

   TGeoParaboloid *paraboloid = new TGeoParaboloid(NameShort(name), rlo, rhi, dz);

   fsolmap[name.Data()] = paraboloid;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, an Arb8 may be declared.
/// when the arb8 keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoArb8 and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Arb8(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString v1xpos = "0";
   TString v1ypos = "0";
   TString v2xpos = "0";
   TString v2ypos = "0";
   TString v3xpos = "0";
   TString v3ypos = "0";
   TString v4xpos = "0";
   TString v4ypos = "0";
   TString v5xpos = "0";
   TString v5ypos = "0";
   TString v6xpos = "0";
   TString v6ypos = "0";
   TString v7xpos = "0";
   TString v7ypos = "0";
   TString v8xpos = "0";
   TString v8ypos = "0";
   TString dzpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "v1x") {
         v1xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v1y") {
         v1ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v2x") {
         v2xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v2y") {
         v2ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v3x") {
         v3xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v3y") {
         v3ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v4x") {
         v4xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v4y") {
         v4ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v5x") {
         v5xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v5y") {
         v5ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v6x") {
         v6xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v6y") {
         v6ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v7x") {
         v7xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v7y") {
         v7ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v8x") {
         v8xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "v8y") {
         v8ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "dz") {
         dzpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t v1x = Value(v1xpos) * retunit;
   Double_t v1y = Value(v1ypos) * retunit;
   Double_t v2x = Value(v2xpos) * retunit;
   Double_t v2y = Value(v2ypos) * retunit;
   Double_t v3x = Value(v3xpos) * retunit;
   Double_t v3y = Value(v3ypos) * retunit;
   Double_t v4x = Value(v4xpos) * retunit;
   Double_t v4y = Value(v4ypos) * retunit;
   Double_t v5x = Value(v5xpos) * retunit;
   Double_t v5y = Value(v5ypos) * retunit;
   Double_t v6x = Value(v6xpos) * retunit;
   Double_t v6y = Value(v6ypos) * retunit;
   Double_t v7x = Value(v7xpos) * retunit;
   Double_t v7y = Value(v7ypos) * retunit;
   Double_t v8x = Value(v8xpos) * retunit;
   Double_t v8y = Value(v8ypos) * retunit;
   Double_t dz = Value(dzpos) * retunit;

   TGeoArb8 *arb8 = new TGeoArb8(NameShort(name), dz);

   arb8->SetVertex(0, v1x, v1y);
   arb8->SetVertex(1, v2x, v2y);
   arb8->SetVertex(2, v3x, v3y);
   arb8->SetVertex(3, v4x, v4y);
   arb8->SetVertex(4, v5x, v5y);
   arb8->SetVertex(5, v6x, v6y);
   arb8->SetVertex(6, v7x, v7y);
   arb8->SetVertex(7, v8x, v8y);

   fsolmap[name.Data()] = arb8;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Tube may be declared.
/// when the tube keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoTubeSeg and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Tube(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString z = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin") {
         rmin = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax") {
         rmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rminline = Value(rmin) * retlunit;
   Double_t rmaxline = Value(rmax) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t startphideg = Value(startphi) * retaunit;
   Double_t deltaphideg = Value(deltaphi) * retaunit;
   Double_t endphideg = startphideg + deltaphideg;

   TGeoShape *tube = 0;
   if (deltaphideg < 360.)
      tube = new TGeoTubeSeg(NameShort(name), rminline, rmaxline, zline / 2, startphideg, endphideg);
   else
      tube = new TGeoTube(NameShort(name), rminline, rmaxline, zline / 2);
   fsolmap[name.Data()] = tube;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Cut Tube may be declared.
/// when the cutTube keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoCtub and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::CutTube(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString z = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString lowX = "0";
   TString lowY = "0";
   TString lowZ = "0";
   TString highX = "0";
   TString highY = "0";
   TString highZ = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin") {
         rmin = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax") {
         rmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "lowx") {
         lowX = gdml->GetAttrValue(attr);
      } else if (tempattr == "lowy") {
         lowY = gdml->GetAttrValue(attr);
      } else if (tempattr == "lowz") {
         lowZ = gdml->GetAttrValue(attr);
      } else if (tempattr == "highx") {
         highX = gdml->GetAttrValue(attr);
      } else if (tempattr == "highy") {
         highY = gdml->GetAttrValue(attr);
      } else if (tempattr == "highz") {
         highZ = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rminline = Value(rmin) * retlunit;
   Double_t rmaxline = Value(rmax) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t startphiline = Value(startphi) * retaunit;
   Double_t deltaphiline = Value(deltaphi) * retaunit + startphiline;
   Double_t lowXline = Value(lowX) * retlunit;
   Double_t lowYline = Value(lowY) * retlunit;
   Double_t lowZline = Value(lowZ) * retlunit;
   Double_t highXline = Value(highX) * retlunit;
   Double_t highYline = Value(highY) * retlunit;
   Double_t highZline = Value(highZ) * retlunit;

   TGeoCtub *cuttube = new TGeoCtub(NameShort(name), rminline, rmaxline, zline / 2, startphiline, deltaphiline,
                                    lowXline, lowYline, lowZline, highXline, highYline, highZline);

   fsolmap[name.Data()] = cuttube;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a cone may be declared.
/// when the cone keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoConSeg and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Cone(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin1 = "0";
   TString rmax1 = "0";
   TString rmin2 = "0";
   TString rmax2 = "0";
   TString z = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin1") {
         rmin1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax1") {
         rmax1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin2") {
         rmin2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax2") {
         rmax2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rmin1line = Value(rmin1) * retlunit;
   Double_t rmax1line = Value(rmax1) * retlunit;
   Double_t rmin2line = Value(rmin2) * retlunit;
   Double_t rmax2line = Value(rmax2) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t sphi = Value(startphi) * retaunit;
   Double_t dphi = Value(deltaphi) * retaunit;
   Double_t ephi = sphi + dphi;

   TGeoShape *cone = 0;
   if (dphi < 360.)
      cone = new TGeoConeSeg(NameShort(name), zline / 2, rmin1line, rmax1line, rmin2line, rmax2line, sphi, ephi);
   else
      cone = new TGeoCone(NameShort(name), zline / 2, rmin1line, rmax1line, rmin2line, rmax2line);

   fsolmap[name.Data()] = cone;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Trap may be declared.
/// when the trap keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoTrap and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Trap(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString x1 = "0";
   TString x2 = "0";
   TString x3 = "0";
   TString x4 = "0";
   TString y1 = "0";
   TString y2 = "0";
   TString z = "0";
   TString phi = "0";
   TString theta = "0";
   TString alpha1 = "0";
   TString alpha2 = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x1") {
         x1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x2") {
         x2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x3") {
         x3 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x4") {
         x4 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y1") {
         y1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y2") {
         y2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "phi") {
         phi = gdml->GetAttrValue(attr);
      } else if (tempattr == "theta") {
         theta = gdml->GetAttrValue(attr);
      } else if (tempattr == "alpha1") {
         alpha1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "alpha2") {
         alpha2 = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t x1line = Value(x1) * retlunit;
   Double_t x2line = Value(x2) * retlunit;
   Double_t x3line = Value(x3) * retlunit;
   Double_t x4line = Value(x4) * retlunit;
   Double_t y1line = Value(y1) * retlunit;
   Double_t y2line = Value(y2) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t philine = Value(phi) * retaunit;
   Double_t thetaline = Value(theta) * retaunit;
   Double_t alpha1line = Value(alpha1) * retaunit;
   Double_t alpha2line = Value(alpha2) * retaunit;

   TGeoTrap *trap = new TGeoTrap(NameShort(name), zline / 2, thetaline, philine, y1line / 2, x1line / 2, x2line / 2,
                                 alpha1line, y2line / 2, x3line / 2, x4line / 2, alpha2line);

   fsolmap[name.Data()] = trap;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Trd may be declared.
/// when the trd keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoTrd2 and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Trd(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString x1 = "0";
   TString x2 = "0";
   TString y1 = "0";
   TString y2 = "0";
   TString z = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x1") {
         x1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x2") {
         x2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y1") {
         y1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y2") {
         y2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);

   Double_t x1line = Value(x1) * retlunit;
   Double_t x2line = Value(x2) * retlunit;
   Double_t y1line = Value(y1) * retlunit;
   Double_t y2line = Value(y2) * retlunit;
   Double_t zline = Value(z) * retlunit;

   TGeoTrd2 *trd = new TGeoTrd2(NameShort(name), x1line / 2, x2line / 2, y1line / 2, y2line / 2, zline / 2);

   fsolmap[name.Data()] = trd;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Polycone may be declared.
/// when the polycone keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoPCon and stored in fsolmap map using the name
/// as its key. Polycone has Zplanes, planes along the z axis specifying
/// the rmin, rmax dimensions at that point along z.

XMLNodePointer_t TGDMLParse::Polycone(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString z = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   // START TO LOOK THRU CHILD (ZPLANE) NODES...

   XMLNodePointer_t child = gdml->GetChild(node);
   int numplanes = 0;

   while (child != 0) {
      numplanes = numplanes + 1;
      child = gdml->GetNext(child);
   }
   if (numplanes < 2) {
      Fatal("Polycone", "Found less than 2 planes for polycone %s", name.Data());
      return child;
   }

   int cols;
   int i;
   cols = 3;
   double **table = new double *[numplanes];
   for (i = 0; i < numplanes; i++) {
      table[i] = new double[cols];
   }

   child = gdml->GetChild(node);
   int planeno = 0;

   while (child != 0) {
      if (strcmp(gdml->GetNodeName(child), "zplane") == 0) {
         // removed original dec
         Double_t rminline = 0;
         Double_t rmaxline = 0;
         Double_t zline = 0;

         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();

            if (tempattr == "rmin") {
               rmin = gdml->GetAttrValue(attr);
               rminline = Value(rmin) * retlunit;
               table[planeno][0] = rminline;
            } else if (tempattr == "rmax") {
               rmax = gdml->GetAttrValue(attr);
               rmaxline = Value(rmax) * retlunit;
               table[planeno][1] = rmaxline;
            } else if (tempattr == "z") {
               z = gdml->GetAttrValue(attr);
               zline = Value(z) * retlunit;
               table[planeno][2] = zline;
            }
            attr = gdml->GetNextAttr(attr);
         }
      }
      planeno = planeno + 1;
      child = gdml->GetNext(child);
   }

   Double_t startphiline = Value(startphi) * retaunit;
   Double_t deltaphiline = Value(deltaphi) * retaunit;

   TGeoPcon *poly = new TGeoPcon(NameShort(name), startphiline, deltaphiline, numplanes);
   Int_t zno = 0;

   for (int j = 0; j < numplanes; j++) {
      poly->DefineSection(zno, table[j][2], table[j][0], table[j][1]);
      zno = zno + 1;
   }

   fsolmap[name.Data()] = poly;
   for (i = 0; i < numplanes; i++) {
      delete[] table[i];
   }
   delete[] table;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Polyhedra may be declared.
/// when the polyhedra keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoPgon and stored in fsolmap map using the name
/// as its key. Polycone has Zplanes, planes along the z axis specifying
/// the rmin, rmax dimensions at that point along z.

XMLNodePointer_t TGDMLParse::Polyhedra(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString z = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString numsides = "1";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "numsides") {
         numsides = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   // START TO LOOK THRU CHILD (ZPLANE) NODES...

   XMLNodePointer_t child = gdml->GetChild(node);
   int numplanes = 0;

   while (child != 0) {
      numplanes = numplanes + 1;
      child = gdml->GetNext(child);
   }
   if (numplanes < 2) {
      Fatal("Polyhedra", "Found less than 2 planes for polyhedra %s", name.Data());
      return child;
   }

   int cols;
   int i;
   cols = 3;
   double **table = new double *[numplanes];
   for (i = 0; i < numplanes; i++) {
      table[i] = new double[cols];
   }

   child = gdml->GetChild(node);
   int planeno = 0;

   while (child != 0) {
      if (strcmp(gdml->GetNodeName(child), "zplane") == 0) {

         Double_t rminline = 0;
         Double_t rmaxline = 0;
         Double_t zline = 0;
         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);
            tempattr.ToLower();

            if (tempattr == "rmin") {
               rmin = gdml->GetAttrValue(attr);
               rminline = Value(rmin) * retlunit;
               table[planeno][0] = rminline;
            } else if (tempattr == "rmax") {
               rmax = gdml->GetAttrValue(attr);
               rmaxline = Value(rmax) * retlunit;
               table[planeno][1] = rmaxline;
            } else if (tempattr == "z") {
               z = gdml->GetAttrValue(attr);
               zline = Value(z) * retlunit;
               table[planeno][2] = zline;
            }

            attr = gdml->GetNextAttr(attr);
         }
      }
      planeno = planeno + 1;
      child = gdml->GetNext(child);
   }

   Double_t startphiline = Value(startphi) * retaunit;
   Double_t deltaphiline = Value(deltaphi) * retaunit;
   Int_t numsidesline = (int)Value(numsides);

   TGeoPgon *polyg = new TGeoPgon(NameShort(name), startphiline, deltaphiline, numsidesline, numplanes);
   Int_t zno = 0;

   for (int j = 0; j < numplanes; j++) {
      polyg->DefineSection(zno, table[j][2], table[j][0], table[j][1]);
      zno = zno + 1;
   }

   fsolmap[name.Data()] = polyg;
   for (i = 0; i < numplanes; i++) {
      delete[] table[i];
   }
   delete[] table;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Sphere may be declared.
/// when the sphere keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoSphere and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Sphere(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString starttheta = "0";
   TString deltatheta = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin") {
         rmin = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax") {
         rmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "starttheta") {
         starttheta = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltatheta") {
         deltatheta = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rminline = Value(rmin) * retlunit;
   Double_t rmaxline = Value(rmax) * retlunit;
   Double_t startphiline = Value(startphi) * retaunit;
   Double_t deltaphiline = startphiline + Value(deltaphi) * retaunit;
   Double_t startthetaline = Value(starttheta) * retaunit;
   Double_t deltathetaline = startthetaline + Value(deltatheta) * retaunit;

   TGeoSphere *sphere =
      new TGeoSphere(NameShort(name), rminline, rmaxline, startthetaline, deltathetaline, startphiline, deltaphiline);

   fsolmap[name.Data()] = sphere;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Torus may be declared.
/// when the torus keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoTorus and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Torus(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString rtor = "0";
   TString startphi = "0";
   TString deltaphi = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin") {
         rmin = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax") {
         rmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "rtor") {
         rtor = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "startphi") {
         startphi = gdml->GetAttrValue(attr);
      } else if (tempattr == "deltaphi") {
         deltaphi = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rminline = Value(rmin) * retlunit;
   Double_t rmaxline = Value(rmax) * retlunit;
   Double_t rtorline = Value(rtor) * retlunit;
   Double_t startphiline = Value(startphi) * retaunit;
   Double_t deltaphiline = Value(deltaphi) * retaunit;

   TGeoTorus *torus = new TGeoTorus(NameShort(name), rtorline, rminline, rmaxline, startphiline, deltaphiline);

   fsolmap[name.Data()] = torus;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Hype may be declared.
/// when the hype keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoHype and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Hype(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString rmin = "0";
   TString rmax = "0";
   TString z = "0";
   TString inst = "0";
   TString outst = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmin") {
         rmin = gdml->GetAttrValue(attr);
      } else if (tempattr == "rmax") {
         rmax = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "inst") {
         inst = gdml->GetAttrValue(attr);
      } else if (tempattr == "outst") {
         outst = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t rminline = Value(rmin) * retlunit;
   Double_t rmaxline = Value(rmax) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t instline = Value(inst) * retaunit;
   Double_t outstline = Value(outst) * retaunit;

   TGeoHype *hype = new TGeoHype(NameShort(name), rminline, instline, rmaxline, outstline, zline / 2);

   fsolmap[name.Data()] = hype;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Para may be declared.
/// when the para keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoPara and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Para(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString x = "0";
   TString y = "0";
   TString z = "0";
   TString phi = "0";
   TString theta = "0";
   TString alpha = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x") {
         x = gdml->GetAttrValue(attr);
      } else if (tempattr == "y") {
         y = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "phi") {
         phi = gdml->GetAttrValue(attr);
      } else if (tempattr == "theta") {
         theta = gdml->GetAttrValue(attr);
      } else if (tempattr == "alpha") {
         alpha = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t xline = Value(x) * retlunit;
   Double_t yline = Value(y) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t philine = Value(phi) * retaunit;
   Double_t alphaline = Value(alpha) * retaunit;
   Double_t thetaline = Value(theta) * retaunit;

   TGeoPara *para = new TGeoPara(NameShort(name), xline / 2, yline / 2, zline / 2, alphaline, thetaline, philine);

   fsolmap[name.Data()] = para;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a TwistTrap may be declared.
/// when the twistedtrap keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoGTra and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::TwistTrap(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString aunit = fDefault_aunit.c_str();
   TString x1 = "0";
   TString x2 = "0";
   TString x3 = "0";
   TString x4 = "0";
   TString y1 = "0";
   TString y2 = "0";
   TString z = "0";
   TString phi = "0";
   TString theta = "0";
   TString alpha1 = "0";
   TString alpha2 = "0";
   TString twist = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "x1") {
         x1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x2") {
         x2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x3") {
         x3 = gdml->GetAttrValue(attr);
      } else if (tempattr == "x4") {
         x4 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y1") {
         y1 = gdml->GetAttrValue(attr);
      } else if (tempattr == "y2") {
         y2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "z") {
         z = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "aunit") {
         aunit = gdml->GetAttrValue(attr);
      } else if (tempattr == "phi") {
         phi = gdml->GetAttrValue(attr);
      } else if (tempattr == "theta") {
         theta = gdml->GetAttrValue(attr);
      } else if (tempattr == "alph") { // gdml schema knows only alph attribute
         alpha1 = gdml->GetAttrValue(attr);
         alpha2 = alpha1;
         //} else if (tempattr == "alpha2") {
         //   alpha2 = gdml->GetAttrValue(attr);
      } else if (tempattr == "phitwist") {
         twist = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);
   Double_t retaunit = GetScaleVal(aunit);

   Double_t x1line = Value(x1) * retlunit;
   Double_t x2line = Value(x2) * retlunit;
   Double_t x3line = Value(x3) * retlunit;
   Double_t x4line = Value(x4) * retlunit;
   Double_t y1line = Value(y1) * retlunit;
   Double_t y2line = Value(y2) * retlunit;
   Double_t zline = Value(z) * retlunit;
   Double_t philine = Value(phi) * retaunit;
   Double_t thetaline = Value(theta) * retaunit;
   Double_t alpha1line = Value(alpha1) * retaunit;
   Double_t alpha2line = Value(alpha2) * retaunit;
   Double_t twistline = Value(twist) * retaunit;

   TGeoGtra *twtrap = new TGeoGtra(NameShort(name), zline / 2, thetaline, philine, twistline, y1line / 2, x1line / 2,
                                   x2line / 2, alpha1line, y2line / 2, x3line / 2, x4line / 2, alpha2line);

   fsolmap[name.Data()] = twtrap;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a ElTube may be declared.
/// when the eltube keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoEltu and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::ElTube(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString xpos = "0";
   TString ypos = "0";
   TString zpos = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "dx") {
         xpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "dy") {
         ypos = gdml->GetAttrValue(attr);
      } else if (tempattr == "dz") {
         zpos = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t xline = Value(xpos) * retunit;
   Double_t yline = Value(ypos) * retunit;
   Double_t zline = Value(zpos) * retunit;

   TGeoEltu *eltu = new TGeoEltu(NameShort(name), xline, yline, zline);

   fsolmap[name.Data()] = eltu;

   return node;
}
////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, an Orb may be declared.
/// when the orb keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoSphere and stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Orb(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   TString r = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "r") {
         r = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retunit = GetScaleVal(lunit);

   Double_t rline = Value(r) * retunit;

   TGeoSphere *orb = new TGeoSphere(NameShort(name), 0, rline, 0, 180, 0, 360);

   fsolmap[name.Data()] = orb;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, an Xtru may be declared.
/// when the xtru keyword is found, this function is called, and the
/// dimensions required are taken and stored, these are then bound and
/// converted to type TGeoXtru and stored in fsolmap map using the name
/// as its key. The xtru has child nodes of either 'twoDimVertex'or
/// 'section'.   These two nodes define the real structure of the shape.
/// The twoDimVertex's define the x,y sizes of a vertice. The section links
/// the vertice to a position within the xtru.

XMLNodePointer_t TGDMLParse::Xtru(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString lunit = fDefault_lunit.c_str();
   //   TString aunit = "rad";
   TString x = "0";
   TString y = "0";
   TString zorder = "0";
   TString zpos = "0";
   TString xoff = "0";
   TString yoff = "0";
   TString scale = "0";
   TString name = "";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "lunit") {
         lunit = gdml->GetAttrValue(attr);
      }

      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   Double_t retlunit = GetScaleVal(lunit);

   // START TO LOOK THRU CHILD NODES...

   XMLNodePointer_t child = gdml->GetChild(node);
   int nosects = 0;
   int noverts = 0;

   while (child != 0) {
      tempattr = gdml->GetNodeName(child);

      if (tempattr == "twoDimVertex") {
         noverts = noverts + 1;
      } else if (tempattr == "section") {
         nosects = nosects + 1;
      }

      child = gdml->GetNext(child);
   }

   if (nosects < 2 || noverts < 3) {
      Fatal("Xtru", "Invalid number of sections/vertices found forxtru %s", name.Data());
      return child;
   }

   // Build the dynamic arrays..
   int cols;
   int i;
   double *vertx = new double[noverts];
   double *verty = new double[noverts];
   cols = 5;
   double **section = new double *[nosects];
   for (i = 0; i < nosects; i++) {
      section[i] = new double[cols];
   }

   child = gdml->GetChild(node);
   int sect = 0;
   int vert = 0;

   while (child != 0) {
      if (strcmp(gdml->GetNodeName(child), "twoDimVertex") == 0) {
         Double_t xline = 0;
         Double_t yline = 0;

         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);

            if (tempattr == "x") {
               x = gdml->GetAttrValue(attr);
               xline = Value(x) * retlunit;
               vertx[vert] = xline;
            } else if (tempattr == "y") {
               y = gdml->GetAttrValue(attr);
               yline = Value(y) * retlunit;
               verty[vert] = yline;
            }

            attr = gdml->GetNextAttr(attr);
         }

         vert = vert + 1;
      }

      else if (strcmp(gdml->GetNodeName(child), "section") == 0) {

         Double_t zposline = 0;
         Double_t xoffline = 0;
         Double_t yoffline = 0;

         attr = gdml->GetFirstAttr(child);

         while (attr != 0) {
            tempattr = gdml->GetAttrName(attr);

            if (tempattr == "zOrder") {
               zorder = gdml->GetAttrValue(attr);
               section[sect][0] = Value(zorder);
            } else if (tempattr == "zPosition") {
               zpos = gdml->GetAttrValue(attr);
               zposline = Value(zpos) * retlunit;
               section[sect][1] = zposline;
            } else if (tempattr == "xOffset") {
               xoff = gdml->GetAttrValue(attr);
               xoffline = Value(xoff) * retlunit;
               section[sect][2] = xoffline;
            } else if (tempattr == "yOffset") {
               yoff = gdml->GetAttrValue(attr);
               yoffline = Value(yoff) * retlunit;
               section[sect][3] = yoffline;
            } else if (tempattr == "scalingFactor") {
               scale = gdml->GetAttrValue(attr);
               section[sect][4] = Value(scale);
            }

            attr = gdml->GetNextAttr(attr);
         }

         sect = sect + 1;
      }
      child = gdml->GetNext(child);
   }

   TGeoXtru *xtru = new TGeoXtru(nosects);
   xtru->SetName(NameShort(name));
   xtru->DefinePolygon(vert, vertx, verty);

   for (int j = 0; j < sect; j++) {
      xtru->DefineSection((int)section[j][0], section[j][1], section[j][2], section[j][3], section[j][4]);
   }

   fsolmap[name.Data()] = xtru;
   delete[] vertx;
   delete[] verty;
   for (i = 0; i < nosects; i++) {
      delete[] section[i];
   }
   delete[] section;
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a tessellated shape may be declared.
/// When the tessellated keyword is found, this function is called, and the
/// triangular/quadrangular facets are read, creating the corresponding
/// TGeoTessellated object stored in fsolmap map using the name
/// as its key.

XMLNodePointer_t TGDMLParse::Tessellated(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   TString name, vname, type;
   TString tempattr;

   while (attr != nullptr) {
      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();
      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }

   auto tsl = new TGeoTessellated(NameShort(name));
   TGeoTranslation *pos = nullptr;
   Tessellated::Vertex_t vertices[4];

   auto SetVertex = [&](int i, TGeoTranslation *trans) {
      if (trans == nullptr)
         return;
      const double *tr = trans->GetTranslation();
      vertices[i].Set(tr[0], tr[1], tr[2]);
   };

   auto AddTriangularFacet = [&](bool relative) {
      if (relative) {
         vertices[2] += vertices[0] + vertices[1];
         vertices[1] += vertices[0];
      }
      tsl->AddFacet(vertices[0], vertices[1], vertices[2]);
   };

   auto AddQuadrangularFacet = [&](bool relative) {
      if (relative) {
         vertices[3] += vertices[0] + vertices[1] + vertices[2];
         vertices[2] += vertices[0] + vertices[1];
         vertices[1] += vertices[0];
      }
      tsl->AddFacet(vertices[0], vertices[1], vertices[2], vertices[3]);
   };

   // Get facet attributes
   XMLNodePointer_t child = gdml->GetChild(node);
   while (child != nullptr) {
      tempattr = gdml->GetNodeName(child);
      tempattr.ToLower();
      if (tempattr == "triangular") {
         attr = gdml->GetFirstAttr(child);
         bool relative = false;

         while (attr != nullptr) {
            tempattr = gdml->GetAttrName(attr);

            if (tempattr == "vertex1") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(0, pos);
            }

            else if (tempattr == "vertex2") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(1, pos);
            }

            else if (tempattr == "vertex3") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(2, pos);
            }

            else if (tempattr == "type") {
               type = gdml->GetAttrValue(attr);
               type.ToLower();
               relative = (type == "relative") ? true : false;
            }

            attr = gdml->GetNextAttr(attr);
         }
         AddTriangularFacet(relative);
      }

      else if (tempattr == "quadrangular") {
         attr = gdml->GetFirstAttr(child);
         bool relative = false;

         while (attr != nullptr) {
            tempattr = gdml->GetAttrName(attr);

            if (tempattr == "vertex1") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(0, pos);
            }

            else if (tempattr == "vertex2") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(1, pos);
            }

            else if (tempattr == "vertex3") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(2, pos);
            }

            else if (tempattr == "vertex4") {
               vname = gdml->GetAttrValue(attr);
               if (fposmap.find(vname.Data()) != fposmap.end())
                  pos = fposmap[vname.Data()];
               else
                  Fatal("Tessellated", "Vertex %s not defined", vname.Data());
               SetVertex(3, pos);
            }

            else if (tempattr == "type") {
               type = gdml->GetAttrValue(attr);
               type.ToLower();
               relative = (type == "relative") ? true : false;
            }

            attr = gdml->GetNextAttr(attr);
         }
         AddQuadrangularFacet(relative);
      }
      child = gdml->GetNext(child);
   }
   tsl->CloseShape(false);

   fsolmap[name.Data()] = tsl;

   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// In the solids section of the GDML file, a Reflected Solid may be
/// declared when the ReflectedSolid keyword is found, this function
/// is called. The rotation, position and scale for the reflection are
/// applied to a matrix that is then stored in the class object
/// TGDMLRefl.   This is then stored in the map freflsolidmap, with
/// the reflection name as a reference. also the name of the solid to
/// be reflected is stored in a map called freflectmap with the reflection
/// name as a reference.

XMLNodePointer_t TGDMLParse::Reflection(TXMLEngine *gdml, XMLNodePointer_t node, XMLAttrPointer_t attr)
{
   std::cout << "WARNING! The reflectedSolid is obsolete! Use scale transformation instead!" << std::endl;

   TString sx = "0";
   TString sy = "0";
   TString sz = "0";
   TString rx = "0";
   TString ry = "0";
   TString rz = "0";
   TString dx = "0";
   TString dy = "0";
   TString dz = "0";
   TString name = "0";
   TString solid = "0";
   TString tempattr;

   while (attr != 0) {

      tempattr = gdml->GetAttrName(attr);
      tempattr.ToLower();

      if (tempattr == "name") {
         name = gdml->GetAttrValue(attr);
      } else if (tempattr == "sx") {
         sx = gdml->GetAttrValue(attr);
      } else if (tempattr == "sy") {
         sy = gdml->GetAttrValue(attr);
      } else if (tempattr == "sz") {
         sz = gdml->GetAttrValue(attr);
      } else if (tempattr == "rx") {
         rx = gdml->GetAttrValue(attr);
      } else if (tempattr == "ry") {
         ry = gdml->GetAttrValue(attr);
      } else if (tempattr == "rz") {
         rz = gdml->GetAttrValue(attr);
      } else if (tempattr == "dx") {
         dx = gdml->GetAttrValue(attr);
      } else if (tempattr == "dy") {
         dy = gdml->GetAttrValue(attr);
      } else if (tempattr == "dz") {
         dz = gdml->GetAttrValue(attr);
      } else if (tempattr == "solid") {
         solid = gdml->GetAttrValue(attr);
      }
      attr = gdml->GetNextAttr(attr);
   }

   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      name = TString::Format("%s_%s", name.Data(), fCurrentFile);
   }
   if ((strcmp(fCurrentFile, fStartFile)) != 0) {
      solid = TString::Format("%s_%s", solid.Data(), fCurrentFile);
   }

   TGeoRotation *rot = new TGeoRotation();
   rot->RotateZ(-(Value(rz)));
   rot->RotateY(-(Value(ry)));
   rot->RotateX(-(Value(rx)));

   if (atoi(sx) == -1) {
      rot->ReflectX(kTRUE);
   }
   if (atoi(sy) == -1) {
      rot->ReflectY(kTRUE);
   }
   if (atoi(sz) == -1) {
      rot->ReflectZ(kTRUE);
   }

   TGeoCombiTrans *relf_matx = new TGeoCombiTrans(Value(dx), Value(dy), Value(dz), rot);

   TGDMLRefl *reflsol = new TGDMLRefl(NameShort(name), solid, relf_matx);
   freflsolidmap[name.Data()] = reflsol;
   freflectmap[name.Data()] = solid;

   return node;
}

/** \class TGDMLRefl
\ingroup Geometry_gdml

This class is a helper class for TGDMLParse.   It assists in the
reflection process.   This process takes a previously defined solid
and can reflect the matrix of it. This class stores the name of the
reflected solid, along with the name of the solid that is being
reflected, and finally the reflected solid's matrix.   This is then
recalled when the volume is used in the structure part of the gdml
file.

*/

ClassImp(TGDMLRefl);

////////////////////////////////////////////////////////////////////////////////
/// This constructor method stores the values brought in as params.

TGDMLRefl::TGDMLRefl(const char *name, const char *solid, TGeoMatrix *matrix)
{
   fNameS = name;
   fSolid = solid;
   fMatrix = matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// This accessor method returns the matrix.

TGeoMatrix *TGDMLRefl::GetMatrix()
{
   return fMatrix;
}
