// @(#)root/gdml:$Id$
// Author: Anton Pytel 15/9/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGDMLWrite
\ingroup Geometry_gdml

  This class contains implementation of converting ROOT's gGeoManager
geometry to GDML file. gGeoManager is the instance of TGeoManager class
containing tree of geometries creating resulting geometry. GDML is xml
based format of file mirroring the tree of geometries according to GDML
schema rules. For more information about GDML see http://gdml.web.cern.ch.
Each object in ROOT is represented by xml tag (=xml node/element) in GDML.

  This class is not needed to be instanciated. It should always be called
by gGeoManager->Export("xyz.gdml") method. Export is driven by extenstion
that is why ".gdml" is important in resulting name.

  Whenever a new ROOT geometry object is implemented or there is a change
in GDML schema this class is needed to be updated to ensure proper mapping
between ROOT objects and GDML elements.

  Current status of mapping ROOT -> GDML is implemented in method called
TGDMLWrite::ChooseObject and it contains following "map":

#### Solids:

~~~
TGeoBBox               ->           <box ... >
TGeoParaboloid         ->           <paraboloid ...>
TGeoSphere             ->           <sphere ...>
TGeoArb8               ->           <arb8 ...>
TGeoConeSeg            ->           <cone ...>
TGeoCone               ->           <cone ...>
TGeoPara               ->           <para ...>
TGeoTrap               ->           <trap ...>   or
-                      ->           <arb8 ...>
TGeoGtra               ->           <twistedtrap ...>   or
-                      ->           <trap ...>          or
-                      ->           <arb8 ...>
TGeoTrd1               ->           <trd ...>
TGeoTrd2               ->           <trd ...>
TGeoTubeSeg            ->           <tube ...>
TGeoCtub               ->           <cutTube ...>
TGeoTube               ->           <tube ...>
TGeoPcon               ->           <polycone ...>
TGeoTorus              ->           <torus ...>
TGeoPgon               ->           <polyhedra ...>
TGeoEltu               ->           <eltube ...>
TGeoHype               ->           <hype ...>
TGeoXtru               ->           <xtru ...>
TGeoTessellated        ->           <tessellated ...>
TGeoCompositeShape     ->           <union ...>            or
-                      ->           <subtraction ...>      or
-                      ->           <intersection ...>

Special cases of solids:
TGeoScaledShape        ->           <elcone ...>  if scaled TGeoCone or
-                      ->           element without scale
TGeoCompositeShape     ->           <ellipsoid ...>
-                                   intersection of:
-                                   scaled TGeoSphere and TGeoBBox
~~~

#### Materials:

~~~
TGeoIsotope            ->           <isotope ...>
TGeoElement            ->           <element ...>
TGeoMaterial           ->           <material ...>
TGeoMixture            ->           <material ...>
~~~

#### Structure

~~~
TGeoVolume             ->           <volume ...>   or
-                      ->           <assembly ...>
TGeoNode               ->           <physvol ...>
TGeoPatternFinder      ->           <divisionvol ...>
~~~

There are options that can be set to change resulting document

##### Options:

~~~
g - is set by default in gGeoManager, this option ensures compatibility
-   with Geant4. It means:
-   -> atomic number of material will be changed if <1 to 1
-   -> if polycone is set badly it will try to export it correctly
-   -> if widht * ndiv + offset is more then width of object being divided
-      (in divisions) then it will be rounded so it will not exceed or
-      if kPhi divsion then it will keep range of offset in -360 -> 0
f - if this option is set then names of volumes and solids will have
-   pointer as a suffix to ensure uniqness of names
n - if this option is set then names will not have suffix, but uniqness is
-   of names is not secured
- - if none of this two options (f,n) is set then default behaviour is so
-   that incremental suffix is added to the names.
-   (eg. TGeoBBox_0x1, TGeoBBox_0x2 ...)
~~~

#### USAGE:

~~~
gGeoManager->Export("output.gdml");
gGeoManager->Export("output.gdml","","vg"); //the same as previous just
                                            //options are set explicitly
gGeoManager->Export("output.gdml","","vgf");
gGeoManager->Export("output.gdml","","gn");
gGeoManager->Export("output.gdml","","f");
...
~~~

#### Note:
  Options discussed above are used only for TGDMLWrite class. There are
other options in the TGeoManager::Export(...) method that can be used.
See that function for details.

*/

#include "TGDMLWrite.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TXMLEngine.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoParaboloid.h"
#include "TGeoArb8.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoTrd1.h"
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
#include "TMath.h"
#include "TGeoBoolNode.h"
#include "TGeoMedium.h"
#include "TGeoElement.h"
#include "TGeoShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoOpticalSurface.h"
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <ctime>
#include <sstream>

ClassImp(TGDMLWrite);

TGDMLWrite *TGDMLWrite::fgGDMLWrite = nullptr;

namespace {
  struct MaterialExtractor  {
    std::set<TGeoMaterial*> materials;
    void operator() (const TGeoVolume* v)   {
      materials.insert(v->GetMaterial());
      for(Int_t i=0; i<v->GetNdaughters(); ++i)
        (*this)(v->GetNode(i)->GetVolume());
    }
  };
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGDMLWrite::TGDMLWrite()
   : TObject(),
     fIsotopeList(0),
     fElementList(0),
     fAccPatt(0),
     fRejShape(0),
     fNameList(0),
     fgNamingSpeed(0),
     fgG4Compatibility(0),
     fGdmlFile(0),
     fTopVolumeName(0),
     fGdmlE(0),
     fDefineNode(0),
     fMaterialsNode(0),
     fSolidsNode(0),
     fStructureNode(0),
     fVolCnt(0),
     fPhysVolCnt(0),
     fActNameErr(0),
     fSolCnt(0),
     fFltPrecision(17)   // %.17g
{
   if (fgGDMLWrite) delete fgGDMLWrite;
   fgGDMLWrite = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGDMLWrite::~TGDMLWrite()
{
   delete fIsotopeList;
   delete fElementList;
   delete fAccPatt;
   delete fRejShape;
   delete fNameList;

   fgGDMLWrite = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set convention of naming solids and volumes

void TGDMLWrite::SetNamingSpeed(ENamingType naming)
{
   fgNamingSpeed = naming;
}

////////////////////////////////////////////////////////////////////////////////
//wrapper of all main methods for extraction
void TGDMLWrite::WriteGDMLfile(TGeoManager * geomanager, const char* filename, TString option)
{
  TList* materials = geomanager->GetListOfMaterials();
  TGeoNode* node = geomanager->GetTopNode();
  if ( !node )   {
    Info("WriteGDMLfile", "Top volume does not exist!");
    return;
  }
  fTopVolumeName = "";
  WriteGDMLfile(geomanager, node, materials, filename, option);
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper to only selectively write one branch of the volume hierarchy to file
void TGDMLWrite::WriteGDMLfile(TGeoManager * geomanager, TGeoNode* node, const char* filename, TString option)
{
  TGeoVolume* volume = node->GetVolume();
  TList materials, volumes, nodes;
  MaterialExtractor extract;
  if ( !volume )   {
    Info("WriteGDMLfile", "Invalid Volume reference to extract GDML information!");
    return;
  }
  extract(volume);
  for(TGeoMaterial* m : extract.materials)
    materials.Add(m);
  fTopVolumeName = volume->GetName();
  fSurfaceList.clear();
  fVolumeList.clear();
  fNodeList.clear();
  WriteGDMLfile(geomanager, node, &materials, filename, option);
  materials.Clear("nodelete");
  volumes.Clear("nodelete");
  nodes.Clear("nodelete");
}

////////////////////////////////////////////////////////////////////////////////
/// Wrapper of all exporting methods
/// Creates blank GDML file and fills it with gGeoManager structure converted
/// to GDML structure of xml nodes

void TGDMLWrite::WriteGDMLfile(TGeoManager * geomanager,
                               TGeoNode* node,
                               TList* materialsLst,
                               const char* filename,
                               TString option)
{
   //option processing
   option.ToLower();
   if (option.Contains("g")) {
      SetG4Compatibility(kTRUE);
      Info("WriteGDMLfile", "Geant4 compatibility mode set");
   } else {
      SetG4Compatibility(kFALSE);
   }
   if (option.Contains("f")) {
      SetNamingSpeed(kfastButUglySufix);
      Info("WriteGDMLfile", "Fast naming convention with pointer suffix set");
   } else if (option.Contains("n")) {
      SetNamingSpeed(kwithoutSufixNotUniq);
      Info("WriteGDMLfile", "Naming without prefix set - be careful uniqness of name is not ensured");
   } else {
      SetNamingSpeed(kelegantButSlow);
      Info("WriteGDMLfile", "Potentially slow with incremental suffix naming convention set");
   }
   auto def_units = gGeoManager->GetDefaultUnits();
   switch (def_units) {
      case TGeoManager::kG4Units:
         fDefault_lunit = "mm";
         break;
      case TGeoManager::kRootUnits:
         fDefault_lunit = "cm";
         break;
      default: // G4 units
         fDefault_lunit = "mm";
         break;
   }

   //local variables
   Int_t outputLayout = 1;
   const char * krootNodeName = "gdml";
   const char * knsRefGeneral = "http://www.w3.org/2001/XMLSchema-instance";
   const char * knsNameGeneral = "xsi";
   const char * knsRefGdml = "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd";
   const char * knsNameGdml = "xsi:noNamespaceSchemaLocation";

   // First create engine
   fGdmlE = new TXMLEngine;
   fGdmlE->SetSkipComments(kTRUE);

   //create blank GDML file
   fGdmlFile = fGdmlE->NewDoc();

   //create root node and add it to blank GDML file
   XMLNodePointer_t rootNode = fGdmlE->NewChild(nullptr, nullptr, krootNodeName, nullptr);
   fGdmlE->DocSetRootElement(fGdmlFile, rootNode);

   //add namespaces to root node
   fGdmlE->NewNS(rootNode, knsRefGeneral, knsNameGeneral);
   fGdmlE->NewAttr(rootNode, nullptr, knsNameGdml, knsRefGdml);

   //initialize general lists and <define>, <solids>, <structure> nodes
   fIsotopeList  = new StructLst;
   fElementList  = new StructLst;

   fNameList     = new NameLst;

   fDefineNode = fGdmlE->NewChild(nullptr, nullptr, "define", nullptr);
   fSolidsNode = fGdmlE->NewChild(nullptr, nullptr, "solids", nullptr);
   fStructureNode = fGdmlE->NewChild(nullptr, nullptr, "structure", nullptr);
   //========================

   //initialize list of accepted patterns for divisions (in ExtractVolumes)
   fAccPatt   = new StructLst;
   fAccPatt->fLst["TGeoPatternX"] = kTRUE;
   fAccPatt->fLst["TGeoPatternY"] = kTRUE;
   fAccPatt->fLst["TGeoPatternZ"] = kTRUE;
   fAccPatt->fLst["TGeoPatternCylR"] = kTRUE;
   fAccPatt->fLst["TGeoPatternCylPhi"] = kTRUE;
   //========================

   //initialize list of rejected shapes for divisions (in ExtractVolumes)
   fRejShape     = new StructLst;
   //this shapes are rejected because, it is not possible to divide trd2
   //in Y axis and while only trd2 object is imported from GDML
   //it causes a problem when TGeoTrd1 is divided in Y axis
   fRejShape->fLst["TGeoTrd1"] = kTRUE;
   fRejShape->fLst["TGeoTrd2"] = kTRUE;
   //=========================

   //Initialize global counters
   fActNameErr = 0;
   fVolCnt = 0;
   fPhysVolCnt = 0;
   fSolCnt = 0;

   //calling main extraction functions (with measuring time)
   time_t startT, endT;
   startT = time(nullptr);
   ExtractMatrices(geomanager->GetListOfGDMLMatrices());
   ExtractConstants(geomanager);
   fMaterialsNode = ExtractMaterials(materialsLst);

   Info("WriteGDMLfile", "Extracting volumes");
   ExtractVolumes(node);
   Info("WriteGDMLfile", "%i solids added", fSolCnt);
   Info("WriteGDMLfile", "%i volumes added", fVolCnt);
   Info("WriteGDMLfile", "%i physvolumes added", fPhysVolCnt);
   ExtractSkinSurfaces(geomanager->GetListOfSkinSurfaces());
   ExtractBorderSurfaces(geomanager->GetListOfBorderSurfaces());
   ExtractOpticalSurfaces(geomanager->GetListOfOpticalSurfaces());
   endT = time(nullptr);
   //<gdml>
   fGdmlE->AddChild(rootNode, fDefineNode);                 //  <define>...</define>
   fGdmlE->AddChild(rootNode, fMaterialsNode);              //  <materials>...</materials>
   fGdmlE->AddChild(rootNode, fSolidsNode);                 //  <solids>...</solids>
   fGdmlE->AddChild(rootNode, fStructureNode);              //  <structure>...</structure>
   fGdmlE->AddChild(rootNode, CreateSetupN(fTopVolumeName.Data())); //  <setup>...</setup>
   //</gdml>
   Double_t tdiffI = difftime(endT, startT);
   TString tdiffS = (tdiffI == 0 ? TString("< 1 s") : TString::Format("%.0lf s", tdiffI));
   Info("WriteGDMLfile", "Exporting time: %s", tdiffS.Data());
   //=========================

   //Saving document
   fGdmlE->SaveDoc(fGdmlFile, filename, outputLayout);
   Info("WriteGDMLfile", "File %s saved", filename);
   //cleaning
   fGdmlE->FreeDoc(fGdmlFile);
   //unset processing bits:
   UnsetTemporaryBits(geomanager);
   delete fGdmlE;
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting GDML matrices

void TGDMLWrite::ExtractMatrices(TObjArray* matrixList)
{
   if (!matrixList->GetEntriesFast()) return;
   XMLNodePointer_t matrixN;
   TIter next(matrixList);
   TGDMLMatrix *matrix;
   while ((matrix = (TGDMLMatrix*)next())) {
     matrixN = CreateMatrixN(matrix);
     fGdmlE->AddChild(fDefineNode, matrixN);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting GDML matrices

void TGDMLWrite::ExtractConstants(TGeoManager *geom)
{
   if (!geom->GetNproperties()) return;
   XMLNodePointer_t constantN;
   TString property;
   Double_t value;
   for (Int_t i = 0; i < geom->GetNproperties(); ++i) {
      value = geom->GetProperty(i, property);
      constantN = CreateConstantN(property.Data(), value);
      fGdmlE->AddChild(fDefineNode, constantN);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting optical surfaces

void TGDMLWrite::ExtractOpticalSurfaces(TObjArray *surfaces)
{
   if (!surfaces->GetEntriesFast()) return;
   XMLNodePointer_t surfaceN;
   TIter next(surfaces);
   TGeoOpticalSurface *surf;
   while ((surf = (TGeoOpticalSurface*)next())) {
      if ( fSurfaceList.find(surf) == fSurfaceList.end() ) continue;
      surfaceN = CreateOpticalSurfaceN(surf);
      fGdmlE->AddChild(fSolidsNode, surfaceN);
      // Info("ExtractSkinSurfaces", "Extracted optical surface: %s",surf->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting skin surfaces

void TGDMLWrite::ExtractSkinSurfaces(TObjArray *surfaces)
{
   if (!surfaces->GetEntriesFast()) return;
   XMLNodePointer_t surfaceN;
   TIter next(surfaces);
   TGeoSkinSurface *surf;
   while ((surf = (TGeoSkinSurface*)next())) {
      if ( fVolumeList.find(surf->GetVolume()) == fVolumeList.end() ) continue;
      surfaceN = CreateSkinSurfaceN(surf);
      fGdmlE->AddChild(fStructureNode, surfaceN);
      fSurfaceList.insert(surf->GetSurface());
      // Info("ExtractSkinSurfaces", "Extracted skin surface: %s",surf->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting border surfaces

void TGDMLWrite::ExtractBorderSurfaces(TObjArray *surfaces)
{
   if (!surfaces->GetEntriesFast()) return;
   XMLNodePointer_t surfaceN;
   TIter next(surfaces);
   TGeoBorderSurface *surf;
   while ((surf = (TGeoBorderSurface*)next())) {
      auto ia = fNodeList.find(surf->GetNode1());
      auto ib = fNodeList.find(surf->GetNode2());
      if ( ia == fNodeList.end() && ib == fNodeList.end() )  {
        continue;
      }
      else if ( ia == fNodeList.end() && ib != fNodeList.end() )  {
        Warning("ExtractBorderSurfaces", "Inconsistent border surface extraction %s: Node %s"
             " is not part of GDML!",surf->GetName(), surf->GetNode1()->GetName());
        continue;
      }
      else if ( ia != fNodeList.end() && ib == fNodeList.end() )  {
        Warning("ExtractBorderSurfaces", "Inconsistent border surface extraction %s: Node %s"
             " is not part of GDML!",surf->GetName(), surf->GetNode2()->GetName());
        continue;
      }
      surfaceN = CreateBorderSurfaceN(surf);
      fGdmlE->AddChild(fStructureNode, surfaceN);
      fSurfaceList.insert(surf->GetSurface());
      // Info("ExtractBorderSurfaces", "Extracted border surface: %s",surf->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method exporting materials

XMLNodePointer_t TGDMLWrite::ExtractMaterials(TList* materialsLst)
{
   Info("ExtractMaterials", "Extracting materials");
   //crate main <materials> node
   XMLNodePointer_t materialsN = fGdmlE->NewChild(nullptr, nullptr, "materials", nullptr);
   Int_t matcnt = 0;

   //go through materials  - iterator and object declaration
   TIter next(materialsLst);
   TGeoMaterial *lmaterial;

   while ((lmaterial = (TGeoMaterial *)next())) {
      //generate uniq name
      TString lname = GenName(lmaterial->GetName(), TString::Format("%p", lmaterial));

      if (lmaterial->IsMixture()) {
         TGeoMixture  *lmixture = (TGeoMixture *)lmaterial;
         XMLNodePointer_t mixtureN = CreateMixtureN(lmixture, materialsN, lname);
         fGdmlE->AddChild(materialsN, mixtureN);
      } else {
         XMLNodePointer_t materialN = CreateMaterialN(lmaterial, lname);
         fGdmlE->AddChild(materialsN, materialN);
      }
      matcnt++;
   }
   Info("ExtractMaterials", "%i materials added", matcnt);
   return materialsN;
}

////////////////////////////////////////////////////////////////////////////////
/// Method creating solid to xml file and returning its name

TString TGDMLWrite::ExtractSolid(TGeoShape* volShape)
{
   XMLNodePointer_t solidN;
   TString solname = "";
   solidN = ChooseObject(volShape);  //volume->GetShape()
   fGdmlE->AddChild(fSolidsNode, solidN);
   if (solidN != nullptr) fSolCnt++;
   solname = fNameList->fLst[TString::Format("%p", volShape)];
   if (solname.Contains("missing_")) {
      solname = "-1";
   }
   return solname;
}


////////////////////////////////////////////////////////////////////////////////
/// Method extracting geometry structure recursively

void TGDMLWrite::ExtractVolumes(TGeoNode* node)
{
   XMLNodePointer_t volumeN, childN;
   TGeoVolume * volume = node->GetVolume();
   TString volname, matname, solname, pattClsName, nodeVolNameBak;
   TGeoPatternFinder *pattFinder = 0;
   Bool_t isPattern = kFALSE;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   
   fNodeList.insert(node);
   fVolumeList.insert(volume);
   //create the name for volume/assembly
   if (volume->IsTopVolume()) {
      //not needed a special function for generating name
      volname = volume->GetName();
      fTopVolumeName = volname;
      //register name to the pointer
      fNameList->fLst[TString::Format("%p", volume)] = volname;
   } else {
      volname = GenName(volume->GetName(), TString::Format("%p", volume));
   }
   
   //start to create main volume/assembly node
   if (volume->IsAssembly()) {
      volumeN = StartAssemblyN(volname);
   } else {
      //get reference material and add solid to <solids> + get name
      matname = fNameList->fLst[TString::Format("%p", volume->GetMaterial())];
      solname = ExtractSolid(volume->GetShape());
      //If solid is not supported or corrupted
      if (solname == "-1") {
         Info("ExtractVolumes", "ERROR! %s volume was not added, because solid is either not supported or corrupted",
              volname.Data());
         //set volume as missing volume
         fNameList->fLst[TString::Format("%p", volume)] = "missing_" + volname;
         return;
      }
      volumeN = StartVolumeN(volname, solname, matname);

      //divisionvol can't be in assembly
      pattFinder = volume->GetFinder();
      //if found pattern
      if (pattFinder) {
         pattClsName = TString::Format("%s", pattFinder->ClassName());
         TString shapeCls = TString::Format("%s", volume->GetShape()->ClassName());
         //if pattern in accepted pattern list and not in shape rejected list
         if ((fAccPatt->fLst[pattClsName] == kTRUE) &&
             (fRejShape->fLst[shapeCls] != kTRUE)) {
            isPattern = kTRUE;
         }
      }
   }
   //get all nodes in volume
   TObjArray *nodeLst = volume->GetNodes();
   TIter next(nodeLst);
   TGeoNode *geoNode;
   Int_t nCnt = 0;
   //loop through all nodes
   while ((geoNode = (TGeoNode *) next())) {
      //get volume of current node and if not processed then process it
      TGeoVolume * subvol = geoNode->GetVolume();
      fNodeList.insert(geoNode);
      if (subvol->TestAttBit(fgkProcBitVol) == kFALSE) {
         subvol->SetAttBit(fgkProcBitVol);
         ExtractVolumes(geoNode);
      }

      //volume of this node has to exist because it was processed recursively
      TString nodevolname = fNameList->fLst[TString::Format("%p", geoNode->GetVolume())];
      if (nodevolname.Contains("missing_")) {
         continue;
      }
      if (nCnt == 0) { //save name of the first node for divisionvol
         nodeVolNameBak = nodevolname;
      }

      if (isPattern == kFALSE) {
         //create name for node
         TString nodename, posname, rotname;
         nodename = GenName(geoNode->GetName(), TString::Format("%p", geoNode));
         nodename = nodename + "in" + volname;

         //create name for position and clear rotation
         posname = nodename + "pos";
         rotname = "";

         //position
         const Double_t * pos = geoNode->GetMatrix()->GetTranslation();
         Xyz nodPos;
         nodPos.x = pos[0];
         nodPos.y = pos[1];
         nodPos.z = pos[2];
         childN = CreatePositionN(posname.Data(), nodPos);
         fGdmlE->AddChild(fDefineNode, childN); //adding node to <define> node
         //Deal with reflection
         XMLNodePointer_t scaleN = nullptr;
         Double_t lx, ly, lz;
         Double_t xangle = 0;
         Double_t zangle = 0;
         lx = geoNode->GetMatrix()->GetRotationMatrix()[0];
         ly = geoNode->GetMatrix()->GetRotationMatrix()[4];
         lz = geoNode->GetMatrix()->GetRotationMatrix()[8];
         if (geoNode->GetMatrix()->IsReflection()
             && TMath::Abs(lx) == 1 &&  TMath::Abs(ly) == 1 && TMath::Abs(lz) == 1) {
            scaleN = fGdmlE->NewChild(nullptr, nullptr, "scale", nullptr);
            fGdmlE->NewAttr(scaleN, nullptr, "name", (nodename + "scl").Data());
            fGdmlE->NewAttr(scaleN, nullptr, "x", TString::Format(fltPrecision.Data(), lx));
            fGdmlE->NewAttr(scaleN, nullptr, "y", TString::Format(fltPrecision.Data(), ly));
            fGdmlE->NewAttr(scaleN, nullptr, "z", TString::Format(fltPrecision.Data(), lz));
            //experimentally found out, that rotation should be updated like this
            if (lx == -1) {
               zangle = 180;
            }
            if (lz == -1) {
               xangle = 180;
            }
         }

         //rotation
         TGDMLWrite::Xyz lxyz = GetXYZangles(geoNode->GetMatrix()->GetRotationMatrix());
         lxyz.x -= xangle;
         lxyz.z -= zangle;
         if ((lxyz.x != 0.0) || (lxyz.y != 0.0) || (lxyz.z != 0.0)) {
            rotname = nodename + "rot";
            childN = CreateRotationN(rotname.Data(), lxyz);
            fGdmlE->AddChild(fDefineNode, childN); //adding node to <define> node
         }

         //create physvol for main volume/assembly node
         childN = CreatePhysVolN(geoNode->GetName(), geoNode->GetNumber(), nodevolname.Data(), posname.Data(), rotname.Data(), scaleN);
         fGdmlE->AddChild(volumeN, childN);
      }
      nCnt++;
   }
   //create only one divisionvol node
   if (isPattern && pattFinder) {
      //retrieve attributes of division
      Int_t ndiv, divaxis;
      Double_t offset, width, xlo, xhi;
      TString axis, unit;

      ndiv = pattFinder->GetNdiv();
      width = pattFinder->GetStep();

      divaxis = pattFinder->GetDivAxis();
      volume->GetShape()->GetAxisRange(divaxis, xlo, xhi);

      //compute relative start (not positional)
      offset = pattFinder->GetStart() - xlo;
      axis = GetPattAxis(divaxis, pattClsName, unit);

      //create division node
      childN = CreateDivisionN(offset, width, ndiv, axis.Data(), unit.Data(), nodeVolNameBak.Data());
      fGdmlE->AddChild(volumeN, childN);
   }

   fVolCnt++;
   //add volume/assembly node into the <structure> node
   fGdmlE->AddChild(fStructureNode, volumeN);

}

////////////////////////////////////////////////////////////////////////////////
/// Creates "atom" node for GDML

XMLNodePointer_t TGDMLWrite::CreateAtomN(Double_t atom, const char * unit)
{
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   XMLNodePointer_t atomN = fGdmlE->NewChild(nullptr, nullptr, "atom", nullptr);
   if ( gGeoManager->GetDefaultUnits() != TGeoManager::kRootUnits ) atom /= 1e19; // Correct for G4 unit system
   fGdmlE->NewAttr(atomN, nullptr, "unit", unit);
   fGdmlE->NewAttr(atomN, nullptr, "value", TString::Format(fltPrecision.Data(), atom));
   return atomN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "D" density node for GDML

XMLNodePointer_t TGDMLWrite::CreateDN(Double_t density, const char * unit)
{
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   XMLNodePointer_t densN = fGdmlE->NewChild(nullptr, nullptr, "D", nullptr);
   if ( gGeoManager->GetDefaultUnits() != TGeoManager::kRootUnits ) density /= 1e16; // Correct for G4 unit system
   fGdmlE->NewAttr(densN, nullptr, "unit", unit);
   fGdmlE->NewAttr(densN, nullptr, "value", TString::Format(fltPrecision.Data(), density));
   return densN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "fraction" node for GDML

XMLNodePointer_t TGDMLWrite::CreateFractionN(Double_t percentage, const char * refName)
{
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   XMLNodePointer_t fractN = fGdmlE->NewChild(nullptr, nullptr, "fraction", nullptr);
   fGdmlE->NewAttr(fractN, nullptr, "n", TString::Format(fltPrecision.Data(), percentage));
   fGdmlE->NewAttr(fractN, nullptr, "ref", refName);
   return fractN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "property" node for GDML

XMLNodePointer_t TGDMLWrite::CreatePropertyN(TNamed const &property)
{
  XMLNodePointer_t propertyN = fGdmlE->NewChild(nullptr, nullptr, "property", nullptr);
  fGdmlE->NewAttr(propertyN, nullptr, "name", property.GetName());
  fGdmlE->NewAttr(propertyN, nullptr, "ref", property.GetTitle());
  return propertyN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "isotope" node for GDML

XMLNodePointer_t TGDMLWrite::CreateIsotopN(TGeoIsotope * isotope, const char * name)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "isotope", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   fGdmlE->NewAttr(mainN, nullptr, "N", TString::Format("%i", isotope->GetN()));
   fGdmlE->NewAttr(mainN, nullptr, "Z", TString::Format("%i", isotope->GetZ()));
   fGdmlE->AddChild(mainN, CreateAtomN(isotope->GetA()));
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "element" node for GDML
///element node and attribute

XMLNodePointer_t TGDMLWrite::CreateElementN(TGeoElement * element, XMLNodePointer_t materials, const char * name)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "element", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   //local associative arrays for saving isotopes and their weight
   //inside element
   NameListF wPercentage;
   NameListI wCounter;

   if (element->HasIsotopes()) {
      Int_t nOfIso = element->GetNisotopes();
      //go through isotopes
      for (Int_t idx = 0; idx < nOfIso; idx++) {
         TGeoIsotope *myIsotope = element->GetIsotope(idx);
         if (!myIsotope) {
            Fatal("CreateElementN", "Missing isotopes for element %s", element->GetName());
            return mainN;
         }

         //Get name of the Isotope (
         TString lname = myIsotope->GetName();
         //_iso suffix is added to avoid problems with same names
         //for material, element and isotopes
         lname = TString::Format("%s_iso", lname.Data());

         //cumulates abundance, in case 2 isotopes with same names
         //within one element
         wPercentage[lname] += element->GetRelativeAbundance(idx);
         wCounter[lname]++;

         //check whether isotope name is not in list of isotopes
         if (IsInList(fIsotopeList->fLst, lname)) {
            continue;
         }
         //add isotope to list of isotopes and to main <materials> node
         fIsotopeList->fLst[lname] = kTRUE;
         XMLNodePointer_t isoNode = CreateIsotopN(myIsotope, lname);
         fGdmlE->AddChild(materials, isoNode);
      }
      //loop through asoc array of isotopes
      for (NameListI::iterator itr = wCounter.begin(); itr != wCounter.end(); ++itr) {
         if (itr->second > 1) {
            Info("CreateMixtureN", "WARNING! 2 equal isotopes in one element. Check: %s isotope of %s element",
                 itr->first.Data(), name);
         }
         //add fraction child to element with reference to isotope
         fGdmlE->AddChild(mainN, CreateFractionN(wPercentage[itr->first], itr->first.Data()));
      }
   } else {
      fGdmlE->NewAttr(mainN, nullptr, "formula", element->GetName());
      Int_t valZ = element->Z();
      // Z can't be <1 in Geant4 and Z is optional parameter
      if (valZ >= 1) {
         fGdmlE->NewAttr(mainN, nullptr, "Z", TString::Format("%i", valZ));
      }
      fGdmlE->AddChild(mainN, CreateAtomN(element->A()));
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "material" node for GDML with references to other sub elements

XMLNodePointer_t TGDMLWrite::CreateMixtureN(TGeoMixture * mixture, XMLNodePointer_t materials, TString mname)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "material", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", mname);
   fGdmlE->AddChild(mainN, CreateDN(mixture->GetDensity()));
   //local associative arrays for saving elements and their weight
   //inside mixture
   NameListF wPercentage;
   NameListI wCounter;

   Int_t nOfElm = mixture->GetNelements();
   //go through elements
   for (Int_t idx = 0; idx < nOfElm; idx++) {
      TGeoElement *myElement = mixture->GetElement(idx);

      //Get name of the element
      //NOTE: that for element - GetTitle() returns the "name" tag
      //and GetName() returns "formula" tag (see createElementN)
      TString lname = myElement->GetTitle();
      //_elm suffix is added to avoid problems with same names
      //for material and element
      lname = TString::Format("%s_elm", lname.Data());

      //cumulates percentage, in case 2 elements with same names within one mixture
      wPercentage[lname] += mixture->GetWmixt()[idx];
      wCounter[lname]++;

      //check whether element name is not in list of elements already created
      if (IsInList(fElementList->fLst, lname)) {
         continue;
      }

      //add element to list of elements and to main <materials> node
      fElementList->fLst[lname] = kTRUE;
      XMLNodePointer_t elmNode = CreateElementN(myElement, materials, lname);
      fGdmlE->AddChild(materials, elmNode);
   }
   //loop through asoc array
   for (NameListI::iterator itr = wCounter.begin(); itr != wCounter.end(); ++itr) {
      if (itr->second > 1) {
         Info("CreateMixtureN", "WARNING! 2 equal elements in one material. Check: %s element of %s material",
              itr->first.Data(), mname.Data());
      }
      //add fraction child to material with reference to element
      fGdmlE->AddChild(mainN, CreateFractionN(wPercentage[itr->first], itr->first.Data()));
   }

   // Write properties
   TList const &properties = mixture->GetProperties();
   if (properties.GetSize()) {
      TIter next(&properties);
      TNamed *property;
      while ((property = (TNamed*)next()))
        fGdmlE->AddChild(mainN, CreatePropertyN(*property));
   }
   // Write CONST properties
   TList const &const_properties = mixture->GetConstProperties();
   if (const_properties.GetSize()) {
      TIter next(&const_properties);
      TNamed *property;
      while ((property = (TNamed*)next()))
        fGdmlE->AddChild(mainN, CreatePropertyN(*property));
   }

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "material" node for GDML

XMLNodePointer_t TGDMLWrite::CreateMaterialN(TGeoMaterial * material, TString mname)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "material", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", mname);
   Double_t valZ = material->GetZ();
   //Z can't be zero in Geant4 so this is workaround for vacuum
   TString tmpname = mname;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   tmpname.ToLower();
   if (valZ < 1) {
      if (tmpname == "vacuum") {
         valZ = 1;
      } else {
         if (fgG4Compatibility == kTRUE) {
            Info("CreateMaterialN", "WARNING! value of Z in %s material can't be < 1 in Geant4, that is why it was changed to 1, please check it manually! ",
                 mname.Data());
            valZ = 1;
         } else {
            Info("CreateMaterialN", "WARNING! value of Z in %s material can't be < 1 in Geant4", mname.Data());
         }
      }
   }
   fGdmlE->NewAttr(mainN, nullptr, "Z", TString::Format(fltPrecision.Data(), valZ)); //material->GetZ()));
   fGdmlE->AddChild(mainN, CreateDN(material->GetDensity()));
   fGdmlE->AddChild(mainN, CreateAtomN(material->GetA()));
   // Create properties if any
   TList const &properties = material->GetProperties();
   if (properties.GetSize()) {
      TIter next(&properties);
      TNamed *property;
      while ((property = (TNamed*)next()))
        fGdmlE->AddChild(mainN, CreatePropertyN(*property));
   }
   // Write CONST properties
   TList const &const_properties = material->GetConstProperties();
   if (const_properties.GetSize()) {
      TIter next(&const_properties);
      TNamed *property;
      while ((property = (TNamed*)next()))
        fGdmlE->AddChild(mainN, CreatePropertyN(*property));
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "box" node for GDML

XMLNodePointer_t TGDMLWrite::CreateBoxN(TGeoBBox * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "box", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDX(), "DX", lname) ||
       IsNullParam(geoShape->GetDY(), "DY", lname) ||
       IsNullParam(geoShape->GetDZ(), "DZ", lname)) {
      return nullptr;
   }
   fGdmlE->NewAttr(mainN, nullptr, "x", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDX()));
   fGdmlE->NewAttr(mainN, nullptr, "y", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDY()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDZ()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "paraboloid" node for GDML

XMLNodePointer_t TGDMLWrite::CreateParaboloidN(TGeoParaboloid * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "paraboloid", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRhi(), "Rhi", lname) ||
       IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }
   fGdmlE->NewAttr(mainN, nullptr, "rlo", TString::Format(fltPrecision.Data(), geoShape->GetRlo()));
   fGdmlE->NewAttr(mainN, nullptr, "rhi", TString::Format(fltPrecision.Data(), geoShape->GetRhi()));
   fGdmlE->NewAttr(mainN, nullptr, "dz", TString::Format(fltPrecision.Data(), geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "sphere" node for GDML

XMLNodePointer_t TGDMLWrite::CreateSphereN(TGeoSphere * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "sphere", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRmax(), "Rmax", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi2() - geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "starttheta", TString::Format(fltPrecision.Data(), geoShape->GetTheta1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltatheta", TString::Format(fltPrecision.Data(), geoShape->GetTheta2() - geoShape->GetTheta1()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "arb8" node for GDML

XMLNodePointer_t TGDMLWrite::CreateArb8N(TGeoArb8 * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "arb8", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "v1x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[0]));
   fGdmlE->NewAttr(mainN, nullptr, "v1y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[1]));
   fGdmlE->NewAttr(mainN, nullptr, "v2x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[2]));
   fGdmlE->NewAttr(mainN, nullptr, "v2y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[3]));
   fGdmlE->NewAttr(mainN, nullptr, "v3x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[4]));
   fGdmlE->NewAttr(mainN, nullptr, "v3y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[5]));
   fGdmlE->NewAttr(mainN, nullptr, "v4x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[6]));
   fGdmlE->NewAttr(mainN, nullptr, "v4y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[7]));
   fGdmlE->NewAttr(mainN, nullptr, "v5x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[8]));
   fGdmlE->NewAttr(mainN, nullptr, "v5y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[9]));
   fGdmlE->NewAttr(mainN, nullptr, "v6x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[10]));
   fGdmlE->NewAttr(mainN, nullptr, "v6y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[11]));
   fGdmlE->NewAttr(mainN, nullptr, "v7x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[12]));
   fGdmlE->NewAttr(mainN, nullptr, "v7y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[13]));
   fGdmlE->NewAttr(mainN, nullptr, "v8x", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[14]));
   fGdmlE->NewAttr(mainN, nullptr, "v8y", TString::Format(fltPrecision.Data(), geoShape->GetVertices()[15]));
   fGdmlE->NewAttr(mainN, nullptr, "dz", TString::Format(fltPrecision.Data(), geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "cone" node for GDML from TGeoConeSeg object

XMLNodePointer_t TGDMLWrite::CreateConeN(TGeoConeSeg * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "cone", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "rmin1", TString::Format(fltPrecision.Data(), geoShape->GetRmin1()));
   fGdmlE->NewAttr(mainN, nullptr, "rmin2", TString::Format(fltPrecision.Data(), geoShape->GetRmin2()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax1", TString::Format(fltPrecision.Data(), geoShape->GetRmax1()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax2", TString::Format(fltPrecision.Data(), geoShape->GetRmax2()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi2() - geoShape->GetPhi1()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "cone" node for GDML from TGeoCone object

XMLNodePointer_t TGDMLWrite::CreateConeN(TGeoCone * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "cone", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "rmin1", TString::Format(fltPrecision.Data(), geoShape->GetRmin1()));
   fGdmlE->NewAttr(mainN, nullptr, "rmin2", TString::Format(fltPrecision.Data(), geoShape->GetRmin2()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax1", TString::Format(fltPrecision.Data(), geoShape->GetRmax1()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax2", TString::Format(fltPrecision.Data(), geoShape->GetRmax2()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format("%i", 0));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format("%i", 360));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "para" node for GDML

XMLNodePointer_t TGDMLWrite::CreateParaN(TGeoPara * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "para", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", GenName(geoShape->GetName(), TString::Format("%p", geoShape)));

   fGdmlE->NewAttr(mainN, nullptr, "x", TString::Format(fltPrecision.Data(), 2 * geoShape->GetX()));
   fGdmlE->NewAttr(mainN, nullptr, "y", TString::Format(fltPrecision.Data(), 2 * geoShape->GetY()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetZ()));
   fGdmlE->NewAttr(mainN, nullptr, "alpha", TString::Format(fltPrecision.Data(), geoShape->GetAlpha()));
   fGdmlE->NewAttr(mainN, nullptr, "theta", TString::Format(fltPrecision.Data(), geoShape->GetTheta()));
   fGdmlE->NewAttr(mainN, nullptr, "phi", TString::Format(fltPrecision.Data(), geoShape->GetPhi()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "trap" node for GDML

XMLNodePointer_t TGDMLWrite::CreateTrapN(TGeoTrap * geoShape)
{
   XMLNodePointer_t mainN;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);

   //if one base equals 0 create Arb8 instead of trap
   if ((geoShape->GetBl1() == 0 || geoShape->GetTl1() == 0 || geoShape->GetH1() == 0) ||
       (geoShape->GetBl2() == 0 || geoShape->GetTl2() == 0 || geoShape->GetH2() == 0)) {
      mainN = CreateArb8N(geoShape);
      return mainN;
   }

   //if is twisted then create arb8
   if (geoShape->IsTwisted()) {
      mainN = CreateArb8N((TGeoArb8 *) geoShape);
      return mainN;
   }

   mainN = fGdmlE->NewChild(nullptr, nullptr, "trap", nullptr);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "theta", TString::Format(fltPrecision.Data(), geoShape->GetTheta()));
   fGdmlE->NewAttr(mainN, nullptr, "phi", TString::Format(fltPrecision.Data(), geoShape->GetPhi()));
   fGdmlE->NewAttr(mainN, nullptr, "x1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetBl1()));
   fGdmlE->NewAttr(mainN, nullptr, "x2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetTl1()));
   fGdmlE->NewAttr(mainN, nullptr, "x3", TString::Format(fltPrecision.Data(), 2 * geoShape->GetBl2()));
   fGdmlE->NewAttr(mainN, nullptr, "x4", TString::Format(fltPrecision.Data(), 2 * geoShape->GetTl2()));
   fGdmlE->NewAttr(mainN, nullptr, "y1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetH1()));
   fGdmlE->NewAttr(mainN, nullptr, "y2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetH2()));

   fGdmlE->NewAttr(mainN, nullptr, "alpha1", TString::Format(fltPrecision.Data(), geoShape->GetAlpha1()));
   fGdmlE->NewAttr(mainN, nullptr, "alpha2", TString::Format(fltPrecision.Data(), geoShape->GetAlpha2()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "twistedtrap" node for GDML

XMLNodePointer_t TGDMLWrite::CreateTwistedTrapN(TGeoGtra * geoShape)
{
   XMLNodePointer_t mainN;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);

   //if one base equals 0 create Arb8 instead of twisted trap
   if ((geoShape->GetBl1() == 0 && geoShape->GetTl1() == 0 && geoShape->GetH1() == 0) ||
       (geoShape->GetBl2() == 0 && geoShape->GetTl2() == 0 && geoShape->GetH2() == 0)) {
      mainN = CreateArb8N((TGeoArb8 *) geoShape);
      return mainN;
   }

   //if is twisted then create arb8
   if (geoShape->IsTwisted()) {
      mainN = CreateArb8N((TGeoArb8 *) geoShape);
      return mainN;
   }

   //if parameter twistAngle (PhiTwist) equals zero create trap node
   if (geoShape->GetTwistAngle() == 0) {
      mainN = CreateTrapN((TGeoTrap *) geoShape);
      return mainN;
   }

   mainN = fGdmlE->NewChild(nullptr, nullptr, "twistedtrap", nullptr);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "Theta", TString::Format(fltPrecision.Data(), geoShape->GetTheta()));
   fGdmlE->NewAttr(mainN, nullptr, "Phi", TString::Format(fltPrecision.Data(), geoShape->GetPhi()));
   fGdmlE->NewAttr(mainN, nullptr, "x1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetBl1()));
   fGdmlE->NewAttr(mainN, nullptr, "x2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetTl1()));
   fGdmlE->NewAttr(mainN, nullptr, "x3", TString::Format(fltPrecision.Data(), 2 * geoShape->GetBl2()));
   fGdmlE->NewAttr(mainN, nullptr, "x4", TString::Format(fltPrecision.Data(), 2 * geoShape->GetTl2()));
   fGdmlE->NewAttr(mainN, nullptr, "y1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetH1()));
   fGdmlE->NewAttr(mainN, nullptr, "y2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetH2()));

   fGdmlE->NewAttr(mainN, nullptr, "Alph", TString::Format(fltPrecision.Data(), geoShape->GetAlpha1()));

   //check if alpha1 equals to alpha2 (converting to string - to avoid problems with floats)
   if (TString::Format(fltPrecision.Data(), geoShape->GetAlpha1()) != TString::Format(fltPrecision.Data(), geoShape->GetAlpha2())) {
      Info("CreateTwistedTrapN",
           "ERROR! Object %s is not exported correctly because parameter Alpha2 is not declared in GDML schema",
           lname.Data());
   }
   //fGdmlE->NewAttr(mainN,0, "alpha2", TString::Format(fltPrecision.Data(), geoShape->GetAlpha2()));
   fGdmlE->NewAttr(mainN, nullptr, "PhiTwist", TString::Format(fltPrecision.Data(), geoShape->GetTwistAngle()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "trd" node for GDML from object TGeoTrd1

XMLNodePointer_t TGDMLWrite::CreateTrdN(TGeoTrd1 * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "trd", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "x1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDx1()));
   fGdmlE->NewAttr(mainN, nullptr, "x2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDx2()));
   fGdmlE->NewAttr(mainN, nullptr, "y1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDy()));
   fGdmlE->NewAttr(mainN, nullptr, "y2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDy()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "trd" node for GDML from object TGeoTrd2

XMLNodePointer_t TGDMLWrite::CreateTrdN(TGeoTrd2 * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "trd", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "x1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDx1()));
   fGdmlE->NewAttr(mainN, nullptr, "x2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDx2()));
   fGdmlE->NewAttr(mainN, nullptr, "y1", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDy1()));
   fGdmlE->NewAttr(mainN, nullptr, "y2", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDy2()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "tube" node for GDML  from  object TGeoTubeSeg

XMLNodePointer_t TGDMLWrite::CreateTubeN(TGeoTubeSeg * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "tube", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRmax(), "Rmax", lname) ||
       IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(),  2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi2() - geoShape->GetPhi1()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "cutTube" node for GDML

XMLNodePointer_t TGDMLWrite::CreateCutTubeN(TGeoCtub * geoShape)
{
   XMLNodePointer_t mainN;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);

   mainN = fGdmlE->NewChild(nullptr, nullptr, "cutTube", nullptr);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRmax(), "Rmax", lname) ||
       IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }
   //This is not needed, because cutTube is already supported by Geant4 9.5
   if (fgG4Compatibility == kTRUE && kFALSE) {
      TGeoShape * fakeCtub = CreateFakeCtub(geoShape);
      mainN = ChooseObject(fakeCtub);

      //register name for cuttube shape (so it will be found during volume export)
      lname = fNameList->fLst[TString::Format("%p", fakeCtub)];
      fNameList->fLst[TString::Format("%p", geoShape)] = lname;
      Info("CreateCutTubeN", "WARNING! %s - CutTube was replaced by intersection of TGeoTubSeg and two TGeoBBoxes",
           lname.Data());
      return mainN;
   }
   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi2() - geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "lowX", TString::Format(fltPrecision.Data(), geoShape->GetNlow()[0]));
   fGdmlE->NewAttr(mainN, nullptr, "lowY", TString::Format(fltPrecision.Data(), geoShape->GetNlow()[1]));
   fGdmlE->NewAttr(mainN, nullptr, "lowZ", TString::Format(fltPrecision.Data(), geoShape->GetNlow()[2]));
   fGdmlE->NewAttr(mainN, nullptr, "highX", TString::Format(fltPrecision.Data(), geoShape->GetNhigh()[0]));
   fGdmlE->NewAttr(mainN, nullptr, "highY", TString::Format(fltPrecision.Data(), geoShape->GetNhigh()[1]));
   fGdmlE->NewAttr(mainN, nullptr, "highZ", TString::Format(fltPrecision.Data(), geoShape->GetNhigh()[2]));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "tube" node for GDML from  object TGeoTube

XMLNodePointer_t TGDMLWrite::CreateTubeN(TGeoTube * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "tube", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRmax(), "Rmax", lname) ||
       IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format("%i", 0));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format("%i", 360));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "zplane" node for GDML

XMLNodePointer_t TGDMLWrite::CreateZplaneN(Double_t z, Double_t rmin, Double_t rmax)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "zplane", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);

   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), z));
   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), rmin));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), rmax));

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "polycone" node for GDML

XMLNodePointer_t TGDMLWrite::CreatePolyconeN(TGeoPcon * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "polycone", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);

   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetDphi()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   Int_t nZPlns = geoShape->GetNz();
   for (Int_t it = 0; it < nZPlns; it++) {
      //add zplane child node
      fGdmlE->AddChild(mainN, CreateZplaneN(geoShape->GetZ(it), geoShape->GetRmin(it), geoShape->GetRmax(it)));
      //compare actual plane and next plane
      if ((it < nZPlns - 1) && (geoShape->GetZ(it) == geoShape->GetZ(it + 1))) {
         //rmin of actual is greater then rmax of next one
         //        |   |rmax next
         //  __ ...|   |...  __   < rmin actual
         // |  |            |  |
         if (geoShape->GetRmin(it) > geoShape->GetRmax(it + 1)) {
            //adding plane from rmax next to rmin actual at the same z position
            if (fgG4Compatibility == kTRUE) {
               fGdmlE->AddChild(mainN, CreateZplaneN(geoShape->GetZ(it), geoShape->GetRmax(it + 1), geoShape->GetRmin(it)));
               Info("CreatePolyconeN", "WARNING! One plane was added to %s solid to be compatible with Geant4", lname.Data());
            } else {
               Info("CreatePolyconeN", "WARNING! Solid %s definition seems not contiguous may cause problems in Geant4", lname.Data());
            }

         }
         //rmin of next is greater then rmax of actual
         //  |  |         |  |
         //  |  |...___...|  |  rmin next
         //        |   |     > rmax act
         if (geoShape->GetRmin(it + 1) > geoShape->GetRmax(it)) {
            //adding plane from rmax act to rmin next at the same z position
            if (fgG4Compatibility == kTRUE) {
               fGdmlE->AddChild(mainN, CreateZplaneN(geoShape->GetZ(it), geoShape->GetRmax(it), geoShape->GetRmin(it + 1)));
               Info("CreatePolyconeN", "WARNING! One plane was added to %s solid to be compatible with Geant4", lname.Data());
            } else {
               Info("CreatePolyconeN", "WARNING! Solid %s definition seems not contiguous may cause problems in Geant4", lname.Data());
            }
         }
      }
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "torus" node for GDML

XMLNodePointer_t TGDMLWrite::CreateTorusN(TGeoTorus * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "torus", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetRmax(), "Rmax", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "rtor", TString::Format(fltPrecision.Data(), geoShape->GetR()));
   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetDphi()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "polyhedra" node for GDML

XMLNodePointer_t TGDMLWrite::CreatePolyhedraN(TGeoPgon * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "polyhedra", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", GenName(geoShape->GetName(), TString::Format("%p", geoShape)));

   fGdmlE->NewAttr(mainN, nullptr, "startphi", TString::Format(fltPrecision.Data(), geoShape->GetPhi1()));
   fGdmlE->NewAttr(mainN, nullptr, "deltaphi", TString::Format(fltPrecision.Data(), geoShape->GetDphi()));
   fGdmlE->NewAttr(mainN, nullptr, "numsides", TString::Format("%i", geoShape->GetNedges()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   for (Int_t it = 0; it < geoShape->GetNz(); it++) {
      //add zplane child node
      fGdmlE->AddChild(mainN, CreateZplaneN(geoShape->GetZ(it), geoShape->GetRmin(it), geoShape->GetRmax(it)));
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "eltube" node for GDML

XMLNodePointer_t TGDMLWrite::CreateEltubeN(TGeoEltu * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "eltube", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetA(), "A", lname) ||
       IsNullParam(geoShape->GetB(), "B", lname) ||
       IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }

   fGdmlE->NewAttr(mainN, nullptr, "dx", TString::Format(fltPrecision.Data(), geoShape->GetA()));
   fGdmlE->NewAttr(mainN, nullptr, "dy", TString::Format(fltPrecision.Data(), geoShape->GetB()));
   fGdmlE->NewAttr(mainN, nullptr, "dz", TString::Format(fltPrecision.Data(), geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "hype" node for GDML

XMLNodePointer_t TGDMLWrite::CreateHypeN(TGeoHype * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "hype", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);
   if (IsNullParam(geoShape->GetDz(), "Dz", lname)) {
      return nullptr;
   }


   fGdmlE->NewAttr(mainN, nullptr, "rmin", TString::Format(fltPrecision.Data(), geoShape->GetRmin()));
   fGdmlE->NewAttr(mainN, nullptr, "rmax", TString::Format(fltPrecision.Data(), geoShape->GetRmax()));
   fGdmlE->NewAttr(mainN, nullptr, "inst", TString::Format(fltPrecision.Data(), geoShape->GetStIn()));
   fGdmlE->NewAttr(mainN, nullptr, "outst", TString::Format(fltPrecision.Data(), geoShape->GetStOut()));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), 2 * geoShape->GetDz()));

   fGdmlE->NewAttr(mainN, nullptr, "aunit", "deg");
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "xtru" node for GDML

XMLNodePointer_t TGDMLWrite::CreateXtrusionN(TGeoXtru * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "xtru", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TString lname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   fGdmlE->NewAttr(mainN, nullptr, "name", lname);

   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);
   XMLNodePointer_t childN;
   Int_t vertNum =  geoShape->GetNvert();
   Int_t secNum = geoShape->GetNz();
   if (vertNum < 3 || secNum < 2) {
      Info("CreateXtrusionN", "ERROR! TGeoXtru %s has only %i vertices and %i sections. It was not exported",
           lname.Data(), vertNum, secNum);
      mainN = nullptr;
      return mainN;
   }
   for (Int_t it = 0; it < vertNum; it++) {
      //add twoDimVertex child node
      childN = fGdmlE->NewChild(nullptr, nullptr, "twoDimVertex", nullptr);
      fGdmlE->NewAttr(childN, nullptr, "x", TString::Format(fltPrecision.Data(), geoShape->GetX(it)));
      fGdmlE->NewAttr(childN, nullptr, "y", TString::Format(fltPrecision.Data(), geoShape->GetY(it)));
      fGdmlE->AddChild(mainN, childN);
   }
   for (Int_t it = 0; it < secNum; it++) {
      //add section child node
      childN = fGdmlE->NewChild(nullptr, nullptr, "section", nullptr);
      fGdmlE->NewAttr(childN, nullptr, "zOrder", TString::Format("%i", it));
      fGdmlE->NewAttr(childN, nullptr, "zPosition", TString::Format(fltPrecision.Data(), geoShape->GetZ(it)));
      fGdmlE->NewAttr(childN, nullptr, "xOffset", TString::Format(fltPrecision.Data(), geoShape->GetXOffset(it)));
      fGdmlE->NewAttr(childN, nullptr, "yOffset", TString::Format(fltPrecision.Data(), geoShape->GetYOffset(it)));
      fGdmlE->NewAttr(childN, nullptr, "scalingFactor", TString::Format(fltPrecision.Data(), geoShape->GetScale(it)));
      fGdmlE->AddChild(mainN, childN);
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "ellipsoid" node for GDML
/// this is a special case, because ellipsoid is not defined in ROOT
/// so when intersection of scaled sphere and TGeoBBox is found,
/// it is considered as an ellipsoid

XMLNodePointer_t TGDMLWrite::CreateEllipsoidN(TGeoCompositeShape * geoShape, TString elName)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "ellipsoid", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   TGeoScaledShape *leftS = (TGeoScaledShape *)geoShape->GetBoolNode()->GetLeftShape(); //ScaledShape
   TGeoBBox *rightS = (TGeoBBox *)geoShape->GetBoolNode()->GetRightShape(); //BBox


   fGdmlE->NewAttr(mainN, nullptr, "name", elName.Data());
   Double_t sx = leftS->GetScale()->GetScale()[0];
   Double_t sy = leftS->GetScale()->GetScale()[1];
   Double_t radius = ((TGeoSphere *) leftS->GetShape())->GetRmax();

   Double_t ax, by, cz;
   cz = radius;
   ax = sx * radius;
   by = sy * radius;

   Double_t dz = rightS->GetDZ();
   Double_t zorig = rightS->GetOrigin()[2];
   Double_t zcut2 = dz + zorig;
   Double_t zcut1 = 2 * zorig - zcut2;


   fGdmlE->NewAttr(mainN, nullptr, "ax", TString::Format(fltPrecision.Data(), ax));
   fGdmlE->NewAttr(mainN, nullptr, "by", TString::Format(fltPrecision.Data(), by));
   fGdmlE->NewAttr(mainN, nullptr, "cz", TString::Format(fltPrecision.Data(), cz));
   fGdmlE->NewAttr(mainN, nullptr, "zcut1", TString::Format(fltPrecision.Data(), zcut1));
   fGdmlE->NewAttr(mainN, nullptr, "zcut2", TString::Format(fltPrecision.Data(), zcut2));
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "elcone" (elliptical cone) node for GDML
/// this is a special case, because elliptical cone is not defined in ROOT
/// so when scaled cone is found, it is considered as a elliptical cone

XMLNodePointer_t TGDMLWrite::CreateElConeN(TGeoScaledShape * geoShape)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "elcone", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", GenName(geoShape->GetName(), TString::Format("%p", geoShape)));
   Double_t zcut = ((TGeoCone *) geoShape->GetShape())->GetDz();
   Double_t rx1 = ((TGeoCone *) geoShape->GetShape())->GetRmax1();
   Double_t rx2 = ((TGeoCone *) geoShape->GetShape())->GetRmax2();
   Double_t zmax = zcut * ((rx1 + rx2) / (rx1 - rx2));
   Double_t z = zcut + zmax;

   Double_t sy = geoShape->GetScale()->GetScale()[1];
   Double_t ry1 = sy * rx1;

   std::string format(TString::Format("%s/%s", fltPrecision.Data(), fltPrecision.Data()).Data());
   fGdmlE->NewAttr(mainN, nullptr, "dx", TString::Format(format.c_str(), rx1, z));
   fGdmlE->NewAttr(mainN, nullptr, "dy", TString::Format(format.c_str(), ry1, z));
   fGdmlE->NewAttr(mainN, nullptr, "zmax", TString::Format(fltPrecision.Data(), zmax));
   fGdmlE->NewAttr(mainN, nullptr, "zcut", TString::Format(fltPrecision.Data(), zcut));
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "tessellated" (tessellated shape) node for GDML

XMLNodePointer_t TGDMLWrite::CreateTessellatedN(TGeoTessellated * geoShape)
{
   // add all vertices to the define section
   TString genname = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   for (int i = 0; i < geoShape->GetNvertices(); ++i) {
      auto vertex = geoShape->GetVertex(i);
      TString posName = TString::Format("%s_%d", genname.Data(), i);
      Xyz nodPos;
      nodPos.x = vertex[0];
      nodPos.y = vertex[1];
      nodPos.z = vertex[2];
      auto childN = CreatePositionN(posName.Data(), nodPos);
      fGdmlE->AddChild(fDefineNode, childN); //adding node to <define> node
   }
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "tessellated", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", genname.Data());
   fGdmlE->NewAttr(mainN, nullptr, "lunit", fDefault_lunit);

   XMLNodePointer_t childN;
   for (Int_t it = 0; it < geoShape->GetNfacets(); it++) {
      //add section child node
      auto facet = geoShape->GetFacet(it);
      bool triangular = facet.GetNvert() == 3;
      TString ntype = (triangular) ? "triangular" : "quadrangular";
      childN = fGdmlE->NewChild(nullptr, nullptr, ntype.Data(), nullptr);
      fGdmlE->NewAttr(childN, nullptr, "vertex1", TString::Format("%s_%d", genname.Data(), facet.GetVertexIndex(0)));
      fGdmlE->NewAttr(childN, nullptr, "vertex2", TString::Format("%s_%d", genname.Data(), facet.GetVertexIndex(1)));
      fGdmlE->NewAttr(childN, nullptr, "vertex3", TString::Format("%s_%d", genname.Data(), facet.GetVertexIndex(2)));
      if (!triangular)
         fGdmlE->NewAttr(childN, nullptr, "vertex4", TString::Format("%s_%d", genname.Data(), facet.GetVertexIndex(3)));
      fGdmlE->NewAttr(childN, nullptr, "type", "ABSOLUTE");
      fGdmlE->AddChild(mainN, childN);
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates common part of union intersection and subtraction nodes

XMLNodePointer_t TGDMLWrite::CreateCommonBoolN(TGeoCompositeShape *geoShape)
{
   XMLNodePointer_t mainN, ndR, ndL, childN;
   TString nodeName = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   TString lboolType;
   TGeoBoolNode::EGeoBoolType boolType = geoShape->GetBoolNode()->GetBooleanOperator();
   switch (boolType) {
      case TGeoBoolNode::kGeoUnion:
         lboolType = "union";
         break;
      case TGeoBoolNode::kGeoSubtraction:
         lboolType = "subtraction";
         break;
      case TGeoBoolNode::kGeoIntersection:
         lboolType = "intersection";
         break;
   }

   TGDMLWrite::Xyz lrot = GetXYZangles(geoShape->GetBoolNode()->GetLeftMatrix()->Inverse().GetRotationMatrix());
   const Double_t  *ltr = geoShape->GetBoolNode()->GetLeftMatrix()->GetTranslation();
   TGDMLWrite::Xyz rrot = GetXYZangles(geoShape->GetBoolNode()->GetRightMatrix()->Inverse().GetRotationMatrix());
   const Double_t  *rtr = geoShape->GetBoolNode()->GetRightMatrix()->GetTranslation();

   //specific case!
   //Ellipsoid tag preparing
   //if left == TGeoScaledShape AND right  == TGeoBBox
   //   AND if TGeoScaledShape->GetShape == TGeoSphere
   TGeoShape *leftS = geoShape->GetBoolNode()->GetLeftShape();
   TGeoShape *rightS = geoShape->GetBoolNode()->GetRightShape();
   if (strcmp(leftS->ClassName(), "TGeoScaledShape") == 0 &&
       strcmp(rightS->ClassName(), "TGeoBBox") == 0) {
      if (strcmp(((TGeoScaledShape *)leftS)->GetShape()->ClassName(), "TGeoSphere") == 0) {
         if (lboolType == "intersection") {
            mainN = CreateEllipsoidN(geoShape, nodeName);
            return mainN;
         }
      }
   }

   Xyz translL, translR;
   //translation
   translL.x = ltr[0];
   translL.y = ltr[1];
   translL.z = ltr[2];
   translR.x = rtr[0];
   translR.y = rtr[1];
   translR.z = rtr[2];

   //left and right nodes are created here also their names are created
   ndL = ChooseObject(geoShape->GetBoolNode()->GetLeftShape());
   ndR = ChooseObject(geoShape->GetBoolNode()->GetRightShape());

   //retrieve left and right node names by their pointer to make reference
   TString lname = fNameList->fLst[TString::Format("%p", geoShape->GetBoolNode()->GetLeftShape())];
   TString rname = fNameList->fLst[TString::Format("%p", geoShape->GetBoolNode()->GetRightShape())];

   //left and right nodes appended to main structure of nodes (if they are not already there)
   if (ndL != nullptr) {
      fGdmlE->AddChild(fSolidsNode, ndL);
      fSolCnt++;
   } else {
      if (lname.Contains("missing_") || lname == "") {
         Info("CreateCommonBoolN", "ERROR! Left node is NULL - Boolean Shape will be skipped");
         return nullptr;
      }
   }
   if (ndR != nullptr) {
      fGdmlE->AddChild(fSolidsNode, ndR);
      fSolCnt++;
   } else {
      if (rname.Contains("missing_") || rname == "") {
         Info("CreateCommonBoolN", "ERROR! Right node is NULL - Boolean Shape will be skipped");
         return nullptr;
      }
   }

   //create union node and its child nodes (or intersection or subtraction)
   /* <union name="...">
    *   <first ref="left name" />
    *   <second ref="right name" />
    *   <firstposition .../>
    *   <firstrotation .../>
    *   <position .../>
    *   <rotation .../>
    * </union>
   */
   mainN = fGdmlE->NewChild(nullptr, nullptr, lboolType.Data(), nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", nodeName);

   //<first> (left)
   childN = fGdmlE->NewChild(nullptr, nullptr, "first", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", lname);
   fGdmlE->AddChild(mainN, childN);

   //<second> (right)
   childN = fGdmlE->NewChild(nullptr, nullptr, "second", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", rname);
   fGdmlE->AddChild(mainN, childN);

   //<firstposition> (left)
   if ((translL.x != 0.0) || (translL.y != 0.0) || (translL.z != 0.0)) {
      childN = CreatePositionN((nodeName + lname + "pos").Data(), translL, "firstposition");
      fGdmlE->AddChild(mainN, childN);
   }
   //<firstrotation> (left)
   if ((lrot.x != 0.0) || (lrot.y != 0.0) || (lrot.z != 0.0)) {
      childN = CreateRotationN((nodeName + lname + "rot").Data(), lrot, "firstrotation");
      fGdmlE->AddChild(mainN, childN);
   }
   //<position> (right)
   if ((translR.x != 0.0) || (translR.y != 0.0) || (translR.z != 0.0)) {
      childN = CreatePositionN((nodeName + rname + "pos").Data(), translR, "position");
      fGdmlE->AddChild(mainN, childN);
   }
   //<rotation> (right)
   if ((rrot.x != 0.0) || (rrot.y != 0.0) || (rrot.z != 0.0)) {
      childN = CreateRotationN((nodeName + rname + "rot").Data(), rrot, "rotation");
      fGdmlE->AddChild(mainN, childN);
   }

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "opticalsurface" node for GDML

XMLNodePointer_t TGDMLWrite::CreateOpticalSurfaceN(TGeoOpticalSurface * geoSurf)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "opticalsurface", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", geoSurf->GetName());
   fGdmlE->NewAttr(mainN, nullptr, "model", TGeoOpticalSurface::ModelToString(geoSurf->GetModel()));
   fGdmlE->NewAttr(mainN, nullptr, "finish", TGeoOpticalSurface::FinishToString(geoSurf->GetFinish()));
   fGdmlE->NewAttr(mainN, nullptr, "type", TGeoOpticalSurface::TypeToString(geoSurf->GetType()));
   fGdmlE->NewAttr(mainN, nullptr, "value", TString::Format(fltPrecision.Data(), geoSurf->GetValue()));

   // Write properties
   TList const &properties = geoSurf->GetProperties();
   if (properties.GetSize()) {
      TIter next(&properties);
      TNamed *property;
      while ((property = (TNamed*)next()))
        fGdmlE->AddChild(mainN, CreatePropertyN(*property));
   }
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "skinsurface" node for GDML

XMLNodePointer_t TGDMLWrite::CreateSkinSurfaceN(TGeoSkinSurface * geoSurf)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "skinsurface", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", geoSurf->GetName());
   fGdmlE->NewAttr(mainN, nullptr, "surfaceproperty", geoSurf->GetTitle());
   // Cretate the logical volume reference node
   auto childN = fGdmlE->NewChild(nullptr, nullptr, "volumeref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", geoSurf->GetVolume()->GetName());
   fGdmlE->AddChild(mainN, childN);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "bordersurface" node for GDML

XMLNodePointer_t TGDMLWrite::CreateBorderSurfaceN(TGeoBorderSurface * geoSurf)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "bordersurface", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", geoSurf->GetName());
   fGdmlE->NewAttr(mainN, nullptr, "surfaceproperty", geoSurf->GetTitle());
   // Cretate the logical volume reference node
   auto childN = fGdmlE->NewChild(nullptr, nullptr, "physvolref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", geoSurf->GetNode1()->GetName());
   fGdmlE->NewAttr(childN, nullptr, "ref", geoSurf->GetNode2()->GetName());
   fGdmlE->AddChild(mainN, childN);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "position" kind of node for GDML

XMLNodePointer_t TGDMLWrite::CreatePositionN(const char * name, Xyz position, const char * type, const char * unit)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, type, nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   fGdmlE->NewAttr(mainN, nullptr, "x", TString::Format(fltPrecision.Data(), position.x));
   fGdmlE->NewAttr(mainN, nullptr, "y", TString::Format(fltPrecision.Data(), position.y));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), position.z));
   fGdmlE->NewAttr(mainN, nullptr, "unit", unit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "rotation" kind of node for GDML

XMLNodePointer_t TGDMLWrite::CreateRotationN(const char * name, Xyz rotation, const char * type, const char * unit)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, type, nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   fGdmlE->NewAttr(mainN, nullptr, "x", TString::Format(fltPrecision.Data(), rotation.x));
   fGdmlE->NewAttr(mainN, nullptr, "y", TString::Format(fltPrecision.Data(), rotation.y));
   fGdmlE->NewAttr(mainN, nullptr, "z", TString::Format(fltPrecision.Data(), rotation.z));
   fGdmlE->NewAttr(mainN, nullptr, "unit", unit);
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "matrix" kind of node for GDML

XMLNodePointer_t TGDMLWrite::CreateMatrixN(TGDMLMatrix const *matrix)
{
   std::stringstream vals;
   size_t cols = matrix->GetCols();
   size_t rows = matrix->GetRows();
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "matrix", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", matrix->GetName());
   fGdmlE->NewAttr(mainN, nullptr, "coldim", TString::Format("%zu", cols));
   for(size_t i=0; i<rows; ++i)  {
     for(size_t j=0; j<cols; ++j)  {
       vals << matrix->Get(i,j);
       if ( j < cols-1 ) vals << ' ';
     }
     if ( i < rows-1 ) vals << '\n';
   }
   fGdmlE->NewAttr(mainN, nullptr, "values", vals.str().c_str());
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "constant" kind of node for GDML

XMLNodePointer_t TGDMLWrite::CreateConstantN(const char *name, Double_t value)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "constant", nullptr);
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   fGdmlE->NewAttr(mainN, nullptr, "value", TString::Format(fltPrecision.Data(), value));
   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "setup" node for GDML

XMLNodePointer_t TGDMLWrite::CreateSetupN(const char * topVolName, const char * name, const char * version)
{
   XMLNodePointer_t setupN = fGdmlE->NewChild(nullptr, nullptr, "setup", nullptr);
   fGdmlE->NewAttr(setupN, nullptr, "name", name);
   fGdmlE->NewAttr(setupN, nullptr, "version", version);
   XMLNodePointer_t fworldN = fGdmlE->NewChild(setupN, nullptr, "world", nullptr);
   fGdmlE->NewAttr(fworldN, nullptr, "ref", topVolName);
   return setupN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "volume" node for GDML

XMLNodePointer_t TGDMLWrite::StartVolumeN(const char * name, const char * solid, const char * material)
{
   XMLNodePointer_t childN;
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "volume", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);

   childN = fGdmlE->NewChild(nullptr, nullptr, "materialref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", material);
   fGdmlE->AddChild(mainN, childN);

   childN = fGdmlE->NewChild(nullptr, nullptr, "solidref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", solid);
   fGdmlE->AddChild(mainN, childN);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "assembly" node for GDML

XMLNodePointer_t TGDMLWrite::StartAssemblyN(const char * name)
{
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "assembly", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "physvol" node for GDML

XMLNodePointer_t TGDMLWrite::CreatePhysVolN(const char *name, Int_t copyno, const char * volref, const char * posref, const char * rotref, XMLNodePointer_t scaleN)
{
   fPhysVolCnt++;
   XMLNodePointer_t childN;
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "physvol", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "name", name);
   fGdmlE->NewAttr(mainN, nullptr, "copynumber", TString::Format("%d",copyno));

   childN = fGdmlE->NewChild(nullptr, nullptr, "volumeref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", volref);
   fGdmlE->AddChild(mainN, childN);

   childN = fGdmlE->NewChild(nullptr, nullptr, "positionref", nullptr);
   fGdmlE->NewAttr(childN, nullptr, "ref", posref);
   fGdmlE->AddChild(mainN, childN);

   //if is not empty string add this node
   if (strcmp(rotref, "") != 0) {
      childN = fGdmlE->NewChild(nullptr, nullptr, "rotationref", nullptr);
      fGdmlE->NewAttr(childN, nullptr, "ref", rotref);
      fGdmlE->AddChild(mainN, childN);
   }
   if (scaleN != nullptr) {
      fGdmlE->AddChild(mainN, scaleN);
   }

   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates "divisionvol" node for GDML

XMLNodePointer_t TGDMLWrite::CreateDivisionN(Double_t offset, Double_t width, Int_t number, const char * axis, const char * unit, const char * volref)
{
   XMLNodePointer_t childN = 0;
   XMLNodePointer_t mainN = fGdmlE->NewChild(nullptr, nullptr, "divisionvol", nullptr);
   fGdmlE->NewAttr(mainN, nullptr, "axis", axis);
   fGdmlE->NewAttr(mainN, nullptr, "number", TString::Format("%i", number));
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);
   if (fgG4Compatibility  == kTRUE) {
      //if eg. full length is 20 and width * number = 20,0001 problem in geant4
      //unit is either in cm or degrees nothing else
      width = (floor(width * 1E4)) * 1E-4;
      if ((offset >= 0.) && (strcmp(axis, "kPhi") == 0)) {
         Int_t offsetI = (Int_t) offset;
         Double_t decimals = offset - offsetI;
         //put to range from 0 to 360 add decimals and then put to range 0 -> -360
         offset = (offsetI % 360) + decimals - 360;
      }
   }
   fGdmlE->NewAttr(mainN, nullptr, "width", TString::Format(fltPrecision.Data(), width));

   fGdmlE->NewAttr(mainN, nullptr, "offset", TString::Format(fltPrecision.Data(), offset));
   fGdmlE->NewAttr(mainN, nullptr, "unit", unit);
   if (strcmp(volref, "") != 0) {
      childN = fGdmlE->NewChild(nullptr, nullptr, "volumeref", nullptr);
      fGdmlE->NewAttr(childN, nullptr, "ref", volref);
   }
   fGdmlE->AddChild(mainN, childN);


   return mainN;
}

////////////////////////////////////////////////////////////////////////////////
/// Chooses the object and method that should be used for processing object

XMLNodePointer_t TGDMLWrite::ChooseObject(TGeoShape *geoShape)
{
   const char * clsname = geoShape->ClassName();
   XMLNodePointer_t solidN;

   if (CanProcess((TObject *)geoShape) == kFALSE) {
      return nullptr;
   }

   //process different shapes
   if (strcmp(clsname, "TGeoBBox") == 0) {
      solidN = CreateBoxN((TGeoBBox*) geoShape);
   } else if (strcmp(clsname, "TGeoParaboloid") == 0) {
      solidN = CreateParaboloidN((TGeoParaboloid*) geoShape);
   } else if (strcmp(clsname, "TGeoSphere") == 0) {
      solidN = CreateSphereN((TGeoSphere*) geoShape);
   } else if (strcmp(clsname, "TGeoArb8") == 0) {
      solidN = CreateArb8N((TGeoArb8*) geoShape);
   } else if (strcmp(clsname, "TGeoConeSeg") == 0) {
      solidN = CreateConeN((TGeoConeSeg*) geoShape);
   } else if (strcmp(clsname, "TGeoCone") == 0) {
      solidN = CreateConeN((TGeoCone*) geoShape);
   } else if (strcmp(clsname, "TGeoPara") == 0) {
      solidN = CreateParaN((TGeoPara*) geoShape);
   } else if (strcmp(clsname, "TGeoTrap") == 0) {
      solidN = CreateTrapN((TGeoTrap*) geoShape);
   } else if (strcmp(clsname, "TGeoGtra") == 0) {
      solidN = CreateTwistedTrapN((TGeoGtra*) geoShape);
   } else if (strcmp(clsname, "TGeoTrd1") == 0) {
      solidN = CreateTrdN((TGeoTrd1*) geoShape);
   } else if (strcmp(clsname, "TGeoTrd2") == 0) {
      solidN = CreateTrdN((TGeoTrd2*) geoShape);
   } else if (strcmp(clsname, "TGeoTubeSeg") == 0) {
      solidN = CreateTubeN((TGeoTubeSeg*) geoShape);
   } else if (strcmp(clsname, "TGeoCtub") == 0) {
      solidN = CreateCutTubeN((TGeoCtub*) geoShape);
   } else if (strcmp(clsname, "TGeoTube") == 0) {
      solidN = CreateTubeN((TGeoTube*) geoShape);
   } else if (strcmp(clsname, "TGeoPcon") == 0) {
      solidN = CreatePolyconeN((TGeoPcon*) geoShape);
   } else if (strcmp(clsname, "TGeoTorus") == 0) {
      solidN = CreateTorusN((TGeoTorus*) geoShape);
   } else if (strcmp(clsname, "TGeoPgon") == 0) {
      solidN = CreatePolyhedraN((TGeoPgon*) geoShape);
   } else if (strcmp(clsname, "TGeoEltu") == 0) {
      solidN = CreateEltubeN((TGeoEltu*) geoShape);
   } else if (strcmp(clsname, "TGeoHype") == 0) {
      solidN = CreateHypeN((TGeoHype*) geoShape);
   } else if (strcmp(clsname, "TGeoXtru") == 0) {
      solidN = CreateXtrusionN((TGeoXtru*) geoShape);
   } else if (strcmp(clsname, "TGeoTessellated") == 0) {
      solidN = CreateTessellatedN((TGeoTessellated*) geoShape);
   } else if (strcmp(clsname, "TGeoScaledShape") == 0) {
      TGeoScaledShape * geoscale = (TGeoScaledShape *) geoShape;
      TString scaleObjClsName = geoscale->GetShape()->ClassName();
      if (scaleObjClsName == "TGeoCone") {
         solidN = CreateElConeN((TGeoScaledShape*) geoShape);
      } else {
         Info("ChooseObject",
              "ERROR! TGeoScaledShape object is not possible to process correctly. %s object is processed without scale",
              scaleObjClsName.Data());
         solidN = ChooseObject(geoscale->GetShape());
         //Name has to be propagated to geoscale level pointer
         fNameList->fLst[TString::Format("%p", geoscale)] =
            fNameList->fLst[TString::Format("%p", geoscale->GetShape())];
      }
   } else if (strcmp(clsname, "TGeoCompositeShape") == 0) {
      solidN = CreateCommonBoolN((TGeoCompositeShape*) geoShape);
   } else if (strcmp(clsname, "TGeoUnion") == 0) {
      solidN = CreateCommonBoolN((TGeoCompositeShape*) geoShape);
   } else if (strcmp(clsname, "TGeoIntersection") == 0) {
      solidN = CreateCommonBoolN((TGeoCompositeShape*) geoShape);
   } else if (strcmp(clsname, "TGeoSubtraction") == 0) {
      solidN = CreateCommonBoolN((TGeoCompositeShape*) geoShape);
   } else {
      Info("ChooseObject", "ERROR! %s Solid CANNOT be processed, solid is NOT supported",
           clsname);
      solidN = nullptr;
   }
   if (solidN == nullptr) {
      if (fNameList->fLst[TString::Format("%p", geoShape)] == "") {
         TString missingName = geoShape->GetName();
         GenName("missing_" + missingName, TString::Format("%p", geoShape));
      } else {
         fNameList->fLst[TString::Format("%p", geoShape)] = "missing_" + fNameList->fLst[TString::Format("%p", geoShape)];
      }
   }

   return solidN;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves X Y Z angles from rotation matrix

TGDMLWrite::Xyz TGDMLWrite::GetXYZangles(const Double_t * rotationMatrix)
{
   TGDMLWrite::Xyz lxyz;
   Double_t a, b, c;
   Double_t rad = 180.0 / TMath::ACos(-1.0);
   const Double_t *r = rotationMatrix;
   Double_t cosb = TMath::Sqrt(r[0] * r[0] + r[1] * r[1]);
   if (cosb > 0.00001) {
      a = TMath::ATan2(r[5], r[8]) * rad;
      b = TMath::ATan2(-r[2], cosb) * rad;
      c = TMath::ATan2(r[1], r[0]) * rad;
   } else {
      a = TMath::ATan2(-r[7], r[4]) * rad;
      b = TMath::ATan2(-r[2], cosb) * rad;
      c = 0;
   }
   lxyz.x = a;
   lxyz.y = b;
   lxyz.z = c;
   return lxyz;
}

////////////////////////////////////////////////////////////////////////////////
/// Method creating cutTube as an intersection of tube and two boxes
/// - not used anymore because cutube is supported in Geant4 9.5

TGeoCompositeShape* TGDMLWrite::CreateFakeCtub(TGeoCtub* geoShape)
{
   Double_t rmin = geoShape->GetRmin();
   Double_t rmax = geoShape->GetRmax();
   Double_t z = geoShape->GetDz();
   Double_t startphi = geoShape->GetPhi1();
   Double_t deltaphi = geoShape->GetPhi2();
   Double_t x1 = geoShape->GetNlow()[0];
   Double_t y1 = geoShape->GetNlow()[1];
   Double_t z1 = geoShape->GetNlow()[2];
   Double_t x2 = geoShape->GetNhigh()[0];
   Double_t y2 = geoShape->GetNhigh()[1];
   Double_t z2 = geoShape->GetNhigh()[2];
   TString xname = geoShape->GetName();


   Double_t h0 = 2.*((TGeoBBox*)geoShape)->GetDZ();
   Double_t h1 = 2 * z;
   Double_t h2 = 2 * z;
   Double_t boxdx = 1E8 * (2 * rmax) + (2 * z);

   TGeoTubeSeg *T = new TGeoTubeSeg((xname + "T").Data(), rmin, rmax, h0, startphi, deltaphi);
   TGeoBBox *B1 = new TGeoBBox((xname + "B1").Data(), boxdx, boxdx, h1);
   TGeoBBox *B2 = new TGeoBBox((xname + "B2").Data(), boxdx, boxdx, h2);


   //first box position parameters
   Double_t phi1 = 360 - TMath::ATan2(x1, y1) * TMath::RadToDeg();
   Double_t theta1 = 360 - TMath::ATan2(sqrt(x1 * x1 + y1 * y1), z1) * TMath::RadToDeg();

   Double_t phi11 = TMath::ATan2(y1, x1) * TMath::RadToDeg() ;
   Double_t theta11 = TMath::ATan2(z1, sqrt(x1 * x1 + y1 * y1)) * TMath::RadToDeg() ;

   Double_t xpos1 = h1 * TMath::Cos((theta11) * TMath::DegToRad()) * TMath::Cos((phi11) * TMath::DegToRad()) * (-1);
   Double_t ypos1 = h1 * TMath::Cos((theta11) * TMath::DegToRad()) * TMath::Sin((phi11) * TMath::DegToRad()) * (-1);
   Double_t zpos1 = h1 * TMath::Sin((theta11) * TMath::DegToRad()) * (-1);

   //second box position parameters
   Double_t phi2 = 360 - TMath::ATan2(x2, y2) * TMath::RadToDeg();
   Double_t theta2 = 360 - TMath::ATan2(sqrt(x2 * x2 + y2 * y2), z2) * TMath::RadToDeg();

   Double_t phi21 = TMath::ATan2(y2, x2) * TMath::RadToDeg() ;
   Double_t theta21 = TMath::ATan2(z2, sqrt(x2 * x2 + y2 * y2)) * TMath::RadToDeg() ;

   Double_t xpos2 = h2 * TMath::Cos((theta21) * TMath::DegToRad()) * TMath::Cos((phi21) * TMath::DegToRad()) * (-1);
   Double_t ypos2 = h2 * TMath::Cos((theta21) * TMath::DegToRad()) * TMath::Sin((phi21) * TMath::DegToRad()) * (-1);
   Double_t zpos2 = h2 * TMath::Sin((theta21) * TMath::DegToRad()) * (-1);


   //positioning
   TGeoTranslation *t0 = new TGeoTranslation(0, 0, 0);
   TGeoTranslation *t1 = new TGeoTranslation(0 + xpos1, 0 + ypos1, 0 + (zpos1 - z));
   TGeoTranslation *t2 = new TGeoTranslation(0 + xpos2, 0 + ypos2, 0 + (zpos2 + z));
   TGeoRotation *r0 = new TGeoRotation((xname + "r0").Data());
   TGeoRotation *r1 = new TGeoRotation((xname + "r1").Data());
   TGeoRotation *r2 = new TGeoRotation((xname + "r2").Data());

   r1->SetAngles(phi1, theta1, 0);
   r2->SetAngles(phi2, theta2, 0);

   TGeoMatrix* m0 = new TGeoCombiTrans(*t0, *r0);
   TGeoMatrix* m1 = new TGeoCombiTrans(*t1, *r1);
   TGeoMatrix* m2 = new TGeoCombiTrans(*t2, *r2);

   TGeoCompositeShape *CS1 = new TGeoCompositeShape((xname + "CS1").Data(), new TGeoIntersection(T, B1, m0, m1));
   TGeoCompositeShape *cs = new TGeoCompositeShape((xname + "CS").Data(), new TGeoIntersection(CS1, B2, m0, m2));
   delete t0;
   delete t1;
   delete t2;
   delete r0;
   delete r1;
   delete r2;
   return cs;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks whether name2check is in (NameList) list

Bool_t TGDMLWrite::IsInList(NameList list, TString name2check)
{
   Bool_t isIN = list[name2check];
   return isIN;
}

////////////////////////////////////////////////////////////////////////////////
///NCNAME basic restrictions
///Replace "$" character with empty character etc.

TString TGDMLWrite::GenName(TString oldname)
{
   TString newname = oldname.ReplaceAll("$", "");
   newname = newname.ReplaceAll(" ", "_");
   // :, @, $, %, &, /, +, ,, ;, whitespace characters or different parenthesis
   newname = newname.ReplaceAll(":", "");
   newname = newname.ReplaceAll("@", "");
   newname = newname.ReplaceAll("%", "");
   newname = newname.ReplaceAll("&", "");
   newname = newname.ReplaceAll("/", "");
   newname = newname.ReplaceAll("+", "");
   newname = newname.ReplaceAll(";", "");
   newname = newname.ReplaceAll("{", "");
   newname = newname.ReplaceAll("}", "");
   newname = newname.ReplaceAll("(", "");
   newname = newname.ReplaceAll(")", "");
   newname = newname.ReplaceAll("[", "");
   newname = newname.ReplaceAll("]", "");
   newname = newname.ReplaceAll("_refl", "");
   //workaround if first letter is digit than replace it to "O" (ou character)
   TString fstLet = newname(0, 1);
   if (fstLet.IsDigit()) {
      newname = "O" + newname(1, newname.Length());
   }
   return newname;
}

////////////////////////////////////////////////////////////////////////////////
/// Important function which is responsible for naming volumes, solids and materials

TString TGDMLWrite::GenName(TString oldname, TString objPointer)
{
   TString newname = GenName(oldname);
   if (newname != oldname) {
      if (fgkMaxNameErr > fActNameErr) {
         Info("GenName",
              "WARNING! Name of the object was changed because it failed to comply with NCNAME xml datatype restrictions.");
      } else if ((fgkMaxNameErr == fActNameErr)) {
         Info("GenName",
              "WARNING! Probably more names are going to be changed to comply with NCNAME xml datatype restriction, but it will not be displayed on the screen.");
      }
      fActNameErr++;
   }
   TString nameIter;
   Int_t iter = 0;
   switch (fgNamingSpeed) {
      case kfastButUglySufix:
         newname = newname + "0x" + objPointer;
         break;
      case kelegantButSlow:
         //0 means not in the list
         iter = fNameList->fLstIter[newname];
         if (iter == 0) {
            nameIter = "";
         } else {
            nameIter = TString::Format("0x%i", iter);
         }
         fNameList->fLstIter[newname]++;
         newname = newname + nameIter;
         break;
      case kwithoutSufixNotUniq:
         //no change
         break;
   }
   //store the name (mapped to pointer)
   fNameList->fLst[objPointer] = newname;
   return newname;
}


////////////////////////////////////////////////////////////////////////////////
/// Method which tests whether solids can be processed

Bool_t TGDMLWrite::CanProcess(TObject *pointer)
{
   Bool_t isProcessed = kFALSE;
   isProcessed = pointer->TestBit(fgkProcBit);
   pointer->SetBit(fgkProcBit, kTRUE);
   return !(isProcessed);
}

////////////////////////////////////////////////////////////////////////////////
/// Method that retrieves axis and unit along which object is divided

TString TGDMLWrite::GetPattAxis(Int_t divAxis, const char * pattName, TString& unit)
{
   TString resaxis;
   unit = fDefault_lunit;
   switch (divAxis) {
      case 1:
         if (strcmp(pattName, "TGeoPatternX") == 0) {
            return "kXAxis";
         } else if (strcmp(pattName, "TGeoPatternCylR") == 0) {
            return "kRho";
         }
         break;
      case 2:
         if (strcmp(pattName, "TGeoPatternY") == 0) {
            return "kYAxis";
         } else if (strcmp(pattName, "TGeoPatternCylPhi") == 0) {
            unit = "deg";
            return "kPhi";
         }
         break;
      case 3:
         if (strcmp(pattName, "TGeoPatternZ") == 0) {
            return "kZAxis";
         }
         break;
      default:
         return "kUndefined";
         break;
   }
   return "kUndefined";
}

////////////////////////////////////////////////////////////////////////////////
/// Check for null parameter to skip the NULL objects

Bool_t TGDMLWrite::IsNullParam(Double_t parValue, TString parName, TString objName)
{
   if (parValue == 0.) {
      Info("IsNullParam", "ERROR! %s is NULL due to %s = %.12g, Volume based on this shape will be skipped",
           objName.Data(),
           parName.Data(),
           parValue);
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Unsetting bits that were changed in gGeoManager during export so that export
/// can be run more times with the same instance of gGeoManager.

void TGDMLWrite::UnsetTemporaryBits(TGeoManager * geoMng)
{
   TIter next(geoMng->GetListOfVolumes());
   TGeoVolume *vol;
   while ((vol = (TGeoVolume *)next())) {
      ((TObject *)vol->GetShape())->SetBit(fgkProcBit, kFALSE);
      vol->SetAttBit(fgkProcBitVol, kFALSE);
   }

}


////////////////////////////////////////////////////////////////////////////////
//
// Backwards compatibility for old DD4hep version (to be removed in the future)
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility (to be removed in the future): Wrapper to only selectively write one branch
void TGDMLWrite::WriteGDMLfile(TGeoManager * geomanager, TGeoVolume* volume, const char* filename, TString option)
{
  TList materials, volumes, nodes;
  MaterialExtractor extract;
  if ( !volume )   {
    Info("WriteGDMLfile", "Invalid Volume reference to extract GDML information!");
    return;
  }
  extract(volume);
  for(TGeoMaterial* m : extract.materials)
    materials.Add(m);
  fTopVolumeName = volume->GetName();
  fSurfaceList.clear();
  fVolumeList.clear();
  fNodeList.clear();
  WriteGDMLfile(geomanager, volume, &materials, filename, option);
  materials.Clear("nodelete");
  volumes.Clear("nodelete");
  nodes.Clear("nodelete");
}


////////////////////////////////////////////////////////////////////////////////
/// Wrapper of all exporting methods
/// Creates blank GDML file and fills it with gGeoManager structure converted
/// to GDML structure of xml nodes

void TGDMLWrite::WriteGDMLfile(TGeoManager * geomanager,
                               TGeoVolume* volume,
                               TList* materialsLst,
                               const char* filename,
                               TString option)
{
   //option processing
   option.ToLower();
   if (option.Contains("g")) {
      SetG4Compatibility(kTRUE);
      Info("WriteGDMLfile", "Geant4 compatibility mode set");
   } else {
      SetG4Compatibility(kFALSE);
   }
   if (option.Contains("f")) {
      SetNamingSpeed(kfastButUglySufix);
      Info("WriteGDMLfile", "Fast naming convention with pointer suffix set");
   } else if (option.Contains("n")) {
      SetNamingSpeed(kwithoutSufixNotUniq);
      Info("WriteGDMLfile", "Naming without prefix set - be careful uniqness of name is not ensured");
   } else {
      SetNamingSpeed(kelegantButSlow);
      Info("WriteGDMLfile", "Potentially slow with incremental suffix naming convention set");
   }

   //local variables
   Int_t outputLayout = 1;
   const char * krootNodeName = "gdml";
   const char * knsRefGeneral = "http://www.w3.org/2001/XMLSchema-instance";
   const char * knsNameGeneral = "xsi";
   const char * knsRefGdml = "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd";
   const char * knsNameGdml = "xsi:noNamespaceSchemaLocation";

   // First create engine
   fGdmlE = new TXMLEngine;
   fGdmlE->SetSkipComments(kTRUE);

   //create blank GDML file
   fGdmlFile = fGdmlE->NewDoc();

   //create root node and add it to blank GDML file
   XMLNodePointer_t rootNode = fGdmlE->NewChild(nullptr, nullptr, krootNodeName, nullptr);
   fGdmlE->DocSetRootElement(fGdmlFile, rootNode);

   //add namespaces to root node
   fGdmlE->NewNS(rootNode, knsRefGeneral, knsNameGeneral);
   fGdmlE->NewAttr(rootNode, nullptr, knsNameGdml, knsRefGdml);

   //initialize general lists and <define>, <solids>, <structure> nodes
   fIsotopeList  = new StructLst;
   fElementList  = new StructLst;

   fNameList     = new NameLst;

   fDefineNode = fGdmlE->NewChild(nullptr, nullptr, "define", nullptr);
   fSolidsNode = fGdmlE->NewChild(nullptr, nullptr, "solids", nullptr);
   fStructureNode = fGdmlE->NewChild(nullptr, nullptr, "structure", nullptr);
   //========================

   //initialize list of accepted patterns for divisions (in ExtractVolumes)
   fAccPatt   = new StructLst;
   fAccPatt->fLst["TGeoPatternX"] = kTRUE;
   fAccPatt->fLst["TGeoPatternY"] = kTRUE;
   fAccPatt->fLst["TGeoPatternZ"] = kTRUE;
   fAccPatt->fLst["TGeoPatternCylR"] = kTRUE;
   fAccPatt->fLst["TGeoPatternCylPhi"] = kTRUE;
   //========================

   //initialize list of rejected shapes for divisions (in ExtractVolumes)
   fRejShape     = new StructLst;
   //this shapes are rejected because, it is not possible to divide trd2
   //in Y axis and while only trd2 object is imported from GDML
   //it causes a problem when TGeoTrd1 is divided in Y axis
   fRejShape->fLst["TGeoTrd1"] = kTRUE;
   fRejShape->fLst["TGeoTrd2"] = kTRUE;
   //=========================

   //Initialize global counters
   fActNameErr = 0;
   fVolCnt = 0;
   fPhysVolCnt = 0;
   fSolCnt = 0;

   //calling main extraction functions (with measuring time)
   time_t startT, endT;
   startT = time(nullptr);
   ExtractMatrices(geomanager->GetListOfGDMLMatrices());
   ExtractConstants(geomanager);
   fMaterialsNode = ExtractMaterials(materialsLst);

   Info("WriteGDMLfile", "Extracting volumes");
   ExtractVolumes(volume);
   Info("WriteGDMLfile", "%i solids added", fSolCnt);
   Info("WriteGDMLfile", "%i volumes added", fVolCnt);
   Info("WriteGDMLfile", "%i physvolumes added", fPhysVolCnt);
   ExtractSkinSurfaces(geomanager->GetListOfSkinSurfaces());
   ExtractBorderSurfaces(geomanager->GetListOfBorderSurfaces());
   ExtractOpticalSurfaces(geomanager->GetListOfOpticalSurfaces());
   endT = time(nullptr);
   //<gdml>
   fGdmlE->AddChild(rootNode, fDefineNode);                 //  <define>...</define>
   fGdmlE->AddChild(rootNode, fMaterialsNode);              //  <materials>...</materials>
   fGdmlE->AddChild(rootNode, fSolidsNode);                 //  <solids>...</solids>
   fGdmlE->AddChild(rootNode, fStructureNode);              //  <structure>...</structure>
   fGdmlE->AddChild(rootNode, CreateSetupN(fTopVolumeName.Data())); //  <setup>...</setup>
   //</gdml>
   Double_t tdiffI = difftime(endT, startT);
   TString tdiffS = (tdiffI == 0 ? TString("< 1 s") : TString::Format("%.0lf s", tdiffI));
   Info("WriteGDMLfile", "Exporting time: %s", tdiffS.Data());
   //=========================

   //Saving document
   fGdmlE->SaveDoc(fGdmlFile, filename, outputLayout);
   Info("WriteGDMLfile", "File %s saved", filename);
   //cleaning
   fGdmlE->FreeDoc(fGdmlFile);
   //unset processing bits:
   UnsetTemporaryBits(geomanager);
   delete fGdmlE;
}




////////////////////////////////////////////////////////////////////////////////
/// Method extracting geometry structure recursively

void TGDMLWrite::ExtractVolumes(TGeoVolume* volume)
{
   XMLNodePointer_t volumeN, childN;
   TString volname, matname, solname, pattClsName, nodeVolNameBak;
   TGeoPatternFinder *pattFinder = nullptr;
   Bool_t isPattern = kFALSE;
   const TString fltPrecision = TString::Format("%%.%dg", fFltPrecision);

   //create the name for volume/assembly
   if (volume->IsTopVolume()) {
      //not needed a special function for generating name
      volname = volume->GetName();
      fTopVolumeName = volname;
      //register name to the pointer
      fNameList->fLst[TString::Format("%p", volume)] = volname;
   } else {
      volname = GenName(volume->GetName(), TString::Format("%p", volume));
   }

   //start to create main volume/assembly node
   if (volume->IsAssembly()) {
      volumeN = StartAssemblyN(volname);
   } else {
      //get reference material and add solid to <solids> + get name
      matname = fNameList->fLst[TString::Format("%p", volume->GetMaterial())];
      solname = ExtractSolid(volume->GetShape());
      //If solid is not supported or corrupted
      if (solname == "-1") {
         Info("ExtractVolumes", "ERROR! %s volume was not added, because solid is either not supported or corrupted",
              volname.Data());
         //set volume as missing volume
         fNameList->fLst[TString::Format("%p", volume)] = "missing_" + volname;
         return;
      }
      volumeN = StartVolumeN(volname, solname, matname);

      //divisionvol can't be in assembly
      pattFinder = volume->GetFinder();
      //if found pattern
      if (pattFinder) {
         pattClsName = TString::Format("%s", pattFinder->ClassName());
         TString shapeCls = TString::Format("%s", volume->GetShape()->ClassName());
         //if pattern in accepted pattern list and not in shape rejected list
         if ((fAccPatt->fLst[pattClsName] == kTRUE) &&
             (fRejShape->fLst[shapeCls] != kTRUE)) {
            isPattern = kTRUE;
         }
      }
   }
   //get all nodes in volume
   TObjArray *nodeLst = volume->GetNodes();
   TIter next(nodeLst);
   TGeoNode *geoNode;
   Int_t nCnt = 0;
   //loop through all nodes
   while ((geoNode = (TGeoNode *) next())) {
      //get volume of current node and if not processed then process it
      TGeoVolume * subvol = geoNode->GetVolume();
      if (subvol->TestAttBit(fgkProcBitVol) == kFALSE) {
         subvol->SetAttBit(fgkProcBitVol);
         ExtractVolumes(subvol);
      }

      //volume of this node has to exist because it was processed recursively
      TString nodevolname = fNameList->fLst[TString::Format("%p", geoNode->GetVolume())];
      if (nodevolname.Contains("missing_")) {
         continue;
      }
      if (nCnt == 0) { //save name of the first node for divisionvol
         nodeVolNameBak = nodevolname;
      }

      if (isPattern == kFALSE) {
         //create name for node
         TString nodename, posname, rotname;
         nodename = GenName(geoNode->GetName(), TString::Format("%p", geoNode));
         nodename = nodename + "in" + volname;

         //create name for position and clear rotation
         posname = nodename + "pos";
         rotname = "";

         //position
         const Double_t * pos = geoNode->GetMatrix()->GetTranslation();
         Xyz nodPos;
         nodPos.x = pos[0];
         nodPos.y = pos[1];
         nodPos.z = pos[2];
         childN = CreatePositionN(posname.Data(), nodPos);
         fGdmlE->AddChild(fDefineNode, childN); //adding node to <define> node
         //Deal with reflection
         XMLNodePointer_t scaleN = nullptr;
         Double_t lx, ly, lz;
         Double_t xangle = 0;
         Double_t zangle = 0;
         lx = geoNode->GetMatrix()->GetRotationMatrix()[0];
         ly = geoNode->GetMatrix()->GetRotationMatrix()[4];
         lz = geoNode->GetMatrix()->GetRotationMatrix()[8];
         if (geoNode->GetMatrix()->IsReflection()
             && TMath::Abs(lx) == 1 &&  TMath::Abs(ly) == 1 && TMath::Abs(lz) == 1) {
            scaleN = fGdmlE->NewChild(nullptr, nullptr, "scale", nullptr);
            fGdmlE->NewAttr(scaleN, nullptr, "name", (nodename + "scl").Data());
            fGdmlE->NewAttr(scaleN, nullptr, "x", TString::Format(fltPrecision.Data(), lx));
            fGdmlE->NewAttr(scaleN, nullptr, "y", TString::Format(fltPrecision.Data(), ly));
            fGdmlE->NewAttr(scaleN, nullptr, "z", TString::Format(fltPrecision.Data(), lz));
            //experimentally found out, that rotation should be updated like this
            if (lx == -1) {
               zangle = 180;
            }
            if (lz == -1) {
               xangle = 180;
            }
         }

         //rotation
         TGDMLWrite::Xyz lxyz = GetXYZangles(geoNode->GetMatrix()->GetRotationMatrix());
         lxyz.x -= xangle;
         lxyz.z -= zangle;
         if ((lxyz.x != 0.0) || (lxyz.y != 0.0) || (lxyz.z != 0.0)) {
            rotname = nodename + "rot";
            childN = CreateRotationN(rotname.Data(), lxyz);
            fGdmlE->AddChild(fDefineNode, childN); //adding node to <define> node
         }

         //create physvol for main volume/assembly node
         childN = CreatePhysVolN(geoNode->GetName(), geoNode->GetNumber(), nodevolname.Data(), posname.Data(), rotname.Data(), scaleN);
         fGdmlE->AddChild(volumeN, childN);
      }
      nCnt++;
   }
   //create only one divisionvol node
   if (isPattern && pattFinder) {
      //retrieve attributes of division
      Int_t ndiv, divaxis;
      Double_t offset, width, xlo, xhi;
      TString axis, unit;

      ndiv = pattFinder->GetNdiv();
      width = pattFinder->GetStep();

      divaxis = pattFinder->GetDivAxis();
      volume->GetShape()->GetAxisRange(divaxis, xlo, xhi);

      //compute relative start (not positional)
      offset = pattFinder->GetStart() - xlo;
      axis = GetPattAxis(divaxis, pattClsName, unit);

      //create division node
      childN = CreateDivisionN(offset, width, ndiv, axis.Data(), unit.Data(), nodeVolNameBak.Data());
      fGdmlE->AddChild(volumeN, childN);
   }

   fVolCnt++;
   //add volume/assembly node into the <structure> node
   fGdmlE->AddChild(fStructureNode, volumeN);
}
