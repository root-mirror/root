// @(#)root/gdml:$Id$
// Author: Anton Pytel 15/9/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGDMLWRITE
#define ROOT_TGDMLWRITE

#include "TGeoMatrix.h"
#include "TGeoElement.h"
#include "TXMLEngine.h"
#include "TGeoVolume.h"
#include "TGeoParaboloid.h"
#include "TGeoSphere.h"
#include "TGeoArb8.h"
#include "TGeoCone.h"
#include "TGeoPara.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"
#include "TGeoPcon.h"
#include "TGeoTorus.h"
#include "TGeoPgon.h"
#include "TGeoXtru.h"
#include "TGeoEltu.h"
#include "TGeoHype.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoScaledShape.h"
#include "TGeoTessellated.h"
#include "TGeoManager.h"
#include "TGDMLMatrix.h"

#include <map>
#include <set>
#include <iostream>

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGDMLWrite - Class for exporting geometries From ROOT's gGeoManager    //
//    (instance of TGeoManager class) To GDML file. More about GDML       //
//    see http://gdml.web.cern.ch.                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoOpticalSurface;
class TGeoSkinSurface;
class TGeoBorderSurface;

class TGDMLWrite : public TObject {
public:
   TGDMLWrite();
   virtual ~TGDMLWrite();

   static void StartGDMLWriting(TGeoManager * geomanager, const char* filename, TString option) {
      //static function -
      //options:
      //  g - set by default - geant4 compatibility
      //  f,n - if none of this two is set then naming convention is
      //        with incremental suffix, if "f" then suffix is pointer
      //        if "n" then there is no suffix, but uniqness of names
      //        is not secured.
      TGDMLWrite *writer = new TGDMLWrite;
      writer->SetFltPrecision(TGeoManager::GetExportPrecision());
      writer->WriteGDMLfile(geomanager, filename, option);
      delete writer;
   }
   //wrapper of all main methods for extraction
   void WriteGDMLfile(TGeoManager * geomanager, const char* filename = "test.gdml", TString option = "");
   // Wrapper to only selectively write one branch of the volume hierarchy to file
   void WriteGDMLfile(TGeoManager * geomanager, TGeoNode* top_node, const char* filename = "test.gdml", TString option = "");

   enum ENamingType {
      kelegantButSlow = 0,
      kwithoutSufixNotUniq = 1,
      kfastButUglySufix = 2
   };
   void SetNamingSpeed(ENamingType naming);
   void SetG4Compatibility(Bool_t G4Compatible) {
      fgG4Compatibility = G4Compatible;
   };

private:
   struct Xyz {
      Double_t x;
      Double_t y;
      Double_t z;
   };

   typedef  std::set<const TGeoOpticalSurface*> SurfaceList;
   typedef  std::set<const TGeoVolume*> VolList;
   typedef  std::set<const TGeoNode*>   NodeList;
   typedef  std::map<TString, Bool_t>   NameList;
   typedef  std::map<TString, TString>  NameListS;
   typedef  std::map<TString, Int_t>    NameListI;
   typedef  std::map<TString, Float_t>  NameListF;
   struct StructLst {
      NameList fLst;
   };     //to store pointers
   struct NameLst {
      NameListS fLst;        //to map pointers with names
      NameListI fLstIter;    //to store all the iterators for repeating names
   };

   //General lists
   StructLst *fIsotopeList;   //list of isotopes
   StructLst *fElementList;   //list of elements
   StructLst *fAccPatt;       //list of accepted patterns for division
   StructLst *fRejShape;      //list of rejected shapes
   SurfaceList fSurfaceList;  //list of optical surfaces
   VolList     fVolumeList;   //list of volumes
   NodeList    fNodeList;     //list of placed volumes

   NameLst *fNameList; //list of names (pointer mapped)

   //Data members
   static TGDMLWrite *fgGDMLWrite;                         //pointer to gdml writer
   Int_t  fgNamingSpeed;                                   //input option for volume and solid naming
   Bool_t fgG4Compatibility;                               //input option for Geant4 compatibility
   XMLDocPointer_t  fGdmlFile;                             //pointer storing xml file
   TString fDefault_lunit;                                 //Default unit of length (depends on ROOT unit system)
   TString fTopVolumeName;                                 //name of top volume
   TXMLEngine *fGdmlE;                                     //xml engine pointer

   XMLNodePointer_t fDefineNode;                           //main <define> node...
   XMLNodePointer_t fMaterialsNode;                        //main <materials> node...
   XMLNodePointer_t fSolidsNode;                           //main <solids> node...
   XMLNodePointer_t fStructureNode;                        //main <structure> node...
   Int_t        fVolCnt;                                   //count of volumes
   Int_t        fPhysVolCnt;                               //count of physical volumes
   UInt_t       fActNameErr;                               //count of name errors
   UInt_t       fSolCnt;                                   //count of name solids
   UInt_t       fFltPrecision;                             //! floating point precision when writing

   static const UInt_t fgkProcBit    = BIT(14);    //14th bit is set when solid is processed
   static const UInt_t fgkProcBitVol = BIT(19);    //19th bit is set when volume is processed
   static const UInt_t fgkMaxNameErr = 5;          //maximum number of errors for naming

   //I. Methods processing the gGeoManager geometry object structure
   //1. Main methods to extract everything from ROOT gGeoManager
   XMLNodePointer_t ExtractMaterials(TList* materialsLst); //result <materials>...
   TString          ExtractSolid(TGeoShape* volShape);     //adds <shape> to <solids>
   void             ExtractVolumes(TGeoNode* topNode);    //result <volume> node...  + corresp. shape
   void             ExtractMatrices(TObjArray *matrices);  //adds <matrix> to <define>
   void             ExtractConstants(TGeoManager *geom);   //adds <constant> to <define>
   void             ExtractOpticalSurfaces(TObjArray *surfaces); //adds <opticalsurface> to <solids>
   void             ExtractSkinSurfaces(TObjArray *surfaces);    //adds <skinsurface> to <structure>
   void             ExtractBorderSurfaces(TObjArray *surfaces);  //adds <bordersurface> to <structure>

   // Combined implementation to extract GDML information from the geometry tree
   void WriteGDMLfile(TGeoManager * geomanager, TGeoNode* top_node, TList* materialsLst, const char* filename, TString option);

   //1.1 Materials sub methods - creating Nodes
   XMLNodePointer_t CreateAtomN(Double_t atom, const char * unit = "g/mole");
   XMLNodePointer_t CreateDN(Double_t density, const char * unit = "g/cm3");
   XMLNodePointer_t CreateFractionN(Double_t percentage, const char * refName);
   XMLNodePointer_t CreatePropertyN(TNamed const &property);

   XMLNodePointer_t CreateIsotopN(TGeoIsotope * isotope, const char * name);
   XMLNodePointer_t CreateElementN(TGeoElement * element, XMLNodePointer_t materials, const char * name);
   XMLNodePointer_t CreateMixtureN(TGeoMixture * mixture, XMLNodePointer_t materials, TString mname);
   XMLNodePointer_t CreateMaterialN(TGeoMaterial * material, TString mname);


   //1.2 Solids sub methods
   XMLNodePointer_t ChooseObject(TGeoShape *geoShape);
   XMLNodePointer_t CreateZplaneN(Double_t z, Double_t rmin, Double_t rmax);

   XMLNodePointer_t CreateBoxN(TGeoBBox * geoShape);
   XMLNodePointer_t CreateParaboloidN(TGeoParaboloid * geoShape);
   XMLNodePointer_t CreateSphereN(TGeoSphere * geoShape);
   XMLNodePointer_t CreateArb8N(TGeoArb8 * geoShape);
   XMLNodePointer_t CreateConeN(TGeoConeSeg * geoShape);
   XMLNodePointer_t CreateConeN(TGeoCone * geoShape);
   XMLNodePointer_t CreateParaN(TGeoPara * geoShape);
   XMLNodePointer_t CreateTrapN(TGeoTrap * geoShape);
   XMLNodePointer_t CreateTwistedTrapN(TGeoGtra * geoShape);
   XMLNodePointer_t CreateTrdN(TGeoTrd1 * geoShape);
   XMLNodePointer_t CreateTrdN(TGeoTrd2 * geoShape);
   XMLNodePointer_t CreateTubeN(TGeoTubeSeg * geoShape);
   XMLNodePointer_t CreateCutTubeN(TGeoCtub * geoShape);
   XMLNodePointer_t CreateTubeN(TGeoTube * geoShape);
   XMLNodePointer_t CreatePolyconeN(TGeoPcon * geoShape);
   XMLNodePointer_t CreateTorusN(TGeoTorus * geoShape);
   XMLNodePointer_t CreatePolyhedraN(TGeoPgon * geoShape);
   XMLNodePointer_t CreateEltubeN(TGeoEltu * geoShape);
   XMLNodePointer_t CreateHypeN(TGeoHype * geoShape);
   XMLNodePointer_t CreateXtrusionN(TGeoXtru * geoShape);
   XMLNodePointer_t CreateEllipsoidN(TGeoCompositeShape * geoShape, TString elName);
   XMLNodePointer_t CreateElConeN(TGeoScaledShape * geoShape);
   XMLNodePointer_t CreateTessellatedN(TGeoTessellated * geoShape);
   XMLNodePointer_t CreateOpticalSurfaceN(TGeoOpticalSurface * geoSurf);
   XMLNodePointer_t CreateSkinSurfaceN(TGeoSkinSurface * geoSurf);
   XMLNodePointer_t CreateBorderSurfaceN(TGeoBorderSurface * geoSurf);

   XMLNodePointer_t CreateCommonBoolN(TGeoCompositeShape *geoShape);

   //1.3 Volume sub methods
   XMLNodePointer_t CreatePhysVolN(const char * name, Int_t copyno, const char * volref, const char * posref, const char * rotref, XMLNodePointer_t scaleN);
   XMLNodePointer_t CreateDivisionN(Double_t offset, Double_t width, Int_t number, const char * axis, const char * unit, const char * volref);

   XMLNodePointer_t CreateSetupN(const char * topVolName , const char * name = "default", const char * version = "1.0");
   XMLNodePointer_t StartVolumeN(const char * name, const char * solid, const char * material);
   XMLNodePointer_t StartAssemblyN(const char * name);


   //II. Utility methods
   Xyz GetXYZangles(const Double_t * rotationMatrix);
   //nodes to create position, rotation and similar types first-position/rotation...
   XMLNodePointer_t CreatePositionN(const char * name, Xyz position, const char * type = "position", const char * unit = "cm");
   XMLNodePointer_t CreateRotationN(const char * name, Xyz rotation, const char * type = "rotation", const char * unit = "deg");
   XMLNodePointer_t CreateMatrixN(TGDMLMatrix const *matrix);
   XMLNodePointer_t CreateConstantN(const char *name, Double_t value);
   TGeoCompositeShape* CreateFakeCtub(TGeoCtub * geoShape);  //create fake cut tube as intersection

   //check name (2nd parameter) whether it is in the list (1st parameter)
   Bool_t IsInList(NameList list, TString name2check);
   TString GenName(TString oldname);
   TString GenName(TString oldname, TString objPointer);
   Bool_t CanProcess(TObject *pointer);
   TString GetPattAxis(Int_t divAxis, const char * pattName, TString& unit);
   Bool_t IsNullParam(Double_t parValue, TString parName, TString objName);
   void UnsetTemporaryBits(TGeoManager * geoMng);
   UInt_t GetFltPrecision() const { return fFltPrecision; }
   void SetFltPrecision(UInt_t prec) { fFltPrecision = prec; }

   ////////////////////////////////////////////////////////////////////////////////
   //
   // Backwards compatibility for old DD4hep version (to be removed in the future)
   //
   ////////////////////////////////////////////////////////////////////////////////
public:
   // Backwards compatibility (to be removed in the future): Wrapper to only selectively write one branch
   void WriteGDMLfile(TGeoManager * geomanager, TGeoVolume* top_vol, const char* filename = "test.gdml", TString option = "");
private:
   // Backwards compatibility (to be removed in the future): Combined implementation to extract GDML information from the geometry tree
   void WriteGDMLfile(TGeoManager * geomanager, TGeoVolume* top_vol, TList* materialsLst, const char* filename, TString option);
   void ExtractVolumes(TGeoVolume* topVolume);    //result <volume> node...  + corresp. shape

   ClassDef(TGDMLWrite, 0)    //imports GDML using DOM and binds it to ROOT
};

#endif /* ROOT_TGDMLWRITE */
