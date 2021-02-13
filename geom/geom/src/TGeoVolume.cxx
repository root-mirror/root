// @(#)root/geom:$Id$
// Author: Andrei Gheata   30/05/02
// Divide(), CheckOverlaps() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoVolume
\ingroup Geometry_classes

TGeoVolume, TGeoVolumeMulti, TGeoVolumeAssembly are the volume classes

  Volumes are the basic objects used in building the geometrical hierarchy.
They represent unpositioned objects but store all information about the
placement of the other volumes they may contain. Therefore a volume can
be replicated several times in the geometry. In order to create a volume, one
has to put together a shape and a medium which are already defined. Volumes
have to be named by users at creation time. Every different name may represent a
an unique volume object, but may also represent more general a family (class)
of volume objects having the same shape type and medium, but possibly
different shape parameters. It is the user's task to provide different names
for different volume families in order to avoid ambiguities at tracking time.
A generic family rather than a single volume is created only in two cases :
when a generic shape is provided to the volume constructor or when a division
operation is applied. Each volume in the geometry stores an unique
ID corresponding to its family. In order to ease-up their creation, the manager
class is providing an API that allows making a shape and a volume in a single step.

  Volumes are objects that can be visualized, therefore having visibility,
colour, line and fill attributes that can be defined or modified any time after
the volume creation. It is advisable however to define these properties just
after the first creation of a volume namespace, since in case of volume families
any new member created by the modeler inherits these properties.

   In order to provide navigation features, volumes have to be able to find
the proper container of any point defined in the local reference frame. This
can be the volume itself, one of its positioned daughter volumes or none if
the point is actually outside. On the other hand, volumes have to provide also
other navigation methods such as finding the distances to its shape boundaries
or which daughter will be crossed first. The implementation of these features
is done at shape level, but the local mother-daughters management is handled
by volumes that builds additional optimisation structures upon geometry closure.
In order to have navigation features properly working one has to follow the
general rules for building a valid geometry (see TGeoManager class).

  Now let's make a simple volume representing a copper wire. We suppose that
a medium is already created (see TGeoMedium class on how to create media).
We will create a TUBE shape for our wire, having Rmin=0cm, Rmax=0.01cm
and a half-length dZ=1cm :

~~~ {.cpp}
  TGeoTube *tube = new TGeoTube("wire_tube", 0, 0.01, 1);
~~~

One may omit the name for the shape if no retrieving by name is further needed
during geometry building. The same shape can be shared by different volumes
having different names and materials. Now let's make the volume for our wire.
The prototype for volumes constructor looks like :

  TGeoVolume::TGeoVolume(const char *name, TGeoShape *shape, TGeoMedium *med)

Since TGeoTube derives from the base shape class, we can provide it to the volume
constructor :

~~~ {.cpp}
  TGeoVolume *wire_co = new TGeoVolume("WIRE_CO", tube, ptrCOPPER);
~~~

Do not bother to delete neither the media, shapes or volumes that you have
created since all will be automatically cleaned on exit by the manager class.
If we would have taken a look inside TGeoManager::MakeTube() method, we would
have been able to create our wire with a single line :

~~~ {.cpp}
  TGeoVolume *wire_co = gGeoManager->MakeTube("WIRE_CO", ptrCOPPER, 0, 0.01, 1);
~~~

The same applies for all primitive shapes, for which there can be found
corresponding MakeSHAPE() methods. Their usage is much more convenient unless
a shape has to be shared between more volumes. Let's make now an aluminium wire
having the same shape, supposing that we have created the copper wire with the
line above :

~~~ {.cpp}
  TGeoVolume *wire_al = new TGeoVolume("WIRE_AL", wire_co->GetShape(), ptrAL);
~~~

Now that we have learned how to create elementary volumes, let's see how we
can create a geometrical hierarchy.


### Positioning volumes

  When creating a volume one does not specify if this will contain or not other
volumes. Adding daughters to a volume implies creating those and adding them
one by one to the list of daughters. Since the volume has to know the position
of all its daughters, we will have to supply at the same time a geometrical
transformation with respect to its local reference frame for each of them.
The objects referencing a volume and a transformation are called NODES and
their creation is fully handled by the modeler. They represent the link
elements in the hierarchy of volumes. Nodes are unique and distinct geometrical
objects ONLY from their container point of view. Since volumes can be replicated
in the geometry, the same node may be found on different branches.

\image html geom_t_example.png

  An important observation is that volume objects are owned by the TGeoManager
class. This stores a list of all volumes in the geometry, that is cleaned
upon destruction.

  Let's consider positioning now our wire in the middle of a gas chamber. We
need first to define the gas chamber :

~~~ {.cpp}
  TGeoVolume *chamber = gGeoManager->MakeTube("CHAMBER", ptrGAS, 0, 1, 1);
~~~

Now we can put the wire inside :

~~~ {.cpp}
  chamber->AddNode(wire_co, 1);
~~~

If we inspect now the chamber volume in a browser, we will notice that it has
one daughter. Of course the gas has some container also, but let's keep it like
that for the sake of simplicity. The full prototype of AddNode() is :

~~~ {.cpp}
  TGeoVolume::AddNode(TGeoVolume *daughter, Int_t usernumber,
                      TGeoMatrix *matrix=gGeoIdentity)
~~~

Since we did not supplied the third argument, the wire will be positioned with
an identity transformation inside the chamber. One will notice that the inner
radii of the wire and chamber are both zero - therefore, aren't the two volumes
overlapping ? The answer is no, the modeler is even relaying on the fact that
any daughter is fully contained by its mother. On the other hand, neither of
the nodes positioned inside a volume should overlap with each other. We will
see that there are allowed some exceptions to those rules.

### Overlapping volumes

  Positioning volumes that does not overlap their neighbours nor extrude
their container is sometimes quite strong constraint. Some parts of the geometry
might overlap naturally, e.g. two crossing tubes. The modeller supports such
cases only if the overlapping nodes are declared by the user. In order to do
that, one should use TGeoVolume::AddNodeOverlap() instead of TGeoVolume::AddNode().
  When 2 or more positioned volumes are overlapping, not all of them have to
be declared so, but at least one. A point inside an overlapping region equally
belongs to all overlapping nodes, but the way these are defined can enforce
the modeler to give priorities.
  The general rule is that the deepest node in the hierarchy containing a point
have the highest priority. For the same geometry level, non-overlapping is
prioritised over overlapping. In order to illustrate this, we will consider
few examples. We will designate non-overlapping nodes as ONLY and the others
MANY as in GEANT3, where this concept was introduced:
  1. The part of a MANY node B extruding its container A will never be "seen"
during navigation, as if B was in fact the result of the intersection of A and B.
  2. If we have two nodes A (ONLY) and B (MANY) inside the same container, all
points in the overlapping region of A and B will be designated as belonging to A.
  3. If A an B in the above case were both MANY, points in the overlapping
part will be designated to the one defined first. Both nodes must have the
same medium.
  4. The slices of a divided MANY will be as well MANY.

One needs to know that navigation inside geometry parts MANY nodes is much
slower. Any overlapping part can be defined based on composite shapes - this
is always recommended.

### Replicating volumes

  What can we do if our chamber contains two identical wires instead of one ?
What if then we would need 1000 chambers in our detector ? Should we create
2000 wires and 1000 chamber volumes ? No, we will just need to replicate the
ones that we have already created.

~~~ {.cpp}
  chamber->AddNode(wire_co, 1, new TGeoTranslation(-0.2,0,0));
  chamber->AddNode(wire_co, 2, new TGeoTranslation(0.2,0,0));
~~~

  The 2 nodes that we have created inside chamber will both point to a wire_co
object, but will be completely distinct : WIRE_CO_1 and WIRE_CO_2. We will
want now to place symmetrically 1000 chambers on a pad, following a pattern
of 20 rows and 50 columns. One way to do this will be to replicate our chamber
by positioning it 1000 times in different positions of the pad. Unfortunately,
this is far from being the optimal way of doing what we want.
Imagine that we would like to find out which of the 1000 chambers is containing
a (x,y,z) point defined in the pad reference. You will never have to do that,
since the modeller will take care of it for you, but let's guess what it has
to do. The most simple algorithm will just loop over all daughters, convert
the point from mother to local reference and check if the current chamber
contains the point or not. This might be efficient for pads with few chambers,
but definitely not for 1000. Fortunately the modeler is smarter than that and
create for each volume some optimization structures called voxels (see Voxelization)
to minimize the penalty having too many daughters, but if you have 100 pads like
this in your geometry you will anyway loose a lot in your tracking performance.

  The way out when volumes can be arranged according to simple patterns is the
usage of divisions. We will describe them in detail later on. Let's think now
at a different situation : instead of 1000 chambers of the same type, we may
have several types of chambers. Let's say all chambers are cylindrical and have
a wire inside, but their dimensions are different. However, we would like all
to be represented by a single volume family, since they have the same properties.
*/

/** \class TGeoVolumeMulti
\ingroup Geometry_classes

Volume families

A volume family is represented by the class TGeoVolumeMulti. It represents
a class of volumes having the same shape type and each member will be
identified by the same name and volume ID. Any operation applied to a
TGeoVolume equally affects all volumes in that family. The creation of a
family is generally not a user task, but can be forced in particular cases:

~~~ {.cpp}
     TGeoManager::Volume(const char *vname, const char *shape, Int_t nmed);
~~~

where VNAME is the family name, NMED is the medium number and SHAPE is the
shape type that can be:

~~~ {.cpp}
  box    - for TGeoBBox
  trd1   - for TGeoTrd1
  trd2   - for TGeoTrd2
  trap   - for TGeoTrap
  gtra   - for TGeoGtra
  para   - for TGeoPara
  tube, tubs - for TGeoTube, TGeoTubeSeg
  cone, cons - for TGeoCone, TgeoCons
  eltu   - for TGeoEltu
  ctub   - for TGeoCtub
  pcon   - for TGeoPcon
  pgon   - for TGeoPgon
~~~

Volumes are then added to a given family upon adding the generic name as node
inside other volume:

~~~ {.cpp}
  TGeoVolume *box_family = gGeoManager->Volume("BOXES", "box", nmed);
  ...
  gGeoManager->Node("BOXES", Int_t copy_no, "mother_name",
                    Double_t x, Double_t y, Double_t z, Int_t rot_index,
                    Bool_t is_only, Double_t *upar, Int_t npar);
~~~

here:

~~~ {.cpp}
  BOXES   - name of the family of boxes
  copy_no - user node number for the created node
  mother_name - name of the volume to which we want to add the node
  x,y,z   - translation components
  rot_index   - indx of a rotation matrix in the list of matrices
  upar    - array of actual shape parameters
  npar    - number of parameters
~~~

The parameters order and number are the same as in the corresponding shape
constructors.

  Another particular case where volume families are used is when we want
that a volume positioned inside a container to match one ore more container
limits. Suppose we want to position the same box inside 2 different volumes
and we want the Z size to match the one of each container:

~~~ {.cpp}
  TGeoVolume *container1 = gGeoManager->MakeBox("C1", imed, 10,10,30);
  TGeoVolume *container2 = gGeoManager->MakeBox("C2", imed, 10,10,20);
  TGeoVolume *pvol       = gGeoManager->MakeBox("PVOL", jmed, 3,3,-1);
  container1->AddNode(pvol, 1);
  container2->AddNode(pvol, 1);
~~~

  Note that the third parameter of PVOL is negative, which does not make sense
as half-length on Z. This is interpreted as: when positioned, create a box
replacing all invalid parameters with the corresponding dimensions of the
container. This is also internally handled by the TGeoVolumeMulti class, which
does not need to be instantiated by users.

### Dividing volumes

  Volumes can be divided according a pattern. The most simple division can
be done along one axis, that can be: X, Y, Z, Phi, Rxy or Rxyz. Let's take
the most simple case: we would like to divide a box in N equal slices along X
coordinate, representing a new volume family. Supposing we already have created
the initial box, this can be done like:

~~~ {.cpp}
     TGeoVolume *slicex = box->Divide("SLICEX", 1, N);
~~~

where SLICE is the name of the new family representing all slices and 1 is the
slicing axis. The meaning of the axis index is the following: for all volumes
having shapes like box, trd1, trd2, trap, gtra or para - 1,2,3 means X,Y,Z; for
tube, tubs, cone, cons - 1 means Rxy, 2 means phi and 3 means Z; for pcon and
pgon - 2 means phi and 3 means Z; for spheres 1 means R and 2 means phi.
  In fact, the division operation has the same effect as positioning volumes
in a given order inside the divided container - the advantage being that the
navigation in such a structure is much faster. When a volume is divided, a
volume family corresponding to the slices is created. In case all slices can
be represented by a single shape, only one volume is added to the family and
positioned N times inside the divided volume, otherwise, each slice will be
represented by a distinct volume in the family.
  Divisions can be also performed in a given range of one axis. For that, one
have to specify also the starting coordinate value and the step:

~~~ {.cpp}
     TGeoVolume *slicex = box->Divide("SLICEX", 1, N, start, step);
~~~

A check is always done on the resulting division range : if not fitting into
the container limits, an error message is posted. If we will browse the divided
volume we will notice that it will contain N nodes starting with index 1 upto
N. The first one has the lower X limit at START position, while the last one
will have the upper X limit at START+N*STEP. The resulting slices cannot
be positioned inside an other volume (they are by default positioned inside the
divided one) but can be further divided and may contain other volumes:

~~~ {.cpp}
     TGeoVolume *slicey = slicex->Divide("SLICEY", 2, N1);
     slicey->AddNode(other_vol, index, some_matrix);
~~~

  When doing that, we have to remember that SLICEY represents a family, therefore
all members of the family will be divided on Y and the other volume will be
added as node inside all.
  In the example above all the resulting slices had the same shape as the
divided volume (box). This is not always the case. For instance, dividing a
volume with TUBE shape on PHI axis will create equal slices having TUBESEG
shape. Other divisions can also create slices having shapes with different
dimensions, e.g. the division of a TRD1 volume on Z.
  When positioning volumes inside slices, one can do it using the generic
volume family (e.g. slicey). This should be done as if the coordinate system
of the generic slice was the same as the one of the divided volume. The generic
slice in case of PHI division is centered with respect to X axis. If the
family contains slices of different sizes, any volume positioned inside should
fit into the smallest one.
   Examples for specific divisions according to shape types can be found inside
shape classes.

~~~ {.cpp}
       TGeoVolume::Divide(N, Xmin, Xmax, "X");
~~~

  The GEANT3 option MANY is supported by TGeoVolumeOverlap class. An overlapping
volume is in fact a virtual container that does not represent a physical object.
It contains a list of nodes that are not its daughters but that must be checked
always before the container itself. This list must be defined by users and it
is checked and resolved in a priority order. Note that the feature is non-standard
to geometrical modelers and it was introduced just to support conversions of
GEANT3 geometries, therefore its extensive usage should be avoided.
*/

/** \class TGeoVolumeAssembly
\ingroup Geometry_classes

Volume assemblies

Assemblies a volumes that have neither a shape or a material/medium. Assemblies
behave exactly like normal volumes grouping several daughters together, but
the daughters can never extrude the assembly since this has no shape. However,
a bounding box and a voxelization structure are built for assemblies as for
normal volumes, so that navigation is still optimized. Assemblies are useful
for grouping hierarchically volumes which are otherwise defined in a flat
manner, but also to avoid clashes between container shapes.
To define an assembly one should just input a name, then start adding other
volumes (or volume assemblies) as content.
*/

#include <fstream>
#include <iomanip>

#include "TString.h"
#include "TBuffer.h"
#include "TBrowser.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TMap.h"
#include "TFile.h"
#include "TKey.h"

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TVirtualGeoPainter.h"
#include "TGeoVolume.h"
#include "TGeoShapeAssembly.h"
#include "TGeoScaledShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoVoxelFinder.h"
#include "TGeoExtension.h"

ClassImp(TGeoVolume);

TGeoMedium *TGeoVolume::fgDummyMedium = 0;

////////////////////////////////////////////////////////////////////////////////
/// Create a dummy medium

void TGeoVolume::CreateDummyMedium()
{
   if (fgDummyMedium) return;
   fgDummyMedium = new TGeoMedium();
   fgDummyMedium->SetName("dummy");
   TGeoMaterial *dummyMaterial = new TGeoMaterial();
   dummyMaterial->SetName("dummy");
   fgDummyMedium->SetMaterial(dummyMaterial);
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolume::ClearThreadData() const
{
   if (fFinder) fFinder->ClearThreadData();
   if (fShape)  fShape->ClearThreadData();
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolume::CreateThreadData(Int_t nthreads)
{
   if (fFinder) fFinder->CreateThreadData(nthreads);
   if (fShape)  fShape->CreateThreadData(nthreads);
}

////////////////////////////////////////////////////////////////////////////////

TGeoMedium *TGeoVolume::DummyMedium()
{
   return fgDummyMedium;
}

////////////////////////////////////////////////////////////////////////////////
/// dummy constructor

TGeoVolume::TGeoVolume()
{
   fNodes    = 0;
   fShape    = 0;
   fMedium   = 0;
   fFinder   = 0;
   fVoxels   = 0;
   fGeoManager = gGeoManager;
   fField    = 0;
   fOption   = "";
   fNumber   = 0;
   fNtotal   = 0;
   fRefCount = 0;
   fUserExtension = 0;
   fFWExtension = 0;
   TObject::ResetBit(kVolumeImportNodes);
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoVolume::TGeoVolume(const char *name, const TGeoShape *shape, const TGeoMedium *med)
           :TNamed(name, "")
{
   fName = fName.Strip();
   fNodes    = 0;
   fShape    = (TGeoShape*)shape;
   if (fShape) {
      if (fShape->TestShapeBit(TGeoShape::kGeoBad)) {
         Warning("Ctor", "volume %s has invalid shape", name);
      }
      if (!fShape->IsValid()) {
         Fatal("ctor", "Shape of volume %s invalid. Aborting!", fName.Data());
      }
   }
   fMedium   = (TGeoMedium*)med;
   if (fMedium && fMedium->GetMaterial()) fMedium->GetMaterial()->SetUsed();
   fFinder   = 0;
   fVoxels   = 0;
   fGeoManager = gGeoManager;
   fField    = 0;
   fOption   = "";
   fNumber   = 0;
   fNtotal   = 0;
   fRefCount = 0;
   fUserExtension = 0;
   fFWExtension = 0;
   if (fGeoManager) fNumber = fGeoManager->AddVolume(this);
   TObject::ResetBit(kVolumeImportNodes);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVolume::~TGeoVolume()
{
   if (fNodes) {
      if (!TObject::TestBit(kVolumeImportNodes)) {
         fNodes->Delete();
      }
      delete fNodes;
   }
   if (fFinder && !TObject::TestBit(kVolumeImportNodes | kVolumeClone) ) delete fFinder;
   if (fVoxels) delete fVoxels;
   if (fUserExtension) {fUserExtension->Release(); fUserExtension=0;}
   if (fFWExtension) {fFWExtension->Release(); fFWExtension=0;}
}

////////////////////////////////////////////////////////////////////////////////
/// How to browse a volume

void TGeoVolume::Browse(TBrowser *b)
{
   if (!b) return;

//   if (!GetNdaughters()) b->Add(this, GetName(), IsVisible());
   TGeoVolume *daughter;
   TString title;
   for (Int_t i=0; i<GetNdaughters(); i++) {
      daughter = GetNode(i)->GetVolume();
      if(daughter->GetTitle()[0]) {
         if (daughter->IsAssembly()) title.TString::Format("Assembly with %d daughter(s)",
                                                daughter->GetNdaughters());
         else if (daughter->GetFinder()) {
            TString s1 = daughter->GetFinder()->ClassName();
            s1.ReplaceAll("TGeoPattern","");
            title.TString::Format("Volume having %s shape divided in %d %s slices",
                       daughter->GetShape()->ClassName(),daughter->GetNdaughters(), s1.Data());

         } else title.TString::Format("Volume with %s shape having %d daughter(s)",
                         daughter->GetShape()->ClassName(),daughter->GetNdaughters());
         daughter->SetTitle(title.Data());
      }
      b->Add(daughter, daughter->GetName(), daughter->IsVisible());
//      if (IsVisDaughters())
//      b->AddCheckBox(daughter, daughter->IsVisible());
//      else
//         b->AddCheckBox(daughter, kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the capacity of this [cm^3] as the capacity of its shape.
/// In case of assemblies, the capacity is computed as the sum of daughter's capacities.

Double_t TGeoVolume::Capacity() const
{
   if (!IsAssembly()) return fShape->Capacity();
   Double_t capacity = 0.0;
   Int_t nd = GetNdaughters();
   Int_t i;
   for (i=0; i<nd; i++) capacity += GetNode(i)->GetVolume()->Capacity();
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Shoot nrays with random directions from starting point (startx, starty, startz)
/// in the reference frame of this volume. Track each ray until exiting geometry, then
/// shoot backwards from exiting point and compare boundary crossing points.

void TGeoVolume::CheckGeometry(Int_t nrays, Double_t startx, Double_t starty, Double_t startz) const
{
   TGeoVolume *old_vol = fGeoManager->GetTopVolume();
   if (old_vol!=this) fGeoManager->SetTopVolume((TGeoVolume*)this);
   else old_vol=0;
   fGeoManager->GetTopVolume()->Draw();
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   painter->CheckGeometry(nrays, startx, starty, startz);
}

////////////////////////////////////////////////////////////////////////////////
/// Overlap checking tool. Check for illegal overlaps within a limit OVLP.
/// Use option="s[number]" to force overlap checking by sampling volume with
/// [number] points.
///
/// Ex:
/// ~~~ {.cpp}
///     myVol->CheckOverlaps(0.01, "s10000000"); // shoot 10000000 points
///     myVol->CheckOverlaps(0.01, "s"); // shoot the default value of 1e6 points
/// ~~~

void TGeoVolume::CheckOverlaps(Double_t ovlp, Option_t *option) const
{
   if (!GetNdaughters() || fFinder) return;
   Bool_t sampling = kFALSE;
   TString opt(option);
   opt.ToLower();
   if (opt.Contains("s")) sampling = kTRUE;
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   if (!sampling) fGeoManager->SetNsegments(80);
   if (!fGeoManager->IsCheckingOverlaps()) {
      fGeoManager->ClearOverlaps();
//      Info("CheckOverlaps", "=== Checking overlaps for volume %s ===\n", GetName());
   }
   painter->CheckOverlaps(this, ovlp, option);
//   if (sampling) return;
   if (!fGeoManager->IsCheckingOverlaps()) {
      fGeoManager->SortOverlaps();
      TObjArray *overlaps = fGeoManager->GetListOfOverlaps();
      Int_t novlps = overlaps->GetEntriesFast();
      TNamed *obj;
      TString name;
      for (Int_t i=0; i<novlps; i++) {
         obj = (TNamed*)overlaps->At(i);
         if (novlps<1000) name = TString::Format("ov%03d", i);
         else             name = TString::Format("ov%06d", i);
         obj->SetName(name);
      }
      if (novlps) Info("CheckOverlaps", "Number of illegal overlaps/extrusions for volume %s: %d\n", GetName(), novlps);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tests for checking the shape navigation algorithms. See TGeoShape::CheckShape()

void TGeoVolume::CheckShape(Int_t testNo, Int_t nsamples, Option_t *option)
{
   fShape->CheckShape(testNo,nsamples,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean data of the volume.

void TGeoVolume::CleanAll()
{
   ClearNodes();
   ClearShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the shape of this volume from the list held by the current manager.

void TGeoVolume::ClearShape()
{
   fGeoManager->ClearShape(fShape);
}

////////////////////////////////////////////////////////////////////////////////
/// check for negative parameters in shapes.

void TGeoVolume::CheckShapes()
{
   if (fShape->IsRunTimeShape()) {
      Error("CheckShapes", "volume %s has run-time shape", GetName());
      InspectShape();
      return;
   }
   if (!fNodes) return;
   Int_t nd=fNodes->GetEntriesFast();
   TGeoNode *node = 0;
   TGeoNode *new_node;
   const TGeoShape *shape = 0;
   TGeoVolume *old_vol;
   for (Int_t i=0; i<nd; i++) {
      node=(TGeoNode*)fNodes->At(i);
      // check if node has name
      if (!node->GetName()[0]) printf("Daughter %i of volume %s - NO NAME!!!\n",
                                           i, GetName());
      old_vol = node->GetVolume();
      shape = old_vol->GetShape();
      if (shape->IsRunTimeShape()) {
//         printf("   Node %s/%s has shape with negative parameters. \n",
//                 GetName(), node->GetName());
//         old_vol->InspectShape();
         // make a copy of the node
         new_node = node->MakeCopyNode();
         if (!new_node) {
            Fatal("CheckShapes", "Cannot make copy node for %s", node->GetName());
            return;
         }
         TGeoShape *new_shape = shape->GetMakeRuntimeShape(fShape, node->GetMatrix());
         if (!new_shape) {
            Error("CheckShapes","cannot resolve runtime shape for volume %s/%s\n",
                   GetName(),old_vol->GetName());
            continue;
         }
         TGeoVolume *new_volume = old_vol->MakeCopyVolume(new_shape);
//         printf(" new volume %s shape params :\n", new_volume->GetName());
//         new_volume->InspectShape();
         new_node->SetVolume(new_volume);
         // decouple the old node and put the new one instead
         fNodes->AddAt(new_node, i);
//         new_volume->CheckShapes();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Count total number of subnodes starting from this volume, nlevels down
///  - option = 0 (default) - count only once per volume
///  - option = 1           - count every time
///  - option = 2           - count volumes on visible branches
///  - option = 3           - return maximum level counted already with option = 0

Int_t TGeoVolume::CountNodes(Int_t nlevels, Int_t option)
{
   static Int_t maxlevel = 0;
   static Int_t nlev = 0;

   if (option<0 || option>3) option = 0;
   Int_t visopt = 0;
   Int_t nd = GetNdaughters();
   Bool_t last = (!nlevels || !nd)?kTRUE:kFALSE;
   switch (option) {
      case 0:
         if (fNtotal) return fNtotal;
      case 1:
         fNtotal = 1;
         break;
      case 2:
         visopt = fGeoManager->GetVisOption();
         if (!IsVisDaughters()) last = kTRUE;
         switch (visopt) {
            case TVirtualGeoPainter::kGeoVisDefault:
               fNtotal = (IsVisible())?1:0;
               break;
            case TVirtualGeoPainter::kGeoVisLeaves:
               fNtotal = (IsVisible() && last)?1:0;
         }
         if (!IsVisibleDaughters()) return fNtotal;
         break;
      case 3:
         return maxlevel;
   }
   if (last) return fNtotal;
   if (gGeoManager->GetTopVolume() == this) {
      maxlevel=0;
      nlev = 0;
   }
   if (nlev>maxlevel) maxlevel = nlev;
   TGeoNode *node;
   TGeoVolume *vol;
   nlev++;
   for (Int_t i=0; i<nd; i++) {
      node = GetNode(i);
      vol = node->GetVolume();
      fNtotal += vol->CountNodes(nlevels-1, option);
   }
   nlev--;
   return fNtotal;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if volume and all daughters are invisible.

Bool_t TGeoVolume::IsAllInvisible() const
{
   if (IsVisible()) return kFALSE;
   Int_t nd = GetNdaughters();
   for (Int_t i=0; i<nd; i++) if (GetNode(i)->GetVolume()->IsVisible()) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make volume and each of it daughters (in)visible.

void TGeoVolume::InvisibleAll(Bool_t flag)
{
   SetAttVisibility(!flag);
   Int_t nd = GetNdaughters();
   TObjArray *list = new TObjArray(nd+1);
   list->Add(this);
   TGeoVolume *vol;
   for (Int_t i=0; i<nd; i++) {
      vol = GetNode(i)->GetVolume();
      vol->SetAttVisibility(!flag);
      list->Add(vol);
   }
   TIter next(gROOT->GetListOfBrowsers());
   TBrowser *browser = 0;
   while ((browser=(TBrowser*)next())) {
      for (Int_t i=0; i<nd+1; i++) {
         vol = (TGeoVolume*)list->At(i);
         browser->CheckObjectItem(vol, !flag);
      }
      browser->Refresh();
   }
   delete list;
   fGeoManager->SetVisOption(4);
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if volume contains nodes

Bool_t TGeoVolume::IsFolder() const
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// check if the visibility and attributes are the default ones

Bool_t TGeoVolume::IsStyleDefault() const
{
   if (!IsVisible()) return kFALSE;
   if (GetLineColor() != gStyle->GetLineColor()) return kFALSE;
   if (GetLineStyle() != gStyle->GetLineStyle()) return kFALSE;
   if (GetLineWidth() != gStyle->GetLineWidth()) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// True if this is the top volume of the geometry

Bool_t TGeoVolume::IsTopVolume() const
{
   if (fGeoManager->GetTopVolume() == this) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the painter is currently ray-tracing the content of this volume.

Bool_t TGeoVolume::IsRaytracing() const
{
   return TGeoAtt::IsVisRaytrace();
}

////////////////////////////////////////////////////////////////////////////////
/// Inspect the material for this volume.

void TGeoVolume::InspectMaterial() const
{
   GetMaterial()->Print();
}

////////////////////////////////////////////////////////////////////////////////
/// Import a volume from a file.

TGeoVolume *TGeoVolume::Import(const char *filename, const char *name, Option_t * /*option*/)
{
   if (!gGeoManager) gGeoManager = new TGeoManager("geometry","");
   if (!filename) return 0;
   TGeoVolume *volume = 0;
   if (strstr(filename,".gdml")) {
   // import from a gdml file
   } else {
   // import from a root file
      TDirectory::TContext ctxt;
      TFile *f = TFile::Open(filename);
      if (!f || f->IsZombie()) {
         printf("Error: TGeoVolume::Import : Cannot open file %s\n", filename);
         return 0;
      }
      if (name && name[0]) {
         volume = (TGeoVolume*)f->Get(name);
      } else {
         TIter next(f->GetListOfKeys());
         TKey *key;
         while ((key = (TKey*)next())) {
            if (strcmp(key->GetClassName(),"TGeoVolume") != 0) continue;
            volume = (TGeoVolume*)key->ReadObj();
            break;
         }
      }
      delete f;
   }
   if (!volume) return NULL;
   volume->RegisterYourself();
   return volume;
}

////////////////////////////////////////////////////////////////////////////////
/// Export this volume to a file.
///
///  - Case 1: root file or root/xml file
///    if filename end with ".root". The key will be named name
///    if filename end with ".xml" a root/xml file is produced.
///
///  - Case 2: C++ script
///    if filename end with ".C"
///
///  - Case 3: gdml file
///    if filename end with ".gdml"
///
///  NOTE that to use this option, the PYTHONPATH must be defined like
///      export PYTHONPATH=$ROOTSYS/lib:$ROOTSYS/gdml
///

Int_t TGeoVolume::Export(const char *filename, const char *name, Option_t *option)
{
   TString sfile(filename);
   if (sfile.Contains(".C")) {
      //Save volume as a C++ script
      Info("Export","Exporting volume %s as C++ code", GetName());
      SaveAs(filename, "");
      return 1;
   }
   if (sfile.Contains(".gdml")) {
     //Save geometry as a gdml file
      Info("Export","Exporting %s as gdml code - not implemented yet", GetName());
      return 0;
   }
   if (sfile.Contains(".root") || sfile.Contains(".xml")) {
      //Save volume in a root file
      Info("Export","Exporting %s as root file.", GetName());
      TString opt(option);
      if (!opt.Length()) opt = "recreate";
      TFile *f = TFile::Open(filename,opt.Data());
      if (!f || f->IsZombie()) {
         Error("Export","Cannot open file");
         return 0;
      }
      TString keyname(name);
      if (keyname.IsNull()) keyname = GetName();
      Int_t nbytes = Write(keyname);
      delete f;
      return nbytes;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Actualize matrix of node indexed <inode>

void TGeoVolume::cd(Int_t inode) const
{
   if (fFinder) fFinder->cd(inode-fFinder->GetDivIndex());
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TGeoNode to the list of nodes. This is the usual method for adding
/// daughters inside the container volume.

void TGeoVolume::AddNode(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t * /*option*/)
{
   TGeoMatrix *matrix = mat;
   if (matrix==0) matrix = gGeoIdentity;
   else           matrix->RegisterYourself();
   if (!vol) {
      Error("AddNode", "Volume is NULL");
      return;
   }
   if (!vol->IsValid()) {
      Error("AddNode", "Won't add node with invalid shape");
      printf("### invalid volume was : %s\n", vol->GetName());
      return;
   }
   if (!fNodes) fNodes = new TObjArray();

   if (fFinder) {
      // volume already divided.
      Error("AddNode", "Cannot add node %s_%i into divided volume %s", vol->GetName(), copy_no, GetName());
      return;
   }

   TGeoNodeMatrix *node = 0;
   node = new TGeoNodeMatrix(vol, matrix);
   node->SetMotherVolume(this);
   fNodes->Add(node);
   TString name = TString::Format("%s_%d", vol->GetName(), copy_no);
//   if (fNodes->FindObject(name))
//      Warning("AddNode", "Volume %s : added node %s with same name", GetName(), name.Data());
   node->SetName(name);
   node->SetNumber(copy_no);
   fRefCount++;
   vol->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a division node to the list of nodes. The method is called by
/// TGeoVolume::Divide() for creating the division nodes.

void TGeoVolume::AddNodeOffset(TGeoVolume *vol, Int_t copy_no, Double_t offset, Option_t * /*option*/)
{
   if (!vol) {
      Error("AddNodeOffset", "invalid volume");
      return;
   }
   if (!vol->IsValid()) {
      Error("AddNode", "Won't add node with invalid shape");
      printf("### invalid volume was : %s\n", vol->GetName());
      return;
   }
   if (!fNodes) fNodes = new TObjArray();
   TGeoNode *node = new TGeoNodeOffset(vol, copy_no, offset);
   node->SetMotherVolume(this);
   fNodes->Add(node);
   TString name = TString::Format("%s_%d", vol->GetName(), copy_no+1);
   node->SetName(name);
   node->SetNumber(copy_no+1);
   vol->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TGeoNode to the list of nodes. This is the usual method for adding
/// daughters inside the container volume.

void TGeoVolume::AddNodeOverlap(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t *option)
{
   if (!vol) {
      Error("AddNodeOverlap", "Volume is NULL");
      return;
   }
   if (!vol->IsValid()) {
      Error("AddNodeOverlap", "Won't add node with invalid shape");
      printf("### invalid volume was : %s\n", vol->GetName());
      return;
   }
   if (vol->IsAssembly()) {
      Warning("AddNodeOverlap", "Declaring assembly %s as possibly overlapping inside %s not allowed. Using AddNode instead !",vol->GetName(),GetName());
      AddNode(vol, copy_no, mat, option);
      return;
   }
   TGeoMatrix *matrix = mat;
   if (matrix==0) matrix = gGeoIdentity;
   else           matrix->RegisterYourself();
   if (!fNodes) fNodes = new TObjArray();

   if (fFinder) {
      // volume already divided.
      Error("AddNodeOverlap", "Cannot add node %s_%i into divided volume %s", vol->GetName(), copy_no, GetName());
      return;
   }

   TGeoNodeMatrix *node = new TGeoNodeMatrix(vol, matrix);
   node->SetMotherVolume(this);
   fNodes->Add(node);
   TString name = TString::Format("%s_%d", vol->GetName(), copy_no);
   if (fNodes->FindObject(name))
      Warning("AddNode", "Volume %s : added node %s with same name", GetName(), name.Data());
   node->SetName(name);
   node->SetNumber(copy_no);
   node->SetOverlapping();
   if (vol->GetMedium() == fMedium)
      node->SetVirtual();
   vol->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Division a la G3. The volume will be divided along IAXIS (see shape classes), in NDIV
/// slices, from START with given STEP. The division volumes will have medium number NUMED.
/// If NUMED=0 they will get the medium number of the divided volume (this). If NDIV<=0,
/// all range of IAXIS will be divided and the resulting number of divisions will be centered on
/// IAXIS. If STEP<=0, the real STEP will be computed as the full range of IAXIS divided by NDIV.
/// Options (case insensitive):
///  - N  - divide all range in NDIV cells (same effect as STEP<=0) (GSDVN in G3)
///  - NX - divide range starting with START in NDIV cells          (GSDVN2 in G3)
///  - S  - divide all range with given STEP. NDIV is computed and divisions will be centered
///         in full range (same effect as NDIV<=0)                (GSDVS, GSDVT in G3)
///  - SX - same as DVS, but from START position.                   (GSDVS2, GSDVT2 in G3)

TGeoVolume *TGeoVolume::Divide(const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step, Int_t numed, Option_t *option)
{
   if (fFinder) {
   // volume already divided.
      Fatal("Divide","volume %s already divided", GetName());
      return 0;
   }
   TString opt(option);
   opt.ToLower();
   TString stype = fShape->ClassName();
   if (!fNodes) fNodes = new TObjArray();
   Double_t xlo, xhi, range;
   range = fShape->GetAxisRange(iaxis, xlo, xhi);
   // for phi divisions correct the range
   if (!strcmp(fShape->GetAxisName(iaxis), "PHI")) {
      if ((start-xlo)<-1E-3) start+=360.;
      if (TGeoShape::IsSameWithinTolerance(range,360)) {
         xlo = start;
         xhi = start+range;
      }
   }
   if (range <=0) {
      InspectShape();
      Fatal("Divide", "cannot divide volume %s (%s) on %s axis", GetName(), stype.Data(), fShape->GetAxisName(iaxis));
      return 0;
   }
   if (ndiv<=0 || opt.Contains("s")) {
      if (step<=0) {
         Fatal("Divide", "invalid division type for volume %s : ndiv=%i, step=%g", GetName(), ndiv, step);
         return 0;
      }
      if (opt.Contains("x")) {
         if ((xlo-start)>1E-3 || (xhi-start)<-1E-3) {
            Fatal("Divide", "invalid START=%g for division on axis %s of volume %s. Range is (%g, %g)",
                  start, fShape->GetAxisName(iaxis), GetName(), xlo, xhi);
            return 0;
         }
         xlo = start;
         range = xhi-xlo;
      }
      ndiv = Int_t((range+0.1*step)/step);
      Double_t ddx = range - ndiv*step;
      // always center the division in this case
      if (ddx>1E-3) Warning("Divide", "division of volume %s on %s axis (ndiv=%d) will be centered in the full range",
                            GetName(), fShape->GetAxisName(iaxis), ndiv);
      start = xlo + 0.5*ddx;
   }
   if (step<=0 || opt.Contains("n")) {
      if (opt.Contains("x")) {
         if ((xlo-start)>1E-3 || (xhi-start)<-1E-3) {
            Fatal("Divide", "invalid START=%g for division on axis %s of volume %s. Range is (%g, %g)",
                  start, fShape->GetAxisName(iaxis), GetName(), xlo, xhi);
            return 0;
         }
         xlo = start;
         range = xhi-xlo;
      }
      step  = range/ndiv;
      start = xlo;
   }

   Double_t end = start+ndiv*step;
   if (((start-xlo)<-1E-3) || ((end-xhi)>1E-3)) {
      Fatal("Divide", "division of volume %s on axis %s exceed range (%g, %g)",
            GetName(), fShape->GetAxisName(iaxis), xlo, xhi);
      return 0;
   }
   TGeoVolume *voldiv = fShape->Divide(this, divname, iaxis, ndiv, start, step);
   if (numed) {
      TGeoMedium *medium = fGeoManager->GetMedium(numed);
      if (!medium) {
         Fatal("Divide", "invalid medium number %d for division volume %s", numed, divname);
         return voldiv;
      }
      voldiv->SetMedium(medium);
      if (medium->GetMaterial()) medium->GetMaterial()->SetUsed();
   }
   return voldiv;
}

////////////////////////////////////////////////////////////////////////////////
/// compute the closest distance of approach from point px,py to this volume

Int_t TGeoVolume::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   TVirtualGeoPainter *painter = fGeoManager->GetPainter();
   Int_t dist = 9999;
   if (!painter) return dist;
   dist = painter->DistanceToPrimitiveVol(this, px, py);
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// draw top volume according to option

void TGeoVolume::Draw(Option_t *option)
{
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   fGeoManager->SetUserPaintVolume(this);
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   TGeoAtt::SetVisRaytrace(kFALSE);
   if (!IsVisContainers()) SetVisLeaves();
   if (option && option[0] > 0) {
      painter->DrawVolume(this, option);
   } else {
      painter->DrawVolume(this, gEnv->GetValue("Viewer3D.DefaultDrawOption",""));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// draw only this volume

void TGeoVolume::DrawOnly(Option_t *option)
{
   if (IsAssembly()) {
      Info("DrawOnly", "Volume assemblies do not support this option.");
      return;
   }
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   SetVisOnly();
   TGeoAtt::SetVisRaytrace(kFALSE);
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   if (option && option[0] > 0) {
      painter->DrawVolume(this, option);
   } else {
      painter->DrawVolume(this, gEnv->GetValue("Viewer3D.DefaultDrawOption",""));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform an extensive sampling to find which type of voxelization is
/// most efficient.

Bool_t TGeoVolume::OptimizeVoxels()
{
   printf("Optimizing volume %s ...\n", GetName());
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   return painter->TestVoxels(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Print volume info

void TGeoVolume::Print(Option_t *) const
{
   printf("== Volume: %s type %s positioned %d times\n", GetName(), ClassName(), fRefCount);
   InspectShape();
   InspectMaterial();
}

////////////////////////////////////////////////////////////////////////////////
/// paint volume

void TGeoVolume::Paint(Option_t *option)
{
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   painter->SetTopVolume(this);
//   painter->Paint(option);
   if (option && option[0] > 0) {
      painter->Paint(option);
   } else {
      painter->Paint(gEnv->GetValue("Viewer3D.DefaultDrawOption",""));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the voxels for this volume.

void TGeoVolume::PrintVoxels() const
{
   if (fVoxels) fVoxels->Print();
}

////////////////////////////////////////////////////////////////////////////////
/// Recreate the content of the other volume without pointer copying. Voxels are
/// ignored and supposed to be created in a later step via Voxelize.

void TGeoVolume::ReplayCreation(const TGeoVolume *other)
{
   Int_t nd = other->GetNdaughters();
   if (!nd) return;
   TGeoPatternFinder *finder = other->GetFinder();
   if (finder) {
      Int_t iaxis = finder->GetDivAxis();
      Int_t ndiv = finder->GetNdiv();
      Double_t start = finder->GetStart();
      Double_t step = finder->GetStep();
      Int_t numed = other->GetNode(0)->GetVolume()->GetMedium()->GetId();
      TGeoVolume *voldiv = Divide(other->GetNode(0)->GetVolume()->GetName(), iaxis, ndiv, start, step, numed);
      voldiv->ReplayCreation(other->GetNode(0)->GetVolume());
      return;
   }
   for (Int_t i=0; i<nd; i++) {
      TGeoNode *node = other->GetNode(i);
      if (node->IsOverlapping()) AddNodeOverlap(node->GetVolume(), node->GetNumber(), node->GetMatrix());
      else AddNode(node->GetVolume(), node->GetNumber(), node->GetMatrix());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// print nodes

void TGeoVolume::PrintNodes() const
{
   Int_t nd = GetNdaughters();
   for (Int_t i=0; i<nd; i++) {
      printf("%s\n", GetNode(i)->GetName());
      cd(i);
      GetNode(i)->GetMatrix()->Print();
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Generate a lego plot fot the top volume, according to option.

TH2F *TGeoVolume::LegoPlot(Int_t ntheta, Double_t themin, Double_t themax,
                            Int_t nphi,   Double_t phimin, Double_t phimax,
                            Double_t rmin, Double_t rmax, Option_t *option)
{
   TVirtualGeoPainter *p = fGeoManager->GetGeomPainter();
   TGeoVolume *old_vol = fGeoManager->GetTopVolume();
   if (old_vol!=this) fGeoManager->SetTopVolume(this);
   else old_vol=0;
   TH2F *hist = p->LegoPlot(ntheta, themin, themax, nphi, phimin, phimax, rmin, rmax, option);
   hist->Draw("lego1sph");
   return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// Register the volume and all materials/media/matrices/shapes to the manager.

void TGeoVolume::RegisterYourself(Option_t *option)
{
   if (fGeoManager->GetListOfVolumes()->FindObject(this)) return;
   // Register volume
   fGeoManager->AddVolume(this);
   // Register shape
   if (!fGeoManager->GetListOfShapes()->FindObject(fShape)) {
      if (fShape->IsComposite()) {
         TGeoCompositeShape *comp = (TGeoCompositeShape*)fShape;
         comp->RegisterYourself();
      } else {
         fGeoManager->AddShape(fShape);
      }
   }
   // Register medium/material
   if (fMedium && !fGeoManager->GetListOfMedia()->FindObject(fMedium)) {
      fGeoManager->GetListOfMedia()->Add(fMedium);
      if (!fGeoManager->GetListOfMaterials()->FindObject(fMedium->GetMaterial()))
         fGeoManager->AddMaterial(fMedium->GetMaterial());
   }
   // Register matrices for nodes.
   TGeoMatrix *matrix;
   TGeoNode *node;
   Int_t nd = GetNdaughters();
   Int_t i;
   for (i=0; i<nd; i++) {
      node = GetNode(i);
      matrix = node->GetMatrix();
      if (!matrix->IsRegistered()) matrix->RegisterYourself();
      else if (!fGeoManager->GetListOfMatrices()->FindObject(matrix)) {
         fGeoManager->GetListOfMatrices()->Add(matrix);
      }
   }
   // Call RegisterYourself recursively
   for (i=0; i<nd; i++) GetNode(i)->GetVolume()->RegisterYourself(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw random points in the bounding box of this volume.

void TGeoVolume::RandomPoints(Int_t npoints, Option_t *option)
{
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   TGeoVolume *old_vol = fGeoManager->GetTopVolume();
   if (old_vol!=this) fGeoManager->SetTopVolume(this);
   else old_vol=0;
   fGeoManager->RandomPoints(this, npoints, option);
   if (old_vol) fGeoManager->SetTopVolume(old_vol);
}

////////////////////////////////////////////////////////////////////////////////
/// Random raytracing method.

void TGeoVolume::RandomRays(Int_t nrays, Double_t startx, Double_t starty, Double_t startz, const char *target_vol, Bool_t check_norm)
{
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   TGeoVolume *old_vol = fGeoManager->GetTopVolume();
   if (old_vol!=this) fGeoManager->SetTopVolume(this);
   else old_vol=0;
   fGeoManager->RandomRays(nrays, startx, starty, startz, target_vol, check_norm);
   if (old_vol) fGeoManager->SetTopVolume(old_vol);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this volume with current settings and perform raytracing in the pad.

void TGeoVolume::Raytrace(Bool_t flag)
{
   TGeoAtt::SetVisRaytrace(kFALSE);
   if (gGeoManager != fGeoManager) gGeoManager = fGeoManager;
   TVirtualGeoPainter *painter = fGeoManager->GetGeomPainter();
   Bool_t drawn = (painter->GetDrawnVolume()==this)?kTRUE:kFALSE;
   if (!drawn) {
      painter->DrawVolume(this, "");
      TGeoAtt::SetVisRaytrace(flag);
      painter->ModifiedPad();
      return;
   }
   TGeoAtt::SetVisRaytrace(flag);
   painter->ModifiedPad();
}

////////////////////////////////////////////////////////////////////////////////
///  Save geometry having this as top volume as a C++ macro.

void TGeoVolume::SaveAs(const char *filename, Option_t *option) const
{
   if (!filename) return;
   std::ofstream out;
   out.open(filename, std::ios::out);
   if (out.bad()) {
      Error("SavePrimitive", "Bad file name: %s", filename);
      return;
   }
   if (fGeoManager->GetTopVolume() != this) fGeoManager->SetTopVolume((TGeoVolume*)this);

   TString fname(filename);
   Int_t ind = fname.Index(".");
   if (ind>0) fname.Remove(ind);
   out << "void "<<fname<<"() {" << std::endl;
   out << "   gSystem->Load(\"libGeom\");" << std::endl;
   const UInt_t prec = TGeoManager::GetExportPrecision();
   out << std::setprecision(prec);
   ((TGeoVolume*)this)->SavePrimitive(out,option);
   out << "}" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect user-defined extension to the volume. The volume "grabs" a copy, so
/// the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for user extensions and is guaranteed not
/// to be used by TGeo

void TGeoVolume::SetUserExtension(TGeoExtension *ext)
{
   if (fUserExtension) fUserExtension->Release();
   fUserExtension = 0;
   if (ext) fUserExtension = ext->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect framework defined extension to the volume. The volume "grabs" a copy,
/// so the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for the use by TGeo and the users should
///       NOT connect extensions using this method

void TGeoVolume::SetFWExtension(TGeoExtension *ext)
{
   if (fFWExtension) fFWExtension->Release();
   fFWExtension = 0;
   if (ext) fFWExtension = ext->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the user extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoVolume::GrabUserExtension() const
{
   if (fUserExtension) return fUserExtension->Grab();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the framework extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoVolume::GrabFWExtension() const
{
   if (fFWExtension) return fFWExtension->Grab();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVolume::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   Int_t i,icopy;
   Int_t nd = GetNdaughters();
   TGeoVolume *dvol;
   TGeoNode *dnode;
   TGeoMatrix *matrix;

   // check if we need to save shape/volume
   Bool_t mustDraw = kFALSE;
   if (fGeoManager->GetGeomPainter()->GetTopVolume()==this) mustDraw = kTRUE;
   if (!option[0]) {
      fGeoManager->SetAllIndex();
      out << "   new TGeoManager(\"" << fGeoManager->GetName() << "\", \"" << fGeoManager->GetTitle() << "\");" << std::endl << std::endl;
//      if (mustDraw) out << "   Bool_t mustDraw = kTRUE;" << std::endl;
//      else          out << "   Bool_t mustDraw = kFALSE;" << std::endl;
      out << "   Double_t dx,dy,dz;" << std::endl;
      out << "   Double_t dx1, dx2, dy1, dy2;" << std::endl;
      out << "   Double_t vert[20], par[20];" << std::endl;
      out << "   Double_t theta, phi, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2;" << std::endl;
      out << "   Double_t twist;" << std::endl;
      out << "   Double_t origin[3];" << std::endl;
      out << "   Double_t rmin, rmax, rmin1, rmax1, rmin2, rmax2;" << std::endl;
      out << "   Double_t r, rlo, rhi;" << std::endl;
      out << "   Double_t phi1, phi2;" << std::endl;
      out << "   Double_t a,b;" << std::endl;
      out << "   Double_t point[3], norm[3];" << std::endl;
      out << "   Double_t rin, stin, rout, stout;" << std::endl;
      out << "   Double_t thx, phx, thy, phy, thz, phz;" << std::endl;
      out << "   Double_t alpha, theta1, theta2, phi1, phi2, dphi;" << std::endl;
      out << "   Double_t tr[3], rot[9];" << std::endl;
      out << "   Double_t z, density, radl, absl, w;" << std::endl;
      out << "   Double_t lx,ly,lz,tx,ty,tz;" << std::endl;
      out << "   Double_t xvert[50], yvert[50];" << std::endl;
      out << "   Double_t zsect,x0,y0,scale0;" << std::endl;
      out << "   Int_t nel, numed, nz, nedges, nvert;" << std::endl;
      out << "   TGeoBoolNode *pBoolNode = 0;" << std::endl << std::endl;
      // first save materials/media
      out << "   // MATERIALS, MIXTURES AND TRACKING MEDIA" << std::endl;
      SavePrimitive(out, "m");
      // then, save matrices
      out << std::endl << "   // TRANSFORMATION MATRICES" << std::endl;
      SavePrimitive(out, "x");
      // save this volume and shape
      SavePrimitive(out, "s");
      out << std::endl << "   // SET TOP VOLUME OF GEOMETRY" << std::endl;
      out << "   gGeoManager->SetTopVolume(" << GetPointerName() << ");" << std::endl;
      // save daughters
      out << std::endl << "   // SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY" << std::endl;
      SavePrimitive(out, "d");
      out << std::endl << "   // CLOSE GEOMETRY" << std::endl;
      out << "   gGeoManager->CloseGeometry();" << std::endl;
      if (mustDraw) {
         if (!IsRaytracing()) out << "   gGeoManager->GetTopVolume()->Draw();" << std::endl;
         else                 out << "   gGeoManager->GetTopVolume()->Raytrace();" << std::endl;
      }
      return;
   }
   // check if we need to save shape/volume
   if (!strcmp(option, "s")) {
      // create the shape for this volume
      if (TestAttBit(TGeoAtt::kSavePrimitiveAtt)) return;
      if (!IsAssembly()) {
         fShape->SavePrimitive(out,option);
         out << "   // Volume: " << GetName() << std::endl;
         if (fMedium) out << "   " << GetPointerName() << " = new TGeoVolume(\"" << GetName() << "\"," << fShape->GetPointerName() << ", "<< fMedium->GetPointerName() << ");" << std::endl;
         else out << "   " << GetPointerName() << " = new TGeoVolume(\"" << GetName() << "\"," << fShape->GetPointerName() <<  ");" << std::endl;

      } else {
         out << "   // Assembly: " << GetName() << std::endl;
         out << "   " << GetPointerName() << " = new TGeoVolumeAssembly(\"" << GetName() << "\"" << ");" << std::endl;
      }
      if (fLineColor != 1) out << "   " << GetPointerName() << "->SetLineColor(" << fLineColor << ");" << std::endl;
      if (fLineWidth != 1) out << "   " << GetPointerName() << "->SetLineWidth(" << fLineWidth << ");" << std::endl;
      if (fLineStyle != 1) out << "   " << GetPointerName() << "->SetLineStyle(" << fLineStyle << ");" << std::endl;
      if (!IsVisible() && !IsAssembly()) out << "   " << GetPointerName() << "->SetVisibility(kFALSE);" << std::endl;
      if (!IsVisibleDaughters()) out << "   " << GetPointerName() << "->VisibleDaughters(kFALSE);" << std::endl;
      if (IsVisContainers()) out << "   " << GetPointerName() << "->SetVisContainers(kTRUE);" << std::endl;
      if (IsVisLeaves()) out << "   " << GetPointerName() << "->SetVisLeaves(kTRUE);" << std::endl;
      SetAttBit(TGeoAtt::kSavePrimitiveAtt);
   }
   // check if we need to save the media
   if (!strcmp(option, "m")) {
      if (fMedium) fMedium->SavePrimitive(out,option);
      for (i=0; i<nd; i++) {
         dvol = GetNode(i)->GetVolume();
         dvol->SavePrimitive(out,option);
      }
      return;
   }
   // check if we need to save the matrices
   if (!strcmp(option, "x")) {
      if (fFinder) {
         dvol = GetNode(0)->GetVolume();
         dvol->SavePrimitive(out,option);
         return;
      }
      for (i=0; i<nd; i++) {
         dnode = GetNode(i);
         matrix = dnode->GetMatrix();
         if (!matrix->IsIdentity()) matrix->SavePrimitive(out,option);
         dnode->GetVolume()->SavePrimitive(out,option);
      }
      return;
   }
   // check if we need to save volume daughters
   if (!strcmp(option, "d")) {
      if (!nd) return;
      if (TestAttBit(TGeoAtt::kSaveNodesAtt)) return;
      SetAttBit(TGeoAtt::kSaveNodesAtt);
      if (fFinder) {
         // volume divided: generate volume->Divide()
         dnode = GetNode(0);
         dvol = dnode->GetVolume();
         out << "   TGeoVolume *" << dvol->GetPointerName() << " = ";
         out << GetPointerName() << "->Divide(\"" << dvol->GetName() << "\", ";
         fFinder->SavePrimitive(out,option);
         if (fMedium != dvol->GetMedium()) {
            out << ", " << dvol->GetMedium()->GetId();
         }
         out << ");" << std::endl;
         dvol->SavePrimitive(out,"d");
         return;
      }
      for (i=0; i<nd; i++) {
         dnode = GetNode(i);
         dvol = dnode->GetVolume();
         dvol->SavePrimitive(out,"s");
         matrix = dnode->GetMatrix();
         icopy = dnode->GetNumber();
         // generate AddNode()
         out << "   " << GetPointerName() << "->AddNode";
         if (dnode->IsOverlapping()) out << "Overlap";
         out << "(" << dvol->GetPointerName() << ", " << icopy;
         if (!matrix->IsIdentity()) out << ", " << matrix->GetPointerName();
         out << ");" << std::endl;
      }
      // Recursive loop to daughters
      for (i=0; i<nd; i++) {
         dnode = GetNode(i);
         dvol = dnode->GetVolume();
         dvol->SavePrimitive(out,"d");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset SavePrimitive bits.

void TGeoVolume::UnmarkSaved()
{
   ResetAttBit(TGeoAtt::kSavePrimitiveAtt);
   ResetAttBit(TGeoAtt::kSaveNodesAtt);
   if (fShape) fShape->ResetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute mouse actions on this volume.

void TGeoVolume::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   TVirtualGeoPainter *painter = fGeoManager->GetPainter();
   if (!painter) return;
   painter->ExecuteVolumeEvent(this, event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// search a daughter inside the list of nodes

TGeoNode *TGeoVolume::FindNode(const char *name) const
{
   return ((TGeoNode*)fNodes->FindObject(name));
}

////////////////////////////////////////////////////////////////////////////////
/// Get the index of a daughter within check_list by providing the node pointer.

Int_t TGeoVolume::GetNodeIndex(const TGeoNode *node, Int_t *check_list, Int_t ncheck) const
{
   TGeoNode *current = 0;
   for (Int_t i=0; i<ncheck; i++) {
      current = (TGeoNode*)fNodes->At(check_list[i]);
      if (current==node) return check_list[i];
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// get index number for a given daughter

Int_t TGeoVolume::GetIndex(const TGeoNode *node) const
{
   TGeoNode *current = 0;
   Int_t nd = GetNdaughters();
   if (!nd) return -1;
   for (Int_t i=0; i<nd; i++) {
      current = (TGeoNode*)fNodes->At(i);
      if (current==node) return i;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get volume info for the browser.

char *TGeoVolume::GetObjectInfo(Int_t px, Int_t py) const
{
   TGeoVolume *vol = (TGeoVolume*)this;
   TVirtualGeoPainter *painter = fGeoManager->GetPainter();
   if (!painter) return 0;
   return (char*)painter->GetVolumeInfo(vol, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if cylindrical voxelization is optimal.

Bool_t TGeoVolume::GetOptimalVoxels() const
{
   Int_t nd = GetNdaughters();
   if (!nd) return kFALSE;
   Int_t id;
   Int_t ncyl = 0;
   TGeoNode *node;
   for (id=0; id<nd; id++) {
      node = (TGeoNode*)fNodes->At(id);
      ncyl += node->GetOptimalVoxels();
   }
   if (ncyl>(nd/2)) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Provide a pointer name containing uid.

char *TGeoVolume::GetPointerName() const
{
   static TString name;
   name = TString::Format("p%s_%lx", GetName(), (ULong_t)this);
   return (char*)name.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Getter for optimization structure.

TGeoVoxelFinder *TGeoVolume::GetVoxels() const
{
   if (fVoxels && !fVoxels->IsInvalid()) return fVoxels;
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Move perspective view focus to this volume

void TGeoVolume::GrabFocus()
{
   TVirtualGeoPainter *painter = fGeoManager->GetPainter();
   if (painter) painter->GrabFocus();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if the volume is an assembly or a scaled assembly.

Bool_t TGeoVolume::IsAssembly() const
{
  return fShape->IsAssembly();
}

////////////////////////////////////////////////////////////////////////////////
/// Clone this volume.
/// build a volume with same name, shape and medium

TGeoVolume *TGeoVolume::CloneVolume() const
{
   TGeoVolume *vol = new TGeoVolume(GetName(), fShape, fMedium);
   Int_t i;
   // copy volume attributes
   vol->SetTitle(GetTitle());
   vol->SetLineColor(GetLineColor());
   vol->SetLineStyle(GetLineStyle());
   vol->SetLineWidth(GetLineWidth());
   vol->SetFillColor(GetFillColor());
   vol->SetFillStyle(GetFillStyle());
   // copy other attributes
   Int_t nbits = 8*sizeof(UInt_t);
   for (i=0; i<nbits; i++)
      vol->SetAttBit(1<<i, TGeoAtt::TestAttBit(1<<i));
   for (i=14; i<24; i++)
      vol->SetBit(1<<i, TestBit(1<<i));

   // copy field
   vol->SetField(fField);
   // Set bits
   for (i=0; i<nbits; i++)
      vol->SetBit(1<<i, TObject::TestBit(1<<i));
   vol->SetBit(kVolumeClone);
   // copy nodes
//   CloneNodesAndConnect(vol);
   vol->MakeCopyNodes(this);
   // if volume is divided, copy finder
   vol->SetFinder(fFinder);
   // copy voxels
   TGeoVoxelFinder *voxels = 0;
   if (fVoxels) {
      voxels = new TGeoVoxelFinder(vol);
      vol->SetVoxelFinder(voxels);
   }
   // copy option, uid
   vol->SetOption(fOption);
   vol->SetNumber(fNumber);
   vol->SetNtotal(fNtotal);
   // copy extensions
   vol->SetUserExtension(fUserExtension);
   vol->SetFWExtension(fFWExtension);
   vol->SetOverlappingCandidate(IsOverlappingCandidate());
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone the array of nodes.

void TGeoVolume::CloneNodesAndConnect(TGeoVolume *newmother) const
{
   if (!fNodes) return;
   TGeoNode *node;
   Int_t nd = fNodes->GetEntriesFast();
   if (!nd) return;
   // create new list of nodes
   TObjArray *list = new TObjArray(nd);
   // attach it to new volume
   newmother->SetNodes(list);
//   ((TObject*)newmother)->SetBit(kVolumeImportNodes);
   for (Int_t i=0; i<nd; i++) {
      //create copies of nodes and add them to list
      node = GetNode(i)->MakeCopyNode();
      if (!node) {
         Fatal("CloneNodesAndConnect", "cannot make copy node");
         return;
      }
      node->SetMotherVolume(newmother);
      list->Add(node);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// make a new list of nodes and copy all nodes of other volume inside

void TGeoVolume::MakeCopyNodes(const TGeoVolume *other)
{
   Int_t nd = other->GetNdaughters();
   if (!nd) return;
   if (fNodes) {
      if (!TObject::TestBit(kVolumeImportNodes)) fNodes->Delete();
      delete fNodes;
   }
   fNodes = new TObjArray();
   for (Int_t i=0; i<nd; i++) fNodes->Add(other->GetNode(i));
   TObject::SetBit(kVolumeImportNodes);
}

////////////////////////////////////////////////////////////////////////////////
/// make a copy of this volume
/// build a volume with same name, shape and medium

TGeoVolume *TGeoVolume::MakeCopyVolume(TGeoShape *newshape)
{
   TGeoVolume *vol = new TGeoVolume(GetName(), newshape, fMedium);
   // copy volume attributes
   vol->SetVisibility(IsVisible());
   vol->SetLineColor(GetLineColor());
   vol->SetLineStyle(GetLineStyle());
   vol->SetLineWidth(GetLineWidth());
   vol->SetFillColor(GetFillColor());
   vol->SetFillStyle(GetFillStyle());
   // copy field
   vol->SetField(fField);
   // if divided, copy division object
   if (fFinder) {
//       Error("MakeCopyVolume", "volume %s divided", GetName());
      vol->SetFinder(fFinder);
   }
   // Copy extensions
   vol->SetUserExtension(fUserExtension);
   vol->SetFWExtension(fFWExtension);
   CloneNodesAndConnect(vol);
//   ((TObject*)vol)->SetBit(kVolumeImportNodes);
   ((TObject*)vol)->SetBit(kVolumeClone);
   vol->SetOverlappingCandidate(IsOverlappingCandidate());
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this volume which is reflected with respect to XY plane.

TGeoVolume *TGeoVolume::MakeReflectedVolume(const char *newname) const
{
   static TMap map(100);
   if (!fGeoManager->IsClosed()) {
      Error("MakeReflectedVolume", "Geometry must be closed.");
      return NULL;
   }
   TGeoVolume *vol = (TGeoVolume*)map.GetValue(this);
   if (vol) {
      if (newname && newname[0]) vol->SetName(newname);
      return vol;
   }
//   printf("Making reflection for volume: %s\n", GetName());
   vol = CloneVolume();
   if (!vol) {
      Fatal("MakeReflectedVolume", "Cannot clone volume %s\n", GetName());
      return 0;
   }
   map.Add((TObject*)this, vol);
   if (newname && newname[0]) vol->SetName(newname);
   delete vol->GetNodes();
   vol->SetNodes(NULL);
   vol->SetBit(kVolumeImportNodes, kFALSE);
   CloneNodesAndConnect(vol);
   // The volume is now properly cloned, but with the same shape.
   // Reflect the shape (if any) and connect it.
   if (fShape) {
      TGeoShape *reflected_shape =
         TGeoScaledShape::MakeScaledShape(fShape->GetName(), fShape, new TGeoScale(1.,1.,-1.));
      vol->SetShape(reflected_shape);
   }
   // Reflect the daughters.
   Int_t nd = vol->GetNdaughters();
   if (!nd) return vol;
   TGeoNodeMatrix *node;
   TGeoMatrix *local, *local_cloned;
   TGeoVolume *new_vol;
   if (!vol->GetFinder()) {
      for (Int_t i=0; i<nd; i++) {
         node = (TGeoNodeMatrix*)vol->GetNode(i);
         local = node->GetMatrix();
//         printf("%s before\n", node->GetName());
//         local->Print();
         Bool_t reflected = local->IsReflection();
         local_cloned = new TGeoCombiTrans(*local);
         local_cloned->RegisterYourself();
         node->SetMatrix(local_cloned);
         if (!reflected) {
         // We need to reflect only the translation and propagate to daughters.
            // H' = Sz * H * Sz
            local_cloned->ReflectZ(kTRUE);
            local_cloned->ReflectZ(kFALSE);
//            printf("%s after\n", node->GetName());
//            node->GetMatrix()->Print();
            new_vol = node->GetVolume()->MakeReflectedVolume();
            node->SetVolume(new_vol);
            continue;
         }
         // The next daughter is already reflected, so reflect on Z everything and stop
         local_cloned->ReflectZ(kTRUE); // rot + tr
//         printf("%s already reflected... After:\n", node->GetName());
//         node->GetMatrix()->Print();
      }
      if (vol->GetVoxels()) vol->GetVoxels()->Voxelize();
      return vol;
   }
   // Volume is divided, so we have to reflect the division.
//   printf("   ... divided %s\n", fFinder->ClassName());
   TGeoPatternFinder *new_finder = fFinder->MakeCopy(kTRUE);
   if (!new_finder) {
      Fatal("MakeReflectedVolume", "Could not copy finder for volume %s", GetName());
      return 0;
   }
   new_finder->SetVolume(vol);
   vol->SetFinder(new_finder);
   TGeoNodeOffset *nodeoff;
   new_vol = 0;
   for (Int_t i=0; i<nd; i++) {
      nodeoff = (TGeoNodeOffset*)vol->GetNode(i);
      nodeoff->SetFinder(new_finder);
      new_vol = nodeoff->GetVolume()->MakeReflectedVolume();
      nodeoff->SetVolume(new_vol);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Set this volume as the TOP one (the whole geometry starts from here)

void TGeoVolume::SetAsTopVolume()
{
   fGeoManager->SetTopVolume(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current tracking point.

void TGeoVolume::SetCurrentPoint(Double_t x, Double_t y, Double_t z)
{
   fGeoManager->SetCurrentPoint(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////
/// set the shape associated with this volume

void TGeoVolume::SetShape(const TGeoShape *shape)
{
   if (!shape) {
      Error("SetShape", "No shape");
      return;
   }
   fShape = (TGeoShape*)shape;
}

////////////////////////////////////////////////////////////////////////////////
/// sort nodes by decreasing volume of the bounding box. ONLY nodes comes first,
/// then overlapping nodes and finally division nodes.

void TGeoVolume::SortNodes()
{
   if (!Valid()) {
      Error("SortNodes", "Bounding box not valid");
      return;
   }
   Int_t nd = GetNdaughters();
//   printf("volume : %s, nd=%i\n", GetName(), nd);
   if (!nd) return;
   if (fFinder) return;
//   printf("Nodes for %s\n", GetName());
   Int_t id = 0;
   TGeoNode *node = 0;
   TObjArray *nodes = new TObjArray(nd);
   Int_t inode = 0;
   // first put ONLY's
   for (id=0; id<nd; id++) {
      node = GetNode(id);
      if (node->InheritsFrom(TGeoNodeOffset::Class()) || node->IsOverlapping()) continue;
      nodes->Add(node);
//      printf("inode %i ONLY\n", inode);
      inode++;
   }
   // second put overlapping nodes
   for (id=0; id<nd; id++) {
      node = GetNode(id);
      if (node->InheritsFrom(TGeoNodeOffset::Class()) || (!node->IsOverlapping())) continue;
      nodes->Add(node);
//      printf("inode %i MANY\n", inode);
      inode++;
   }
   // third put the divided nodes
   if (fFinder) {
      fFinder->SetDivIndex(inode);
      for (id=0; id<nd; id++) {
         node = GetNode(id);
         if (!node->InheritsFrom(TGeoNodeOffset::Class())) continue;
         nodes->Add(node);
//         printf("inode %i DIV\n", inode);
         inode++;
      }
   }
   if (inode != nd) printf(" volume %s : number of nodes does not match!!!\n", GetName());
   delete fNodes;
   fNodes = nodes;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGeoVolume.

void TGeoVolume::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoVolume::Class(), this);
      if (fVoxels && fVoxels->IsInvalid()) Voxelize("");
   } else {
      if (!fVoxels) {
         R__b.WriteClassBuffer(TGeoVolume::Class(), this);
      } else {
         if (!fGeoManager->IsStreamingVoxels()) {
            TGeoVoxelFinder *voxels = fVoxels;
            fVoxels = 0;
            R__b.WriteClassBuffer(TGeoVolume::Class(), this);
            fVoxels = voxels;
         } else {
            R__b.WriteClassBuffer(TGeoVolume::Class(), this);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current options (none implemented)

void TGeoVolume::SetOption(const char *option)
{
   fOption = option;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line color.

void TGeoVolume::SetLineColor(Color_t lcolor)
{
   TAttLine::SetLineColor(lcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line style.

void TGeoVolume::SetLineStyle(Style_t lstyle)
{
   TAttLine::SetLineStyle(lstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line width.

void TGeoVolume::SetLineWidth(Style_t lwidth)
{
   TAttLine::SetLineWidth(lwidth);
}

////////////////////////////////////////////////////////////////////////////////
/// get the pointer to a daughter node

TGeoNode *TGeoVolume::GetNode(const char *name) const
{
   if (!fNodes) return 0;
   TGeoNode *node = (TGeoNode *)fNodes->FindObject(name);
   return node;
}

////////////////////////////////////////////////////////////////////////////////
/// get the total size in bytes for this volume

Int_t TGeoVolume::GetByteCount() const
{
   Int_t count = 28+2+6+4+0;    // TNamed+TGeoAtt+TAttLine+TAttFill+TAtt3D
   count += fName.Capacity() + fTitle.Capacity(); // name+title
   count += 7*sizeof(char*); // fShape + fMedium + fFinder + fField + fNodes + 2 extensions
   count += fOption.Capacity(); // fOption
   if (fShape)  count += fShape->GetByteCount();
   if (fFinder) count += fFinder->GetByteCount();
   if (fNodes) {
      count += 32 + 4*fNodes->GetEntries(); // TObjArray
      TIter next(fNodes);
      TGeoNode *node;
      while ((node=(TGeoNode*)next())) count += node->GetByteCount();
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// loop all nodes marked as overlaps and find overlapping brothers

void TGeoVolume::FindOverlaps() const
{
   if (!Valid()) {
      Error("FindOverlaps","Bounding box not valid");
      return;
   }
   if (!fVoxels) return;
   Int_t nd = GetNdaughters();
   if (!nd) return;
   TGeoNode *node=0;
   Int_t inode = 0;
   for (inode=0; inode<nd; inode++) {
      node = GetNode(inode);
      if (!node->IsOverlapping()) continue;
      fVoxels->FindOverlaps(inode);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an existing daughter.

void TGeoVolume::RemoveNode(TGeoNode *node)
{
   if (!fNodes || !fNodes->GetEntriesFast()) return;
   if (!fNodes->Remove(node)) return;
   fNodes->Compress();
   if (fVoxels) fVoxels->SetNeedRebuild();
   if (IsAssembly()) fShape->ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Replace an existing daughter with a new volume having the same name but
/// possibly a new shape, position or medium. Not allowed for positioned assemblies.
/// For division cells, the new shape/matrix are ignored.

TGeoNode *TGeoVolume::ReplaceNode(TGeoNode *nodeorig, TGeoShape *newshape, TGeoMatrix *newpos, TGeoMedium *newmed)
{
   Int_t ind = GetIndex(nodeorig);
   if (ind < 0) return NULL;
   TGeoVolume *oldvol = nodeorig->GetVolume();
   if (oldvol->IsAssembly()) {
      Error("ReplaceNode", "Cannot replace node %s since it is an assembly", nodeorig->GetName());
      return NULL;
   }
   TGeoShape  *shape = oldvol->GetShape();
   if (newshape && !nodeorig->IsOffset()) shape = newshape;
   TGeoMedium *med = oldvol->GetMedium();
   if (newmed) med = newmed;
   // Make a new volume
   TGeoVolume *vol = new TGeoVolume(oldvol->GetName(), shape, med);
   // copy volume attributes
   vol->SetVisibility(oldvol->IsVisible());
   vol->SetLineColor(oldvol->GetLineColor());
   vol->SetLineStyle(oldvol->GetLineStyle());
   vol->SetLineWidth(oldvol->GetLineWidth());
   vol->SetFillColor(oldvol->GetFillColor());
   vol->SetFillStyle(oldvol->GetFillStyle());
   // copy field
   vol->SetField(oldvol->GetField());
   // Make a copy of the node
   TGeoNode *newnode = nodeorig->MakeCopyNode();
   if (!newnode) {
      Fatal("ReplaceNode", "Cannot make copy node for %s", nodeorig->GetName());
      return 0;
   }
   // Change the volume for the new node
   newnode->SetVolume(vol);
   // Replace the matrix
   if (newpos && !nodeorig->IsOffset()) {
      TGeoNodeMatrix *nodemat = (TGeoNodeMatrix*)newnode;
      nodemat->SetMatrix(newpos);
   }
   // Replace nodeorig with new one
   fNodes->RemoveAt(ind);
   fNodes->AddAt(newnode, ind);
   if (fVoxels) fVoxels->SetNeedRebuild();
   if (IsAssembly()) fShape->ComputeBBox();
   return newnode;
}

////////////////////////////////////////////////////////////////////////////////
/// Select this volume as matching an arbitrary criteria. The volume is added to
/// a static list and the flag TGeoVolume::kVolumeSelected is set. All flags need
/// to be reset at the end by calling the method with CLEAR=true. This will also clear
/// the list.

void TGeoVolume::SelectVolume(Bool_t clear)
{
   static TObjArray array(256);
   static Int_t len = 0;
   Int_t i;
   TObject *vol;
   if (clear) {
      for (i=0; i<len; i++) {
         vol = array.At(i);
         vol->ResetBit(TGeoVolume::kVolumeSelected);
      }
      array.Clear();
      len = 0;
      return;
   }
   SetBit(TGeoVolume::kVolumeSelected);
   array.AddAtAndExpand(this, len++);
}

////////////////////////////////////////////////////////////////////////////////
/// set visibility of this volume

void TGeoVolume::SetVisibility(Bool_t vis)
{
   TGeoAtt::SetVisibility(vis);
   if (fGeoManager->IsClosed()) SetVisTouched(kTRUE);
   fGeoManager->SetVisOption(4);
   TSeqCollection *brlist = gROOT->GetListOfBrowsers();
   TIter next(brlist);
   TBrowser *browser = 0;
   while ((browser=(TBrowser*)next())) {
      browser->CheckObjectItem(this, vis);
      browser->Refresh();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility for containers.

void TGeoVolume::SetVisContainers(Bool_t flag)
{
   TGeoAtt::SetVisContainers(flag);
   if (fGeoManager && fGeoManager->IsClosed()) {
      if (flag) fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisDefault);
      else      fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisLeaves);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility for leaves.

void TGeoVolume::SetVisLeaves(Bool_t flag)
{
   TGeoAtt::SetVisLeaves(flag);
   if (fGeoManager && fGeoManager->IsClosed()) {
      if (flag) fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisLeaves);
      else      fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisDefault);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set visibility for leaves.

void TGeoVolume::SetVisOnly(Bool_t flag)
{
   if (IsAssembly()) return;
   TGeoAtt::SetVisOnly(flag);
   if (fGeoManager && fGeoManager->IsClosed()) {
      if (flag) fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisOnly);
      else      fGeoManager->SetVisOption(TVirtualGeoPainter::kGeoVisLeaves);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the shape of this volume is valid.

Bool_t TGeoVolume::Valid() const
{
   return fShape->IsValidBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Find a daughter node having VOL as volume and fill TGeoManager::fHMatrix
/// with its global matrix.

Bool_t TGeoVolume::FindMatrixOfDaughterVolume(TGeoVolume *vol) const
{
   if (vol == this) return kTRUE;
   Int_t nd = GetNdaughters();
   if (!nd) return kFALSE;
   TGeoHMatrix *global = fGeoManager->GetHMatrix();
   if (!global) return kFALSE;
   TGeoNode *dnode;
   TGeoVolume *dvol;
   TGeoMatrix *local;
   Int_t i;
   for (i=0; i<nd; i++) {
      dnode = GetNode(i);
      dvol = dnode->GetVolume();
      if (dvol == vol) {
         local = dnode->GetMatrix();
         global->MultiplyLeft(local);
         return kTRUE;
      }
   }
   for (i=0; i<nd; i++) {
      dnode = GetNode(i);
      dvol = dnode->GetVolume();
      if (dvol->FindMatrixOfDaughterVolume(vol)) return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// set visibility for daughters

void TGeoVolume::VisibleDaughters(Bool_t vis)
{
   SetVisDaughters(vis);
   if (fGeoManager->IsClosed()) SetVisTouched(kTRUE);
   fGeoManager->SetVisOption(4);
}

////////////////////////////////////////////////////////////////////////////////
/// build the voxels for this volume

void TGeoVolume::Voxelize(Option_t *option)
{
   if (!Valid()) {
      Error("Voxelize", "Bounding box not valid");
      return;
   }
   // do not voxelize divided volumes
   if (fFinder) return;
   // or final leaves
   Int_t nd = GetNdaughters();
   if (!nd) return;
   // If this is an assembly, re-compute bounding box
   if (IsAssembly()) fShape->ComputeBBox();
   // delete old voxelization if any
   if (fVoxels) {
      if (!TObject::TestBit(kVolumeClone)) delete fVoxels;
      fVoxels = 0;
   }
   // Create the voxels structure
   fVoxels = new TGeoVoxelFinder(this);
   fVoxels->Voxelize(option);
   if (fVoxels) {
      if (fVoxels->IsInvalid()) {
         delete fVoxels;
         fVoxels = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate the weight of a volume (in kg) with SIGMA(M)/M better than PRECISION.
/// Option can contain : v - verbose, a - analytical  (default)

Double_t TGeoVolume::Weight(Double_t precision, Option_t *option)
{
   TGeoVolume *top = fGeoManager->GetTopVolume();
   if (top != this) fGeoManager->SetTopVolume(this);
   else top = 0;
   Double_t weight =  fGeoManager->Weight(precision, option);
   if (top) fGeoManager->SetTopVolume(top);
   return weight;
}

////////////////////////////////////////////////////////////////////////////////
/// Analytical computation of the weight.

Double_t TGeoVolume::WeightA() const
{
   Double_t capacity = Capacity();
   Double_t weight = 0.0;
   Int_t i;
   Int_t nd = GetNdaughters();
   TGeoVolume *daughter;
   for (i=0; i<nd; i++) {
      daughter = GetNode(i)->GetVolume();
      weight += daughter->WeightA();
      capacity -= daughter->Capacity();
   }
   Double_t density = 0.0;
   if (!IsAssembly()) {
      if (fMedium) density = fMedium->GetMaterial()->GetDensity();
      if (density<0.01) density = 0.0; // do not weight gases
   }
   weight += 0.001*capacity * density; //[kg]
   return weight;
}

ClassImp(TGeoVolumeMulti);


////////////////////////////////////////////////////////////////////////////////
/// dummy constructor

TGeoVolumeMulti::TGeoVolumeMulti()
{
   fVolumes   = 0;
   fDivision = 0;
   fNumed = 0;
   fNdiv = 0;
   fAxis = 0;
   fStart = 0;
   fStep = 0;
   fAttSet = kFALSE;
   TObject::SetBit(kVolumeMulti);
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoVolumeMulti::TGeoVolumeMulti(const char *name, TGeoMedium *med)
{
   fVolumes = new TObjArray();
   fDivision = 0;
   fNumed = 0;
   fNdiv = 0;
   fAxis = 0;
   fStart = 0;
   fStep = 0;
   fAttSet = kFALSE;
   TObject::SetBit(kVolumeMulti);
   SetName(name);
   SetMedium(med);
   fGeoManager->AddVolume(this);
//   printf("--- volume multi %s created\n", name);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVolumeMulti::~TGeoVolumeMulti()
{
   if (fVolumes) delete fVolumes;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a volume with valid shape to the list of volumes. Copy all existing nodes
/// to this volume

void TGeoVolumeMulti::AddVolume(TGeoVolume *vol)
{
   Int_t idx = fVolumes->GetEntriesFast();
   fVolumes->AddAtAndExpand(vol,idx);
   vol->SetUniqueID(idx+1);
   TGeoVolumeMulti *div;
   TGeoVolume *cell;
   if (fDivision) {
      div = (TGeoVolumeMulti*)vol->Divide(fDivision->GetName(), fAxis, fNdiv, fStart, fStep, fNumed, fOption.Data());
      if (!div) {
         Fatal("AddVolume", "Cannot divide volume %s", vol->GetName());
         return;
      }
      for (Int_t i=0; i<div->GetNvolumes(); i++) {
         cell = div->GetVolume(i);
         fDivision->AddVolume(cell);
      }
   }
   if (fNodes) {
      Int_t nd = fNodes->GetEntriesFast();
      for (Int_t id=0; id<nd; id++) {
         TGeoNode *node = (TGeoNode*)fNodes->At(id);
         Bool_t many = node->IsOverlapping();
         if (many) vol->AddNodeOverlap(node->GetVolume(), node->GetNumber(), node->GetMatrix());
         else      vol->AddNode(node->GetVolume(), node->GetNumber(), node->GetMatrix());
      }
   }
//      vol->MakeCopyNodes(this);
}


////////////////////////////////////////////////////////////////////////////////
/// Add a new node to the list of nodes. This is the usual method for adding
/// daughters inside the container volume.

void TGeoVolumeMulti::AddNode(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t *option)
{
   TGeoVolume::AddNode(vol, copy_no, mat, option);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *volume = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      volume = GetVolume(ivo);
      volume->SetLineColor(GetLineColor());
      volume->SetLineStyle(GetLineStyle());
      volume->SetLineWidth(GetLineWidth());
      volume->SetVisibility(IsVisible());
      volume->AddNode(vol, copy_no, mat, option);
   }
//   printf("--- vmulti %s : node %s added to %i components\n", GetName(), vol->GetName(), nvolumes);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new node to the list of nodes, This node is possibly overlapping with other
/// daughters of the volume or extruding the volume.

void TGeoVolumeMulti::AddNodeOverlap(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t *option)
{
   TGeoVolume::AddNodeOverlap(vol, copy_no, mat, option);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *volume = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      volume = GetVolume(ivo);
      volume->SetLineColor(GetLineColor());
      volume->SetLineStyle(GetLineStyle());
      volume->SetLineWidth(GetLineWidth());
      volume->SetVisibility(IsVisible());
      volume->AddNodeOverlap(vol, copy_no, mat, option);
   }
//   printf("--- vmulti %s : node ovlp %s added to %i components\n", GetName(), vol->GetName(), nvolumes);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the last shape.

TGeoShape *TGeoVolumeMulti::GetLastShape() const
{
   TGeoVolume *vol = GetVolume(fVolumes->GetEntriesFast()-1);
   if (!vol) return 0;
   return vol->GetShape();
}

////////////////////////////////////////////////////////////////////////////////
/// division of multiple volumes

TGeoVolume *TGeoVolumeMulti::Divide(const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step, Int_t numed, const char *option)
{
   if (fDivision) {
      Error("Divide", "volume %s already divided", GetName());
      return 0;
   }
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoMedium *medium = fMedium;
   if (numed) {
      medium = fGeoManager->GetMedium(numed);
      if (!medium) {
         Error("Divide", "Invalid medium number %d for division volume %s", numed, divname);
         medium = fMedium;
      }
   }
   if (!nvolumes) {
      // this is a virtual volume
      fDivision = new TGeoVolumeMulti(divname, medium);
      fNumed = medium->GetId();
      fOption = option;
      fAxis = iaxis;
      fNdiv = ndiv;
      fStart = start;
      fStep = step;
      // nothing else to do at this stage
      return fDivision;
   }
   TGeoVolume *vol = 0;
   fDivision = new TGeoVolumeMulti(divname, medium);
   if (medium) fNumed = medium->GetId();
   fOption = option;
   fAxis = iaxis;
   fNdiv = ndiv;
   fStart = start;
   fStep = step;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetLineColor(GetLineColor());
      vol->SetLineStyle(GetLineStyle());
      vol->SetLineWidth(GetLineWidth());
      vol->SetVisibility(IsVisible());
      fDivision->AddVolume(vol->Divide(divname,iaxis,ndiv,start,step, numed, option));
   }
//   printf("--- volume multi %s (%i volumes) divided\n", GetName(), nvolumes);
   if (numed) fDivision->SetMedium(medium);
   return fDivision;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a copy of this volume
/// build a volume with same name, shape and medium

TGeoVolume *TGeoVolumeMulti::MakeCopyVolume(TGeoShape *newshape)
{
   TGeoVolume *vol = new TGeoVolume(GetName(), newshape, fMedium);
   Int_t i=0;
   // copy volume attributes
   vol->SetVisibility(IsVisible());
   vol->SetLineColor(GetLineColor());
   vol->SetLineStyle(GetLineStyle());
   vol->SetLineWidth(GetLineWidth());
   vol->SetFillColor(GetFillColor());
   vol->SetFillStyle(GetFillStyle());
   // copy field
   vol->SetField(fField);
   // Copy extensions
   vol->SetUserExtension(fUserExtension);
   vol->SetFWExtension(fFWExtension);
   // if divided, copy division object
//    if (fFinder) {
//       Error("MakeCopyVolume", "volume %s divided", GetName());
//       vol->SetFinder(fFinder);
//    }
   if (fDivision) {
      TGeoVolume *cell;
      TGeoVolumeMulti *div = (TGeoVolumeMulti*)vol->Divide(fDivision->GetName(), fAxis, fNdiv, fStart, fStep, fNumed, fOption.Data());
      if (!div) {
         Fatal("MakeCopyVolume", "Cannot divide volume %s", vol->GetName());
         return 0;
      }
      for (i=0; i<div->GetNvolumes(); i++) {
         cell = div->GetVolume(i);
         fDivision->AddVolume(cell);
      }
   }

   if (!fNodes) return vol;
   TGeoNode *node;
   Int_t nd = fNodes->GetEntriesFast();
   if (!nd) return vol;
   // create new list of nodes
   TObjArray *list = new TObjArray();
   // attach it to new volume
   vol->SetNodes(list);
   ((TObject*)vol)->SetBit(kVolumeImportNodes);
   for (i=0; i<nd; i++) {
      //create copies of nodes and add them to list
      node = GetNode(i)->MakeCopyNode();
      if (!node) {
         Fatal("MakeCopyNode", "cannot make copy node for daughter %d of %s", i, GetName());
         return 0;
      }
      node->SetMotherVolume(vol);
      list->Add(node);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line color for all components.

void TGeoVolumeMulti::SetLineColor(Color_t lcolor)
{
   TGeoVolume::SetLineColor(lcolor);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *vol = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetLineColor(lcolor);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line style for all components.

void TGeoVolumeMulti::SetLineStyle(Style_t lstyle)
{
   TGeoVolume::SetLineStyle(lstyle);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *vol = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetLineStyle(lstyle);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the line width for all components.

void TGeoVolumeMulti::SetLineWidth(Width_t lwidth)
{
   TGeoVolume::SetLineWidth(lwidth);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *vol = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetLineWidth(lwidth);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set medium for a multiple volume.

void TGeoVolumeMulti::SetMedium(TGeoMedium *med)
{
   TGeoVolume::SetMedium(med);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *vol = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetMedium(med);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set visibility for all components.

void TGeoVolumeMulti::SetVisibility(Bool_t vis)
{
   TGeoVolume::SetVisibility(vis);
   Int_t nvolumes = fVolumes->GetEntriesFast();
   TGeoVolume *vol = 0;
   for (Int_t ivo=0; ivo<nvolumes; ivo++) {
      vol = GetVolume(ivo);
      vol->SetVisibility(vis);
   }
}

ClassImp(TGeoVolumeAssembly);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoVolumeAssembly::ThreadData_t::ThreadData_t() :
   fCurrent(-1), fNext(-1)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVolumeAssembly::ThreadData_t::~ThreadData_t()
{
}

////////////////////////////////////////////////////////////////////////////////

TGeoVolumeAssembly::ThreadData_t& TGeoVolumeAssembly::GetThreadData() const
{
   Int_t tid = TGeoManager::ThreadId();
   return *fThreadData[tid];
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolumeAssembly::ClearThreadData() const
{
   std::lock_guard<std::mutex> guard(fMutex);
   TGeoVolume::ClearThreadData();
   std::vector<ThreadData_t*>::iterator i = fThreadData.begin();
   while (i != fThreadData.end())
   {
      delete *i;
      ++i;
   }
   fThreadData.clear();
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolumeAssembly::CreateThreadData(Int_t nthreads)
{
   std::lock_guard<std::mutex> guard(fMutex);
   // Create assembly thread data here
   fThreadData.resize(nthreads);
   fThreadSize = nthreads;
   for (Int_t tid=0; tid<nthreads; tid++) {
      if (fThreadData[tid] == 0) {
         fThreadData[tid] = new ThreadData_t;
      }
   }
   TGeoVolume::CreateThreadData(nthreads);
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGeoVolumeAssembly::GetCurrentNodeIndex() const
{
   return fThreadData[TGeoManager::ThreadId()]->fCurrent;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGeoVolumeAssembly::GetNextNodeIndex() const
{
   return fThreadData[TGeoManager::ThreadId()]->fNext;
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolumeAssembly::SetCurrentNodeIndex(Int_t index)
{
   fThreadData[TGeoManager::ThreadId()]->fCurrent = index;
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVolumeAssembly::SetNextNodeIndex(Int_t index)
{
   fThreadData[TGeoManager::ThreadId()]->fNext = index;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVolumeAssembly::TGeoVolumeAssembly()
                   :TGeoVolume()
{
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Just the name has to be provided. Assemblies does not have their own
/// shape or medium.

TGeoVolumeAssembly::TGeoVolumeAssembly(const char *name)
                   :TGeoVolume()
{
   fName = name;
   fName = fName.Strip();
   fShape = new TGeoShapeAssembly(this);
   if (fGeoManager) fNumber = fGeoManager->AddVolume(this);
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. The assembly is owner of its "shape".

TGeoVolumeAssembly::~TGeoVolumeAssembly()
{
   ClearThreadData();
   if (fShape) delete fShape;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a component to the assembly.

void TGeoVolumeAssembly::AddNode(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t *option)
{
   TGeoVolume::AddNode(vol,copy_no,mat,option);
//   ((TGeoShapeAssembly*)fShape)->RecomputeBoxLast();
   ((TGeoShapeAssembly*)fShape)->NeedsBBoxRecompute();
}

////////////////////////////////////////////////////////////////////////////////
/// Add an overlapping node - not allowed for assemblies.

void TGeoVolumeAssembly::AddNodeOverlap(TGeoVolume *vol, Int_t copy_no, TGeoMatrix *mat, Option_t *option)
{
   Warning("AddNodeOverlap", "Declaring assembly %s as possibly overlapping inside %s not allowed. Using AddNode instead !",vol->GetName(),GetName());
   AddNode(vol, copy_no, mat, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Clone this volume.
/// build a volume with same name, shape and medium

TGeoVolume *TGeoVolumeAssembly::CloneVolume() const
{
   TGeoVolume *vol = new TGeoVolumeAssembly(GetName());
   Int_t i;
   // copy other attributes
   Int_t nbits = 8*sizeof(UInt_t);
   for (i=0; i<nbits; i++)
      vol->SetAttBit(1<<i, TGeoAtt::TestAttBit(1<<i));
   for (i=14; i<24; i++)
      vol->SetBit(1<<i, TestBit(1<<i));

   // copy field
   vol->SetField(fField);
   // Set bits
   for (i=0; i<nbits; i++)
      vol->SetBit(1<<i, TObject::TestBit(1<<i));
   vol->SetBit(kVolumeClone);
   // make copy nodes
   vol->MakeCopyNodes(this);
//   CloneNodesAndConnect(vol);
   ((TGeoShapeAssembly*)vol->GetShape())->NeedsBBoxRecompute();
   // copy voxels
   TGeoVoxelFinder *voxels = 0;
   if (fVoxels) {
      voxels = new TGeoVoxelFinder(vol);
      vol->SetVoxelFinder(voxels);
   }
   // copy option, uid
   vol->SetOption(fOption);
   vol->SetNumber(fNumber);
   vol->SetNtotal(fNtotal);
   vol->SetTitle(GetTitle());
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Division makes no sense for assemblies.

TGeoVolume *TGeoVolumeAssembly::Divide(const char *, Int_t, Int_t, Double_t, Double_t, Int_t, Option_t *)
{
   Error("Divide","Assemblies cannot be divided");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign to the assembly a collection of identical volumes positioned according
/// a predefined pattern. The option can be spaced out or touching depending on the empty
/// space between volumes.

TGeoVolume *TGeoVolumeAssembly::Divide(TGeoVolume *cell, TGeoPatternFinder *pattern, Option_t *option)
{
   if (fNodes) {
      Error("Divide", "Cannot divide assembly %s since it has nodes", GetName());
      return NULL;
   }
   if (fFinder) {
      Error("Divide", "Assembly %s already divided", GetName());
      return NULL;
   }
   Int_t ncells = pattern->GetNdiv();
   if (!ncells || pattern->GetStep()<=0) {
      Error("Divide", "Pattern finder for dividing assembly %s not initialized. Use SetRange() method.", GetName());
      return NULL;
   }
   fFinder = pattern;
   TString opt(option);
   opt.ToLower();
   if (opt.Contains("spacedout")) fFinder->SetSpacedOut(kTRUE);
   else fFinder->SetSpacedOut(kFALSE);
   // Position volumes
   for (Int_t i=0; i<ncells; i++) {
      fFinder->cd(i);
      TGeoNodeOffset *node = new TGeoNodeOffset(cell, i, 0.);
      node->SetFinder(fFinder);
      fNodes->Add(node);
   }
   return cell;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of volume VOL but which is an assembly.

TGeoVolumeAssembly *TGeoVolumeAssembly::MakeAssemblyFromVolume(TGeoVolume *volorig)
{
   if (volorig->IsAssembly() || volorig->IsVolumeMulti()) return 0;
   Int_t nd = volorig->GetNdaughters();
   if (!nd) return 0;
   TGeoVolumeAssembly *vol = new TGeoVolumeAssembly(volorig->GetName());
   Int_t i;
   // copy other attributes
   Int_t nbits = 8*sizeof(UInt_t);
   for (i=0; i<nbits; i++)
      vol->SetAttBit(1<<i, volorig->TestAttBit(1<<i));
   for (i=14; i<24; i++)
      vol->SetBit(1<<i, volorig->TestBit(1<<i));

   // copy field
   vol->SetField(volorig->GetField());
   // Set bits
   for (i=0; i<nbits; i++)
      vol->SetBit(1<<i, volorig->TestBit(1<<i));
   vol->SetBit(kVolumeClone);
   // make copy nodes
   vol->MakeCopyNodes(volorig);
//   volorig->CloneNodesAndConnect(vol);
   vol->GetShape()->ComputeBBox();
   // copy voxels
   TGeoVoxelFinder *voxels = 0;
   if (volorig->GetVoxels()) {
      voxels = new TGeoVoxelFinder(vol);
      vol->SetVoxelFinder(voxels);
   }
   // copy option, uid
   vol->SetOption(volorig->GetOption());
   vol->SetNumber(volorig->GetNumber());
   vol->SetNtotal(volorig->GetNtotal());
   return vol;
}
